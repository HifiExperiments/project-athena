//
//  AssimpSerializer.cpp
//  libraries/model-serializers/src
//
//  Created by HifiExperiments on 4/16/21
//  Copyright 2021 Vircadia contributors
//
//  Distributed under the Apache License, Version 2.0.
//  See the accompanying file LICENSE or http://www.apache.org/licenses/LICENSE-2.0.html
//

#include "AssimpSerializer.h"

#include <assimp/Importer.hpp>
#include <assimp/pbrmaterial.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <QThreadStorage>
#include <QImage>

#include <ResourceCache.h>
#include <PathUtils.h>
#include <hfm/ModelFormatLogging.h>

MediaType AssimpSerializer::getMediaType() const {
    MediaType mediaType;

    Assimp::Importer importer;
    aiString aiExts;
    importer.GetExtensionList(aiExts);

    QStringList exts = QString(aiExts.C_Str()).split(";");
    foreach (QString ext, exts) {
        mediaType.extensions.push_back(ext.remove("*.").toStdString());
    }
    return mediaType;
}

std::unique_ptr<hfm::Serializer::Factory> AssimpSerializer::getFactory() const {
    return std::make_unique<hfm::Serializer::SimpleFactory<AssimpSerializer>>();
}

HFMTexture AssimpSerializer::getTexture(const aiString& aiURL, const aiScene* scene, const hifi::URL& modelURL) {
    std::string url = aiURL.C_Str();

    auto itr = textures.find(url);
    if (itr != textures.end()) {
        return itr->second;
    } else {
        HFMTexture texture = HFMTexture();

        QRegExp indexFinder = QRegExp("^\\*([0-9]+)$");
        if (indexFinder.exactMatch(QString(url.c_str()))) {
            auto index = indexFinder.cap(1).toUInt();

            if (index < scene->mNumTextures) {
                auto aiTexture = scene->mTextures[index];
                auto width = aiTexture->mWidth;
                auto height = aiTexture->mHeight;

                texture.name = aiTexture->mFilename.C_Str();
                texture.filename = QString::number(index).toLocal8Bit();

                if (height == 0) {
                    texture.content = QByteArray((char*)aiTexture->pcData, width);
                } else {
                    texture.content = QByteArray((char*)aiTexture->pcData, width * height);
                }
            }
        } else {
            QUrl textureURL = modelURL.resolved(QString(url.c_str()));
            texture.name = textureURL.fileName();
            texture.filename = textureURL.toEncoded();
        }

        textures[url] = texture;
        return texture;
    }
}

void AssimpSerializer::parseMaterials(const aiScene* scene, const HFMModel::Pointer& hfmModel, const hifi::URL& modelURL) {
    for(unsigned int materialIndex = 0; materialIndex < scene->mNumMaterials; materialIndex++) {
        aiMaterial* mMaterial = scene->mMaterials[materialIndex];

        hfm::Material material;

        material.name = material.materialID = mMaterial->GetName().C_Str();
        material._material = std::make_shared<graphics::Material>();

        // Normal materials

        aiColor3D diffuseColor;
        if (mMaterial->Get(AI_MATKEY_COLOR_DIFFUSE, diffuseColor) == AI_SUCCESS) {
            material.diffuseColor = glm::vec3(diffuseColor.r, diffuseColor.g, diffuseColor.b);
            material._material->setAlbedo(material.diffuseColor);
        }

        aiColor3D specularColor;
        if (mMaterial->Get(AI_MATKEY_COLOR_SPECULAR, specularColor) == AI_SUCCESS) {
            material.specularColor = glm::vec3(specularColor.r, specularColor.g, specularColor.b);
            float metallic = std::max(material.specularColor.r, std::max(material.specularColor.g, material.specularColor.b));
            material._material->setMetallic(metallic);
        }

        aiColor3D emissiveColor;
        if (mMaterial->Get(AI_MATKEY_COLOR_EMISSIVE, emissiveColor) == AI_SUCCESS) {
            material.emissiveColor = glm::vec3(emissiveColor.r, emissiveColor.g, emissiveColor.b);
            material._material->setEmissive(material.emissiveColor);
        }

        int cullNone;
        if (mMaterial->Get(AI_MATKEY_TWOSIDED, cullNone) == AI_SUCCESS && cullNone) {
            material._material->setCullFaceMode(graphics::MaterialKey::CullFaceMode::CULL_NONE);
        }

        float opacity;
        if (mMaterial->Get(AI_MATKEY_OPACITY, opacity) == AI_SUCCESS) {
            material.opacity = opacity;
            material._material->setOpacity(opacity);
        }

        float shininess;
        if (mMaterial->Get(AI_MATKEY_SHININESS, shininess) == AI_SUCCESS) {
            material.shininess = shininess;
            material._material->setRoughness(graphics::Material::shininessToRoughness(material.shininess));
        }

        aiString diffuseTexture;
        if (mMaterial->GetTexture(aiTextureType_DIFFUSE, 0, &diffuseTexture) == AI_SUCCESS) {
            material.albedoTexture = getTexture(diffuseTexture, scene, modelURL);
            material.opacityTexture = getTexture(diffuseTexture, scene, modelURL);
            material.useAlbedoMap = true;
        }

        aiString emissiveTexture;
        if (mMaterial->GetTexture(aiTextureType_EMISSIVE, 0, &emissiveTexture) == AI_SUCCESS) {
            material.emissiveTexture = getTexture(emissiveTexture, scene, modelURL);
            material.useEmissiveMap = true;
        }

        aiString normalTexture;
        if (mMaterial->GetTexture(aiTextureType_NORMALS, 0, &normalTexture) == AI_SUCCESS) {
            material.normalTexture = getTexture(normalTexture, scene, modelURL);
            material.useNormalMap = true;
        }

        aiString lightmapTexture;
        if (mMaterial->GetTexture(aiTextureType_LIGHTMAP, 0, &lightmapTexture) == AI_SUCCESS) {
            material.lightmapTexture = getTexture(lightmapTexture, scene, modelURL);
        }

        // GLTF/PBR materials

        aiColor4D baseColorFactor;
        if (mMaterial->Get(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_BASE_COLOR_FACTOR, baseColorFactor) == AI_SUCCESS) {
            material.diffuseColor = glm::vec3(baseColorFactor.r, baseColorFactor.g, baseColorFactor.b);
            material.opacity = baseColorFactor.a;
            material._material->setAlbedo(material.diffuseColor);
            material._material->setOpacity(material.opacity);
        }

        float metallicFactor;
        if (mMaterial->Get(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_METALLIC_FACTOR, metallicFactor) == AI_SUCCESS) {
            material.metallic = metallicFactor;
            material._material->setMetallic(metallicFactor);
        }

        float roughnessFactor;
        if (mMaterial->Get(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_ROUGHNESS_FACTOR, roughnessFactor) == AI_SUCCESS) {
            material.roughness = roughnessFactor;
            material._material->setRoughness(roughnessFactor);
        }

        aiString alphaMode;
        if (mMaterial->Get(AI_MATKEY_GLTF_ALPHAMODE, alphaMode) == AI_SUCCESS) {
            QString alphaModeString = alphaMode.C_Str();
            if (alphaModeString == "BLEND") {
                material.alphaMode = graphics::MaterialKey::OPACITY_MAP_BLEND;
                 material._material->setOpacityMapMode(graphics::MaterialKey::OPACITY_MAP_BLEND);
            } else if (alphaModeString == "MASK") {
                material.alphaMode = graphics::MaterialKey::OPACITY_MAP_MASK;
                 material._material->setOpacityMapMode(graphics::MaterialKey::OPACITY_MAP_MASK);
            } else {
                material.alphaMode = graphics::MaterialKey::OPACITY_MAP_OPAQUE;
                 material._material->setOpacityMapMode(graphics::MaterialKey::OPACITY_MAP_OPAQUE);
            }
        }

        float alphaCutoff;
        if (mMaterial->Get(AI_MATKEY_GLTF_ALPHACUTOFF, alphaCutoff) == AI_SUCCESS) {
            material.alphaCutoff = alphaCutoff;
            material._material->setOpacityCutoff(alphaCutoff);
        }

        int unlit;
        if (mMaterial->Get(AI_MATKEY_GLTF_UNLIT, unlit) == AI_SUCCESS) {
            material._material->setUnlit(unlit);
        }

        aiString baseColorTexture;
        if (mMaterial->GetTexture(aiTextureType_BASE_COLOR, 0, &baseColorTexture) == AI_SUCCESS) {
            material.albedoTexture = getTexture(baseColorTexture, scene, modelURL);
            material.opacityTexture = material.albedoTexture;
            material.useAlbedoMap = true;
        }

        aiString metallicTexture;
        if (mMaterial->GetTexture(aiTextureType_METALNESS, 0, &metallicTexture) == AI_SUCCESS) {
            material.metallicTexture = getTexture(metallicTexture, scene, modelURL);
            material.useMetallicMap = true;
        }

        aiString diffuseRoughnessTexture;
        if (mMaterial->GetTexture(aiTextureType_DIFFUSE_ROUGHNESS, 0, &diffuseRoughnessTexture) == AI_SUCCESS) {
            material.albedoTexture = getTexture(diffuseRoughnessTexture, scene, modelURL);
            material.roughnessTexture = getTexture(diffuseRoughnessTexture, scene, modelURL);
            material.roughnessTexture.sourceChannel = image::ColorChannel::ALPHA;
            material.useAlbedoMap = true;
        }

        aiString occlusionTexture;
        if (mMaterial->GetTexture(aiTextureType_AMBIENT_OCCLUSION, 0, &occlusionTexture) == AI_SUCCESS) {
            material.occlusionTexture = getTexture(occlusionTexture, scene, modelURL);
            material.useOcclusionMap = true;
        }

        aiString metallicRoughnessTexture;
        if (mMaterial->GetTexture(AI_MATKEY_GLTF_PBRMETALLICROUGHNESS_METALLICROUGHNESS_TEXTURE, &metallicRoughnessTexture) == AI_SUCCESS) {          
            material.roughnessTexture = getTexture(metallicRoughnessTexture, scene, modelURL);
            material.roughnessTexture.sourceChannel = image::ColorChannel::GREEN;
            material.useRoughnessMap = true;

            material.metallicTexture = material.roughnessTexture;
            material.metallicTexture.sourceChannel = image::ColorChannel::BLUE;
            material.useMetallicMap = true;
        }

        hfmModel->materials.insert(material.name, material);
    }
}


void AssimpSerializer::parseMeshes(const aiScene* scene, const HFMModel::Pointer& hfmModel) {
    for (size_t meshIndex = 0; meshIndex < scene->mNumMeshes; meshIndex++) {
        aiMesh* mMesh = scene->mMeshes[meshIndex];

        hfm::Mesh mesh;
        mesh.meshIndex = (unsigned int)meshIndex;
        hfm::MeshPart meshPart;

        for (size_t vertexIndex = 0; vertexIndex < mMesh->mNumVertices; vertexIndex++) {
            aiVector3D& v = mMesh->mVertices[vertexIndex];
            mesh.vertices.push_back(glm::vec3(v.x, v.y, v.z));

            if (mMesh->HasNormals()) {
                aiVector3D& n = mMesh->mNormals[vertexIndex];
                mesh.normals.push_back(glm::vec3(n.x, n.y, n.z));
            }

            if (mMesh->HasTangentsAndBitangents()) {
                aiVector3D& t = mMesh->mTangents[vertexIndex];
                mesh.tangents.push_back(glm::vec3(t.x, t.y, t.z));
            }

            if (mMesh->HasVertexColors(0)) {
                aiColor4D& c = mMesh->mColors[0][vertexIndex];
                mesh.colors.push_back(glm::vec3(c.r, c.g, c.b));
            }

            if (mMesh->HasTextureCoords(0)) {
                aiVector3D& uv = mMesh->mTextureCoords[0][vertexIndex];
                mesh.texCoords.push_back(glm::vec2(uv.x, uv.y));
            }

            if (mMesh->HasTextureCoords(1)) {
                aiVector3D& uv = mMesh->mTextureCoords[1][vertexIndex];
                mesh.texCoords1.push_back(glm::vec2(uv.x, uv.y));
            }
        }

        for (size_t faceIndex = 0; faceIndex < mMesh->mNumFaces; faceIndex++) {
            aiFace& face = mMesh->mFaces[faceIndex];
            for (size_t i = 0; i < face.mNumIndices; i++) {
                if (face.mNumIndices == 3) {
                    meshPart.triangleIndices.push_back(face.mIndices[i]);
                } else if (face.mNumIndices == 4) {
                    // We shouldn't have any quads because we have Assimp turn them into triangles,
                    // but just in case...
                    meshPart.quadIndices.push_back(face.mIndices[i]);
                }
            }
        }

        //hifi::VariantHash blendshapeMappings = mapping.value("bs").toHash();

        //if (mMesh->mNumAnimMeshes > 0) {
        //    bool isBlendshapeChannelNamesEmpty = hfmModel->blendshapeChannelNames.size() == 0;

        //    for (size_t i = 0; i < (int)Blendshapes::BlendshapeCount; i++) {
        //        const QString fromName = BLENDSHAPE_NAMES[i];
        //        
        //        if (isBlendshapeChannelNamesEmpty) {
        //            hfmModel->blendshapeChannelNames.push_back(fromName);
        //        }

        //        meshPtr->blendshapes.push_back(hfm::Blendshape());
        //        hfm::Blendshape* blendshapePtr = &meshPtr->blendshapes.back();

        //        QString toName;
        //        float weight = 1.0f;

        //        if (blendshapeMappings.contains(fromName)) {
        //            auto mapping = blendshapeMappings[fromName].toList();
        //            if (mapping.count() >= 2) {
        //                toName = mapping.at(0).toString();
        //                weight = mapping.at(1).toFloat();
        //            }
        //        } else {
        //            toName = fromName;
        //        }

        //        if (weight == 0.0f && fromName == toName) {
        //            continue;
        //        }

        //        aiAnimMesh* animMesh = getBlendshapeByName(toName, mMesh);
        //        if (animMesh == nullptr) {
        //            continue;
        //        }
        //        
        //        for (size_t vertexIndex = 0; vertexIndex < animMesh->mNumVertices; vertexIndex++) {
        //            auto v = (asVec3(animMesh->mVertices[vertexIndex]) - asVec3(mMesh->mVertices[vertexIndex])) * weight;

        //            auto n = animMesh->HasNormals() ?
        //                (asVec3(animMesh->mNormals[vertexIndex]) - asVec3(mMesh->mNormals[vertexIndex])) * weight :
        //                vec3(0.0f);

        //            auto t = animMesh->HasTangentsAndBitangents() ?
        //                (asVec3(animMesh->mTangents[vertexIndex]) - asVec3(mMesh->mTangents[vertexIndex])) * weight :
        //                vec3(0.0f);
        //            
        //            if (isEmptyVec3(v) && isEmptyVec3(n) && isEmptyVec3(t)) {
        //                continue;
        //            }

        //            blendshapePtr->indices.push_back(vertexIndex);
        //            
        //            blendshapePtr->vertices.push_back(v);
        //            if (animMesh->HasNormals()) blendshapePtr->normals.push_back(n);
        //            if (animMesh->HasTangentsAndBitangents()) blendshapePtr->tangents.push_back(t);
        //        }
        //    }
        //}

        mesh.parts.push_back(meshPart);
        hfmModel->meshes.push_back(mesh);
    }
}

void AssimpSerializer::parseNode(const HFMModel::Pointer& hfmModel, aiNode *node, int parentIndex) {
    auto nodeIndex = hfmModel->joints.size();

    hfm::Joint joint;

    joint.name = node->mName.C_Str();
    hfmModel->jointIndices[joint.name] = nodeIndex;

    joint.parentIndex = parentIndex;
    joint.isSkeletonJoint = false;

    joint.transform = glm::mat4(node->mTransformation.a1, node->mTransformation.b1, node->mTransformation.c1, node->mTransformation.d1,
                                node->mTransformation.a2, node->mTransformation.b2, node->mTransformation.c2, node->mTransformation.d2,
                                node->mTransformation.a3, node->mTransformation.b3, node->mTransformation.c3, node->mTransformation.d3,
                                node->mTransformation.a4, node->mTransformation.b4, node->mTransformation.c4, node->mTransformation.d4);
    joint.translation = extractTranslation(joint.transform);
    joint.rotation = glmExtractRotation(joint.transform);
    glm::vec3 scale = extractScale(joint.transform);
    joint.postTransform = glm::scale(glm::mat4(), scale);

    if (joint.parentIndex == -1) {
        joint.transform = hfmModel->offset * joint.transform;
    } else {
        const auto& parentNode = hfmModel->joints[joint.parentIndex];
        joint.transform = parentNode.transform * joint.transform;
    }

    for (size_t i = 0; i < node->mNumMeshes; i++) {
        auto meshIndex = node->mMeshes[i];
        hfm::Mesh& mesh = hfmModel->meshes[meshIndex];

    //    hfm::Shape shape;

    //    // no verticies and indices causes segfault
    //    if (mesh->vertices.size() == 0) continue;
    //    bool hasIndicies = false;
    //    for (auto part : mesh->parts) {
    //        if (part.triangleIndices.size() > 0) {
    //            hasIndicies = true;
    //            break;
    //        }
    //    }
    //    if (!hasIndicies) continue;

    //    hfmModel->shapes.emplace_back();
    //    hfm::Shape* shapePtr = &hfmModel->shapes.back();

    //    shapePtr->joint = nodeIndex;
    //    shapePtr->mesh = meshIndex;
    //    shapePtr->meshPart = 0;

    //    // TODO: fixes bounding box? seems not for https://files.tivolicloud.com/caitlyn/fun/cu-cat/cu-cat.glb
    //    // mesh->modelTransform = joint.globalTransform;
    //    hfmModel->meshIndicesToModelNames[meshIndex] = joint.name;

    //    for (auto vertex : mesh->vertices) {
    //        glm::vec3 transformedVertex = glm::vec3(joint.globalTransform * glm::vec4(vertex, 1.0f));
    //        mesh->meshExtents.addPoint(transformedVertex);
    //    }
    //    shapePtr->transformedExtents = mesh->meshExtents;
    //    hfmModel->meshExtents.addExtents(mesh->meshExtents);

    //    auto materialIndex = scene->mMeshes[meshIndex]->mMaterialIndex;
    //    if (materialIndex < hfmModel->materials.size()) {
    //        shapePtr->material = materialIndex;
    //    }
    }

    hfmModel->joints.push_back(joint);

    for (size_t i = 0; i < node->mNumChildren; i++) {
        parseNode(hfmModel, node->mChildren[i], nodeIndex);
    }
}

HFMModel::Pointer AssimpSerializer::read(const hifi::ByteArray& data, const hifi::VariantHash& mapping, const hifi::URL& inputUrl) {
	hifi::URL url = inputUrl;

    hifi::URL normalizeUrl = DependencyManager::get<ResourceManager>()->normalizeURL(url);
    if (normalizeUrl.scheme().isEmpty() || (normalizeUrl.scheme() == "file")) {
        QString localFileName = PathUtils::expandToLocalDataAbsolutePath(normalizeUrl).toLocalFile();
        url = hifi::URL(QFileInfo(localFileName).absoluteFilePath());
    }

    QString extension = QFileInfo(url.path()).completeSuffix();

    qCDebug(modelformat) << "AssimpSerializer::read" << url;

    //static QThreadStorage<Assimp::Importer> importers;
    ////importer.SetIOHandler(new TivoliIOSystem(url));

    //if (!importers.hasLocalData()) {
    //    importers.setLocalData(Assimp::Importer());
    //}

    auto &importer = Assimp::Importer();//importers.localData();

    const aiScene *scene = importer.ReadFileFromMemory(
        data.constData(), data.size(),
        aiProcess_CalcTangentSpace |
        aiProcess_JoinIdenticalVertices |
        aiProcess_Triangulate |
        aiProcess_GenNormals |
        aiProcess_SplitLargeMeshes |
        aiProcess_LimitBoneWeights |
        aiProcess_ImproveCacheLocality |
        aiProcess_PopulateArmatureData |
        aiProcess_GenUVCoords |
        aiProcess_FlipUVs |
        aiProcess_OptimizeMeshes,
        extension.toLocal8Bit()
    );

    if (!scene) {
        qCDebug(modelformat) << "AssimpSerializer::read Error parsing model file" << importer.GetErrorString();
        return nullptr;
    }

    HFMModel::Pointer hfmModel = std::make_shared<HFMModel>();
    hfmModel->originalURL = url.toString();

    parseMaterials(scene, hfmModel, url);
    parseMeshes(scene, hfmModel);

    hfmModel->jointIndices["x"] = 0;
    parseNode(hfmModel, scene->mRootNode);

    importer.FreeScene();

    //hfmModel->debugDump();
    return hfmModel;
}
