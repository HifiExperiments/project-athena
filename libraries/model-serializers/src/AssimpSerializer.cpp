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
#include <assimp/postprocess.h>
#include <assimp/IOStream.hpp>
#include <assimp/IOSystem.hpp>
#include <assimp/pbrmaterial.h>

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

HFMModel::Pointer AssimpSerializer::read(const hifi::ByteArray& data, const hifi::VariantHash& mapping, const hifi::URL& inputUrl) {
	_url = inputUrl;

    hifi::URL normalizeUrl = DependencyManager::get<ResourceManager>()->normalizeURL(_url);
    if (normalizeUrl.scheme().isEmpty() || (normalizeUrl.scheme() == "file")) {
        QString localFileName = PathUtils::expandToLocalDataAbsolutePath(normalizeUrl).toLocalFile();
        _url = hifi::URL(QFileInfo(localFileName).absoluteFilePath());
    }

    _extension = QFileInfo(_url.path()).completeSuffix();

    qCDebug(modelformat) << "AssimpSerializer::read url"<< _url << _extension;

    Assimp::Importer importer;
    //importer.SetIOHandler(new TivoliIOSystem(url));

    //scene = importer.ReadFileFromMemory(
    //    data.constData(), data.size(),
    //    aiProcess_CalcTangentSpace |
    //    aiProcess_JoinIdenticalVertices | 
    //    aiProcess_Triangulate |
    //    aiProcess_GenNormals | 
    //    // aiProcess_SplitLargeMeshes |
    //    aiProcess_RemoveRedundantMaterials |
    //    // aiProcess_GenUVCoords | 
    //    // aiProcess_OptimizeMeshes |
    //    // aiProcess_OptimizeGraph,
    //    aiProcess_FlipUVs,
    //    _extension.toLocal8Bit()
    //);

    //if (!scene) {
    //    qCDebug(modelformat) << "AssimpSerializer::read Error parsing model file" << importer.GetErrorString();
    //    return nullptr;
    //}

    //processScene();

    // hfmModel->debugDump();

    return nullptr;//hfmModel;
}
