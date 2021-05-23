//
//  AssimpSerializer.h
//  libraries/model-serializers/src
//
//  Created by HifiExperiments on 4/16/21.
//  Copyright 2021 Vircadia contributors
//
//  Distributed under the Apache License, Version 2.0.
//  See the accompanying file LICENSE or http://www.apache.org/licenses/LICENSE-2.0.html
//

#ifndef hifi_AssimpSerializer_h
#define hifi_AssimpSerializer_h

#include <hfm/HFMSerializer.h>

struct aiNode;
struct aiScene;
struct aiString;

class AssimpSerializer : public QObject, public HFMSerializer {
    Q_OBJECT

public:
    MediaType getMediaType() const override;
    std::unique_ptr<hfm::Serializer::Factory> getFactory() const override;

    HFMModel::Pointer read(const hifi::ByteArray& data, const hifi::VariantHash& mapping, const hifi::URL& url = hifi::URL()) override;

private:
    std::unordered_map<std::string, HFMTexture> textures;

    HFMTexture getTexture(const aiString& aiURL, const aiScene *scene, const hifi::URL& modelURL);
    void parseMaterials(const aiScene *scene, const HFMModel::Pointer& hfmModel, const hifi::URL& modelURL);
    void parseMeshes(const aiScene* scene, const HFMModel::Pointer& hfmModel);
    void parseNode(const HFMModel::Pointer& hfmModel, aiNode *node, int parentIndex = -1);
};

#endif // hifi_AssimpSerializer_h
