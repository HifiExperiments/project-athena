//
//  FBXToJSON.h
//  libraries/fbx/src
//
//  Created by Simon Walton on 5/4/2013.
//  Copyright 2018 High Fidelity, Inc.
//
//  Distributed under the Apache License, Version 2.0.
//  See the accompanying file LICENSE or http://www.apache.org/licenses/LICENSE-2.0.html
//

#ifndef hifi_FBXToJSON_h
#define hifi_FBXToJSON_h

#include <sstream>

// Forward declarations.
class FBXNode;
template<typename T> class QVector;

class FBXToJSON : public std::ostringstream {
public:
    FBXToJSON& operator<<(const FBXNode& fbxNode);

private:
    template<typename T> FBXToJSON& operator<<(QVector<T>& arrayProp);
    int _indentLevel { 0 };
};

#endif  // hifi_FBXToJSON_h
