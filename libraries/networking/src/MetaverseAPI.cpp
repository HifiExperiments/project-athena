//
//  MetaverseAPI.cpp
//  libraries/networking/src
//
//  Created by Kalila (kasenvr) on 2019-12-16.
//  Copyright 2019 High Fidelity, Inc.
//
//  Distributed under the Apache License, Version 2.0.
//  See the accompanying file LICENSE or http://www.apache.org/licenses/LICENSE-2.0.html
//

#include "MetaverseAPI.h"

#include <QUrl>
#include <QDebug>
#include "NetworkingConstants.h"
#include "SettingHandle.h"
#include < time.h>


namespace MetaverseAPI {
    time_t start = time(0);
    double seconds_since_start = difftime(time(0), start);

    // You can change the return of this function if you want to use a custom metaverse URL at compile time
    // or you can pass a custom URL via the env variable
    QUrl getCurrentMetaverseServerURL() {
        QUrl selectedMetaverseURL;
        //QVariant tempURLSettingString;

        Setting::Handle<QUrl> selectedMetaverseURLSetting { "private/selectedMetaverseURL",
                                                       NetworkingConstants::METAVERSE_SERVER_URL_STABLE };

        selectedMetaverseURL = QUrl(selectedMetaverseURLSetting.get());
        // QVariant tempURLSetting = Setting::Handle<QVariant>("private/selectedMetaverseURL").get();

        qDebug() << selectedMetaverseURL << "LOOK HERE!";
        const QString HIFI_METAVERSE_URL_ENV = "HIFI_METAVERSE_URL";
        QUrl serverURL;
        if (QProcessEnvironment::systemEnvironment().contains(HIFI_METAVERSE_URL_ENV)) {
            serverURL = QUrl(QProcessEnvironment::systemEnvironment().value(HIFI_METAVERSE_URL_ENV));

        } else if (!selectedMetaverseURLSetting.isSet()) {
            serverURL = selectedMetaverseURL;

        } else {
            serverURL = NetworkingConstants::METAVERSE_SERVER_URL_STABLE;

        }

        return serverURL;
    };
}
