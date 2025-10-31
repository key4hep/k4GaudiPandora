/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 *  @file   DDMarlinPandora/include/DDBFieldPlugin.h
 *
 *  @brief  Header file for the BField plugin class.
 *
 *  $Log: $
 */

#ifndef K4GAUDIPANDORA_DDBFIELD_PLUGIN_H
#define K4GAUDIPANDORA_DDBFIELD_PLUGIN_H 1

#include "Objects/CartesianVector.h"
#include "Plugins/BFieldPlugin.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDBFieldPlugin class
 */
class DDBFieldPlugin : public pandora::BFieldPlugin {
public:
  DDBFieldPlugin(const dd4hep::Detector& detector);
  float GetBField(const pandora::CartesianVector& positionVector) const override;

private:
  pandora::StatusCode Initialize() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  dd4hep::OverlayedField m_field; ///< The field instance from DD4hep
};

#endif // #ifndef K4GAUDIPANDORA_DDBFIELD_PLUGIN_H
