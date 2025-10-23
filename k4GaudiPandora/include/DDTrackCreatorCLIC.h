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
 *  @file   DDMarlinPandora/include/DDTrackCreatorCLIC.h
 *
 *  @brief  Header file for the ILD implementation of the track creator class.
 *
 *  $Log: $
 */

#ifndef DDTRACK_CREATOR_CLIC_H
#define DDTRACK_CREATOR_CLIC_H 1

#include "DDTrackCreatorBase.h"

#include <Gaudi/Algorithm.h>

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDTrackCreatorCLIC class
 */
class DDTrackCreatorCLIC : public DDTrackCreatorBase {
public:
  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pPandora address of the relevant pandora instance
   */
  DDTrackCreatorCLIC(const Gaudi::Algorithm*, const Settings& settings, const pandora::Pandora* const pPandora);

  /**
   *  @brief  Create tracks implementation, insert user code here. Detector specific
   *
   *  @param  pLCEvent the lcio event
   */
  pandora::StatusCode CreateTracks(const std::vector<edm4hep::Track>& tracks) override;

protected:
  // Detector-specific configuration variables
  float m_trackerInnerR; ///< Inner radius of the barrel tracker
  float m_trackerOuterR; ///< Outer radius of the barrel tracker
  float m_trackerZmax;   ///< max extent of the tracker along z
  float m_cosTracker;

  DoubleVector m_endcapDiskInnerRadii; ///< List of endcapDisk inner radii
  DoubleVector m_endcapDiskOuterRadii; ///< List of endcapDisk outer radii
  DoubleVector m_endcapDiskZPositions; ///< List of endcapDisk z positions
  unsigned int m_nEndcapDiskLayers;    ///< Number of endcapDisk layers
  unsigned int m_barrelTrackerLayers;  ///< Number of barrel tracker layers

  float m_tanLambdaEndcapDisk; ///< Tan lambda for first endcapDisk layer

  /**
   *  @brief  Whether track passes the quality cuts required in order to be used to form a pfo. Detector specific
   *
   *  @param  pTrack the lcio track
   *  @param  trackParameters the track parameters
   *
   *  @return boolean
   */

  bool PassesQualityCuts(const edm4hep::Track& pTrack,
                         const PandoraApi::Track::Parameters& trackParameters) const override;

  /**
   *  @brief  Decide whether track reaches the ecal surface. Detector specific
   *
   *  @param  pTrack the lcio track
   *  @param  trackParameters the track parameters
   */
  void TrackReachesECAL(const edm4hep::Track& pTrack, PandoraApi::Track::Parameters& trackParameters) const override;

  /**
   *  @brief  Determine whether a track can be used to form a pfo under the following conditions:
   *          1) if the track proves to be associated with a cluster, OR
   *          2) if the track proves to have no cluster associations
   *          Detector specific
   *
   *  @param  pTrack the lcio track
   *  @param  trackParameters the track parameters
   */
  void DefineTrackPfoUsage(const edm4hep::Track& pTrack, PandoraApi::Track::Parameters& trackParameters) const override;
};

#endif // #ifndef DDTRACK_CREATOR_CLIC_H
