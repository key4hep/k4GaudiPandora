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
 *  @file   DDMarlinPandora/src/DDTrackCreatorCLIC.cc
 *
 *  @brief  Implementation of the track creator class for a CLIC all silicon tracker.
 *
 *  $Log: $
 */

#include "DDTrackCreatorCLIC.h"

#include <LCObjects/LCTrack.h>
#include <Pandora/PdgTable.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetectorSelector.h>
#include <DDRec/DetectorData.h>

#include <edm4hep/TrackCollection.h>

#include <k4Interface/IGeoSvc.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// forward declarations. See in DDPandoraPFANewProcessor.cc
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorCLIC::DDTrackCreatorCLIC(const Settings& settings, pandora::Pandora& pandora,
                                       const Gaudi::Algorithm* thisAlg, SmartIF<IGeoSvc> geoSvc)
    : DDTrackCreatorBase(settings, pandora, thisAlg), m_trackerInnerR(0.f), m_trackerOuterR(0.f), m_trackerZmax(0.f),
      m_cosTracker(0.f), m_endcapDiskInnerRadii(DoubleVector()), m_endcapDiskOuterRadii(DoubleVector()),
      m_endcapDiskZPositions(DoubleVector()), m_nEndcapDiskLayers(0), m_barrelTrackerLayers(0),
      m_tanLambdaEndcapDisk(0.f), m_geoSvc(geoSvc)

{
  m_trackerInnerR = getTrackingRegionExtent()[0];
  m_trackerOuterR = getTrackingRegionExtent()[1];
  m_trackerZmax = getTrackingRegionExtent()[2];

  /// FIXME: Probably need to be something relating to last disk inner radius
  m_cosTracker = m_trackerZmax / std::sqrt(m_trackerZmax * m_trackerZmax + m_trackerInnerR * m_trackerInnerR);

  dd4hep::Detector* mainDetector = m_geoSvc->getDetector();

  // Maybe we need to veto the vertex? That was done in the ILD case
  const std::vector<dd4hep::DetElement>& barrelDets =
      dd4hep::DetectorSelector(*mainDetector).detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL));

  m_barrelTrackerLayers = 0;
  for (const auto& theDetector : barrelDets) {
    try {
      dd4hep::rec::ZPlanarData* theExtension = nullptr;

      theExtension = theDetector.extension<dd4hep::rec::ZPlanarData>();

      size_t N = theExtension->layers.size();
      m_barrelTrackerLayers = m_barrelTrackerLayers + N;

      m_algorithm.debug() << " Adding layers for barrel tracker from DD4hep for " << theDetector.name()
                          << "- n layers: " << N << " sum up to now: " << m_barrelTrackerLayers << endmsg;
    } catch (std::runtime_error& exception) {
      m_algorithm.warning() << "DDTrackCreatorCLIC exception during Barrel Tracker layer sum for " << theDetector.name()
                            << " : " << exception.what() << endmsg;
    }
  }

  m_nEndcapDiskLayers = 0;
  m_endcapDiskInnerRadii.clear();
  m_endcapDiskOuterRadii.clear();

  // Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD

  const std::vector<dd4hep::DetElement>& endcapDets =
      dd4hep::DetectorSelector(*mainDetector).detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP));

  for (auto iter = endcapDets.begin(), iterEnd = endcapDets.end(); iter != iterEnd; ++iter) {
    try {
      dd4hep::rec::ZDiskPetalsData* theExtension = nullptr;

      const dd4hep::DetElement& theDetector = *iter;
      theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();

      size_t N = theExtension->layers.size();

      m_algorithm.debug() << " Filling FTD-like parameters from DD4hep for " << theDetector.name()
                          << "- n layers: " << N << endmsg;

      for (unsigned int i = 0; i < N; ++i) {
        dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer = theExtension->layers[i];

        // Create a disk to represent even number petals front side
        // FIXME! VERIFY THAT TIS MAKES SENSE!
        m_endcapDiskInnerRadii.push_back(thisLayer.distanceSensitive / dd4hep::mm);
        m_endcapDiskOuterRadii.push_back(thisLayer.distanceSensitive / dd4hep::mm +
                                         thisLayer.lengthSensitive / dd4hep::mm);

        // Take the mean z position of the staggered petals
        const double zpos(thisLayer.zPosition / dd4hep::mm);
        m_endcapDiskZPositions.push_back(zpos);

        m_algorithm.debug() << "     layer " << i << " - mean z position = " << zpos << endmsg;
      }

      m_nEndcapDiskLayers = m_endcapDiskZPositions.size();

    } catch (std::runtime_error& exception) {
      m_algorithm.warning()
          << "DDTrackCreatorCLIC exception during Forward Tracking Disk parameter construction for detector "
          << const_cast<dd4hep::DetElement&>(*iter).name() << " : " << exception.what() << endmsg;
    }
  }

  for (unsigned int iEndcapDiskLayer = 0; iEndcapDiskLayer < m_nEndcapDiskLayers; ++iEndcapDiskLayer) {
    if (std::fabs(m_endcapDiskOuterRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon() ||
        std::fabs(m_endcapDiskInnerRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon()) {
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
  }

  m_tanLambdaEndcapDisk = m_endcapDiskZPositions[0] / m_endcapDiskOuterRadii[0];
}

pandora::StatusCode DDTrackCreatorCLIC::CreateTracks(const std::vector<edm4hep::Track>& tracks) {
  /// FIXME: Should really move to using surfaces
  for (const auto& pTrack : tracks) {
    m_algorithm.debug()
        << " Warning! Ignoring expected number of hits and other hit number cuts. Should eventually change!" << endmsg;

    // Take the first track state for the parameters
    if (pTrack.getTrackStates().empty()) {
      continue;
      // throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
    const auto& trackState = pTrack.getTrackStates()[0];

    // Proceed to create the pandora track
    lc_content::LCTrackParameters trackParameters;
    trackParameters.m_d0 = trackState.D0;
    trackParameters.m_z0 = trackState.Z0;
    trackParameters.m_pParentAddress = &pTrack;

    // By default, assume tracks are charged pions
    const float signedCurvature = trackState.omega;
    trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

    // Use particle id information from V0 and Kink finders
    // TODO: fix ID?
    auto trackPIDiter = m_trackToPidMap.find(GetTrackID(pTrack));

    if (trackPIDiter != m_trackToPidMap.end()) {
      trackParameters.m_particleId = trackPIDiter->second;
      trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackPIDiter->second);
    }

    if (signedCurvature != 0.f)
      trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

    try {
      // FIXME: Should consider adding a try-catch block to check
      // against cases where a track has invalid parameters like
      // very small omega. Especially for omega<epsilon=1e-7,
      // pandora throws a pandora::StatusCodeException
      //&statusCodeException which is not caught.
      // FIXED: if an exception happens we ignore the track

      GetTrackStates(pTrack, trackParameters);
      TrackReachesECAL(pTrack, trackParameters);
      GetTrackStatesAtCalo(pTrack, trackParameters);
      DefineTrackPfoUsage(pTrack, trackParameters);

      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                              PandoraApi::Track::Create(m_pandora, trackParameters, *m_lcTrackFactory))
      m_trackVector.push_back(pTrack);
    } catch (pandora::StatusCodeException& statusCodeException) {
      m_algorithm.error() << "Failed to extract a track: " << statusCodeException.ToString() << endmsg;
      m_algorithm.debug() << " failed track : " << pTrack << endmsg;
    }
  }

  m_algorithm.debug() << "After treating the input with " << tracks.size() << " tracks, the track vector size is "
                      << m_trackVector.size() << endmsg;

  return pandora::STATUS_CODE_SUCCESS;
}

bool DDTrackCreatorCLIC::PassesQualityCuts(const edm4hep::Track& pTrack,
                                           const PandoraApi::Track::Parameters& trackParameters) const {
  // First simple sanity checks
  if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() <
      m_settings.m_minTrackECalDistanceFromIp) {
    m_algorithm.warning() << " Dropping track! Distance at ECAL: "
                          << trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() << endmsg;
    m_algorithm.debug() << " track : " << pTrack << endmsg;
    return false;
  }

  const auto& firstTrackState = pTrack.getTrackStates(0);

  if (std::fabs(firstTrackState.omega) < std::numeric_limits<float>::epsilon()) {
    m_algorithm.error() << "Track has Omega = 0 " << endmsg;
    return false;
  }

  // TODO: Add this? It's in the DDTrackCreatorILD version
  // if (pTrack.getNdf() < 0) {
  //   m_algorithm.error() << "Track is unconstrained - ndf = " << pTrack.getNdf() << endmsg;
  //   return false;
  // }

  // Check momentum uncertainty is reasonable to use track
  const pandora::CartesianVector& momentumAtDca = trackParameters.m_momentumAtDca.Get();
  const float sigmaPOverP =
      std::sqrt(firstTrackState.getCovMatrix(edm4hep::TrackParams::omega, edm4hep::TrackParams::omega)) /
      std::fabs(firstTrackState.omega);

  if (sigmaPOverP > m_settings.m_maxTrackSigmaPOverP) {
    m_algorithm.warning() << " Dropping track : " << momentumAtDca.GetMagnitude() << "+-"
                          << sigmaPOverP * (momentumAtDca.GetMagnitude()) << " chi2 = " << pTrack.getChi2() << " "
                          << pTrack.getNdf() << " from " << pTrack.getTrackerHits().size() << endmsg;

    m_algorithm.debug() << " track : " << pTrack << endmsg;
    return false;
  }

  m_algorithm.debug() << " TEMPORARILY ACCEPT TRACK WITHOUT CUTS (should change!)" << pTrack << endmsg;
  return true;

  // Require reasonable number of Tracker hits
  if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks) {
    const float pZ = std::abs(momentumAtDca.GetZ());
    const float pT = std::hypot(std::abs(momentumAtDca.GetX()), std::abs(momentumAtDca.GetY()));
    // TODO: Compute
    const float rInnermostHit = 0;
    throw std::runtime_error("DDTrackCreatorCLIC::PassesQualityCuts: rInnermostHit not implemented");

    if (0.f == pT || 0.f == pZ || rInnermostHit == m_trackerOuterR) {
      m_algorithm.error() << "Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit "
                          << rInnermostHit << endmsg;
      return false;
    }

    float nExpectedTrackerHits = 0.;

    if (pZ < m_trackerZmax / m_trackerOuterR * pT) {
      const float innerExpectedHitRadius = std::max(m_trackerInnerR, rInnermostHit);
      const float frac = (m_trackerOuterR - innerExpectedHitRadius) / (m_trackerOuterR - m_trackerInnerR);
      nExpectedTrackerHits = m_barrelTrackerLayers * frac;
    }

    if (pZ <= m_trackerZmax / m_trackerInnerR * pT && pZ >= m_trackerZmax / m_trackerOuterR * pT) {
      const float innerExpectedHitRadius = std::max(m_trackerInnerR, rInnermostHit);
      const float frac =
          (m_trackerZmax * pT / pZ - innerExpectedHitRadius) / (m_trackerOuterR - innerExpectedHitRadius);
      nExpectedTrackerHits = frac * m_barrelTrackerLayers;
    }

    const auto& hitsBySubdetector(pTrack.getSubdetectorHitNumbers());

    // ---- use hitsInFit :
    // Use flags to access detectors and then access their ids to get hits in fit by ID
    //(new way of storing hits in MarlinTrkUtils)

    // Initialize hits to 0
    int nBarrelTrackerHits = 0;
    dd4hep::Detector* mainDetector = m_geoSvc->getDetector();
    const std::vector<dd4hep::DetElement>& barrelDets =
        dd4hep::DetectorSelector(*mainDetector).detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL));
    for (std::vector<dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();
         iter != iterEnd; ++iter) {
      const dd4hep::DetElement& theDetector = *iter;
      int detId = theDetector.id();
      nBarrelTrackerHits += hitsBySubdetector[2 * detId - 2]; // Offset is 2 for hits in fit
    }

    int nEndcapTrackerHits = 0;
    const std::vector<dd4hep::DetElement>& endcapDets =
        dd4hep::DetectorSelector(*mainDetector).detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP));
    for (std::vector<dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();
         iter != iterEnd; ++iter) {
      const dd4hep::DetElement& theDetector = *iter;
      int detId = theDetector.id();
      nEndcapTrackerHits += hitsBySubdetector[2 * detId - 2]; // Offset is 2 for hits in fit
    }

    const int minTrackerHits =
        static_cast<int>(nExpectedTrackerHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);

    if (nBarrelTrackerHits < minTrackerHits &&
        nEndcapTrackerHits < m_settings.m_minFtdHitsForBarrelTrackerHitFraction) {
      m_algorithm.warning() << " Dropping track : " << momentumAtDca.GetMagnitude()
                            << " Number of Tracker hits = " << nBarrelTrackerHits << " < " << minTrackerHits
                            << " nendcapDisk = " << nEndcapTrackerHits << endmsg;

      m_algorithm.debug() << " track : " << pTrack << endmsg;
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorCLIC::DefineTrackPfoUsage(const edm4hep::Track& pTrack,
                                             PandoraApi::Track::Parameters& trackParameters) const {
  bool canFormPfo = false;
  bool canFormClusterlessPfo = false;

  if (trackParameters.m_reachesCalorimeter.Get() && !IsParent(pTrack)) {
    const auto& firstTrackState = pTrack.getTrackStates(0);

    double rInner = std::numeric_limits<double>::max();
    double zMin = std::numeric_limits<double>::max();

    for (const auto& hit : pTrack.getTrackerHits()) {
      const auto pPosition = hit.getPosition();
      rInner = std::min(rInner, std::hypot(pPosition[0], pPosition[1]));
      zMin = std::min(zMin, std::fabs(pPosition[2]));
    }

    if (PassesQualityCuts(pTrack, trackParameters)) {
      const pandora::CartesianVector& momentumAtDca = trackParameters.m_momentumAtDca.Get();
      const float pZ = momentumAtDca.GetZ();
      const float pT = std::hypot(momentumAtDca.GetX(), momentumAtDca.GetY());

      const float zCutForNonVertexTracks = m_trackerInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks;
      const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) &&
                                   (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

      const bool isV0 = IsV0(pTrack);
      const bool isDaughter = IsDaughter(pTrack);

      m_algorithm.debug() << " -- track passed quality cuts and has : "
                          << " passRzQualityCuts " << passRzQualityCuts << " isV0 " << isV0 << " isDaughter "
                          << isDaughter << endmsg;

      // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
      if ((std::abs(firstTrackState.D0) < m_settings.m_d0TrackCut) &&
          (std::abs(firstTrackState.Z0) < m_settings.m_z0TrackCut) &&
          (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance)) {
        canFormPfo = true;
      } else if (passRzQualityCuts && m_settings.m_usingNonVertexTracks != 0) {
        canFormPfo = true;
      } else if (isV0 || isDaughter) {
        canFormPfo = true;
      }

      // Decide whether track can be used to form a charged PFO, even if track fails to be associated with a pandora
      // cluster
      const float particleMass = trackParameters.m_mass.Get();
      const float trackEnergy = std::sqrt(momentumAtDca.GetMagnitudeSquared() + particleMass * particleMass);

      if (m_settings.m_usingUnmatchedVertexTracks != 0 && trackEnergy < m_settings.m_unmatchedVertexTrackMaxEnergy) {
        if (std::abs(firstTrackState.D0) < m_settings.m_d0UnmatchedVertexTrackCut &&
            std::abs(firstTrackState.Z0) < m_settings.m_z0UnmatchedVertexTrackCut &&
            rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance) {
          canFormClusterlessPfo = true;
        } else if (passRzQualityCuts && m_settings.m_usingNonVertexTracks != 0 &&
                   m_settings.m_usingUnmatchedNonVertexTracks != 0) {
          canFormClusterlessPfo = true;
        } else if (isV0 || isDaughter) {
          canFormClusterlessPfo = true;
        }
      }
    } else if (IsDaughter(pTrack) || IsV0(pTrack)) {
      m_algorithm.warning() << "Recovering daughter or v0 track "
                            << trackParameters.m_momentumAtDca.Get().GetMagnitude() << endmsg;
      canFormPfo = true;
    }

    m_algorithm.debug() << " -- track canFormPfo = " << canFormPfo
                        << " -  canFormClusterlessPfo = " << canFormClusterlessPfo << endmsg;
  }

  trackParameters.m_canFormPfo = canFormPfo;
  trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorCLIC::TrackReachesECAL(const edm4hep::Track& pTrack,
                                          PandoraApi::Track::Parameters& trackParameters) const {
  // Calculate hit position information
  float hitZMin = std::numeric_limits<float>::max();
  float hitZMax = -std::numeric_limits<float>::max();
  float hitOuterR = -std::numeric_limits<float>::max();

  int nTrackerHits = 0;
  int nEndcapDiskHits = 0;
  int maxOccupiedEndcapDiskLayer = 0;

  for (const auto& hit : pTrack.getTrackerHits()) {
    const float z = static_cast<float>(hit.getPosition()[2]);
    const float r = std::hypot(hit.getPosition()[0], hit.getPosition()[1]);

    hitZMin = std::min(hitZMin, z);
    hitZMax = std::max(hitZMax, z);
    hitOuterR = std::max(hitOuterR, r);

    if (r > m_trackerInnerR && r < m_trackerOuterR && std::fabs(z) <= m_trackerZmax) {
      nTrackerHits++;
      continue;
    }

    for (unsigned int j = 0; j < m_nEndcapDiskLayers; ++j) {
      if (r > m_endcapDiskInnerRadii[j] && r < m_endcapDiskOuterRadii[j] &&
          std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_endcapDiskZPositions[j] &&
          std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_endcapDiskZPositions[j]) {
        if (static_cast<int>(j) > maxOccupiedEndcapDiskLayer)
          maxOccupiedEndcapDiskLayer = static_cast<int>(j);

        nEndcapDiskHits++;
        break;
      }
    }
  }

  // Require sufficient hits in barrel or endcap trackers, then compare extremal hit positions with tracker dimensions
  if (nTrackerHits >= m_settings.m_reachesECalNBarrelTrackerHits ||
      nEndcapDiskHits >= m_settings.m_reachesECalNFtdHits) {
    if (hitOuterR - m_trackerOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance ||
        std::fabs(hitZMax) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance ||
        std::fabs(hitZMin) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance ||
        maxOccupiedEndcapDiskLayer >= m_settings.m_reachesECalMinFtdLayer) {
      trackParameters.m_reachesCalorimeter = true;
      return;
    }
  }

  // If track is lowpt, it may curl up and end inside tpc inner radius
  const pandora::CartesianVector& momentumAtDca = trackParameters.m_momentumAtDca.Get();
  const float cosAngleAtDca = std::fabs(momentumAtDca.GetZ()) / momentumAtDca.GetMagnitude();
  const float pT = std::hypot(momentumAtDca.GetX(), momentumAtDca.GetY());

  if (cosAngleAtDca > m_cosTracker ||
      pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField * m_trackerOuterR) {
    trackParameters.m_reachesCalorimeter = true;
    return;
  }

  trackParameters.m_reachesCalorimeter = false;
}
