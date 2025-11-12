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
 *  @file   MarlinPandora/src/DDTrackCreatorILD.cc
 *
 *  @brief  Implementation of the track creator class for ILD.
 *
 *  $Log: $
 */

#include "DDTrackCreatorILD.h"
#include "Pandora/PdgTable.h"

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetectorSelector.h>
#include <DDRec/DetectorData.h>

#include <LCObjects/LCTrack.h>
#include <edm4hep/TrackCollection.h>

#include <k4Interface/IGeoSvc.h>

#include <Gaudi/Algorithm.h>

#include <algorithm>
#include <cmath>
#include <limits>

DDTrackCreatorILD::DDTrackCreatorILD(const Settings& settings, pandora::Pandora& pPandora,
                                     const Gaudi::Algorithm* thisAlg, SmartIF<IGeoSvc> geoSvc)
    : DDTrackCreatorBase(settings, pPandora, thisAlg), m_cosTpc(0.f), m_tpcInnerR(0.f), m_tpcOuterR(0.f),
      m_tpcMaxRow(0), m_tpcZmax(0.f), m_tpcMembraneMaxZ(0.f), m_ftdInnerRadii(DoubleVector()),
      m_ftdOuterRadii(DoubleVector()), m_ftdZPositions(DoubleVector()), m_nFtdLayers(0), m_tanLambdaFtd(0.f),
      m_minEtdZPosition(0.f), m_minSetRadius(0.f), m_geoSvc(geoSvc)

{
  m_nFtdLayers = 0;
  m_ftdInnerRadii.clear();
  m_ftdOuterRadii.clear();

  // Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD
  /// FIXME: Should we use surfaces instead?
  dd4hep::Detector* mainDetector = m_geoSvc->getDetector();
  const std::vector<dd4hep::DetElement>& endcapDets =
      dd4hep::DetectorSelector(*mainDetector).detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP));
  for (std::vector<dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();
       iter != iterEnd; ++iter) {
    try {
      dd4hep::rec::ZDiskPetalsData* theExtension = nullptr;

      const dd4hep::DetElement& theDetector = *iter;
      theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();

      unsigned int N = theExtension->layers.size();

      m_algorithm.debug() << " Filling FTD-like parameters from DD4hep for " << theDetector.name()
                          << "- n layers: " << N << endmsg;

      for (unsigned int i = 0; i < N; ++i) {
        dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer = theExtension->layers[i];

        // Create a disk to represent even number petals front side
        // FIXME! VERIFY THAT TIS MAKES SENSE!
        m_ftdInnerRadii.push_back(thisLayer.distanceSensitive / dd4hep::mm);
        m_ftdOuterRadii.push_back(thisLayer.distanceSensitive / dd4hep::mm + thisLayer.lengthSensitive / dd4hep::mm);

        // Take the mean z position of the staggered petals
        const double zpos(thisLayer.zPosition / dd4hep::mm);
        m_ftdZPositions.push_back(zpos);

        m_algorithm.debug() << "     layer " << i << " - mean z position = " << zpos << endmsg;
      }

      m_nFtdLayers = m_ftdZPositions.size();

    } catch (std::runtime_error& exception) {
      m_algorithm.warning()
          << "DDTrackCreatorILD exception during Forward Tracking Disk parameter construction for detector "
          << const_cast<dd4hep::DetElement&>(*iter).name() << " : " << exception.what() << endmsg;
    }
  }

  try {
    dd4hep::rec::FixedPadSizeTPCData* theExtension = nullptr;
    // Get the TPC, make sure not to get the vertex
    const std::vector<dd4hep::DetElement>& tpcDets =
        dd4hep::DetectorSelector(*mainDetector)
            .detectors((dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL | dd4hep::DetType::GASEOUS),
                       dd4hep::DetType::VERTEX);

    // There should only be one TPC
    theExtension = tpcDets[0].extension<dd4hep::rec::FixedPadSizeTPCData>();

    m_tpcInnerR = theExtension->rMin / dd4hep::mm;
    m_tpcOuterR = theExtension->rMax / dd4hep::mm;
    m_tpcMaxRow = theExtension->maxRow;
    m_tpcZmax = theExtension->zHalf / dd4hep::mm;

    // FIXME! Add to reco structrure and access
    m_tpcMembraneMaxZ = 10;
    std::cout << "WARNING! DO NOT MASK! HANDLE m_tpcMembraneMaxZ (currently hardcoded to 10)!" << std::endl;

  } catch (std::runtime_error&) {
    m_algorithm.warning() << "DDTrackCreatorILD exception during TPC parameter construction." << endmsg;
  }

  // Check tpc parameters
  if (std::fabs(m_tpcZmax) < std::numeric_limits<float>::epsilon() ||
      std::fabs(m_tpcInnerR) < std::numeric_limits<float>::epsilon() ||
      std::fabs(m_tpcOuterR - m_tpcInnerR) < std::numeric_limits<float>::epsilon()) {
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  m_cosTpc = m_tpcZmax / std::sqrt(m_tpcZmax * m_tpcZmax + m_tpcInnerR * m_tpcInnerR);

  // Check ftd parameters
  if (0 == m_nFtdLayers || m_nFtdLayers != m_ftdInnerRadii.size() || m_nFtdLayers != m_ftdOuterRadii.size()) {
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }

  for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer) {
    if (std::fabs(m_ftdOuterRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon() ||
        std::fabs(m_ftdInnerRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon()) {
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
  }

  m_tanLambdaFtd = m_ftdZPositions[0] / m_ftdOuterRadii[0];

  // Calculate etd and set parameters
  // fg: make SET and ETD optional - as they might not be in the model ...
  // FIXME: THINK OF A UNIVERSAL WAY TO HANDLE EXISTENCE OF ADDITIONAL DETECTORS

  m_algorithm.warning() << " ETDLayerZ or SETLayerRadius parameters Not being handled!" << endmsg
                        << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() << endmsg;
  m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
  m_minSetRadius = std::numeric_limits<float>::quiet_NaN();

  //     try
  //     {
  //         const DoubleVector
  //         &etdZPositions(marlin::Global::GEAR->getGearParameters("ETD").getDoubleVals("ETDLayerZ")); const
  //         DoubleVector
  //         &setInnerRadii(marlin::Global::GEAR->getGearParameters("SET").getDoubleVals("SETLayerRadius"));
  //
  //         if (etdZPositions.empty() || setInnerRadii.empty())
  //             throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  //
  //         m_minEtdZPosition = *(std::min_element(etdZPositions.begin(), etdZPositions.end()));
  //         m_minSetRadius = *(std::min_element(setInnerRadii.begin(), setInnerRadii.end()));
  //     }
  //     catch(gear::UnknownParameterException &)
  //     {
  //         m_algorithm.warning() << " ETDLayerZ or SETLayerRadius parameters missing from GEAR parameters!" <<
  //         endmsg
  //                                << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() <<
  //                                endmsg;
  //
  //         //fg: Set them to NAN, so that they cannot be used to set   trackParameters.m_reachesCalorimeter = true;
  //         m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
  //         m_minSetRadius = std::numeric_limits<float>::quiet_NaN();
  //     }
}

pandora::StatusCode DDTrackCreatorILD::CreateTracks(const std::vector<edm4hep::Track>& tracks) {
  for (const auto& pTrack : tracks) {

    if (pTrack.getTrackStates().empty()) {
      continue;
      // throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    const auto& trackState = pTrack.getTrackStates()[0];

    int minTrackHits = m_settings.m_minTrackHits;
    const float tanLambda = std::fabs(trackState.tanLambda);

    if (tanLambda > m_tanLambdaFtd) {
      int expectedFtdHits = 0;

      for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer) {
        if ((tanLambda > m_ftdZPositions[iFtdLayer] / m_ftdOuterRadii[iFtdLayer]) &&
            (tanLambda < m_ftdZPositions[iFtdLayer] / m_ftdInnerRadii[iFtdLayer])) {
          expectedFtdHits++;
        }
      }

      minTrackHits = std::max(m_settings.m_minFtdTrackHits, expectedFtdHits);
    }

    const int nTrackHits = static_cast<int>(pTrack.getTrackerHits().size());

    if ((nTrackHits < minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
      continue;

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

    if (std::numeric_limits<float>::epsilon() < std::fabs(signedCurvature))
      trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

    try // fg: include the next calls in the try block to catch tracks that are yet not fitted properly as ERROR and
        // not exceptions
    {
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

  return pandora::STATUS_CODE_SUCCESS;
}

bool DDTrackCreatorILD::PassesQualityCuts(const edm4hep::Track& pTrack,
                                          const PandoraApi::Track::Parameters& trackParameters) const {
  // First simple sanity checks
  if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() <
      m_settings.m_minTrackECalDistanceFromIp) {
    m_algorithm.warning() << " Dropping track! Distance at ECAL: "
                          << trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() << endmsg;
    m_algorithm.debug() << " track : " << pTrack << endmsg;
    return false;
  }

  const auto& firstTrackState = pTrack.getTrackStates()[0];

  if (std::fabs(firstTrackState.omega) < std::numeric_limits<float>::epsilon()) {
    m_algorithm.error() << "Track has Omega = 0 " << endmsg;
    return false;
  }

  if (pTrack.getNdf() < 0) {
    m_algorithm.error() << "Track is unconstrained - ndf = " << pTrack.getNdf() << endmsg;
    return false;
  }

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

  // Require reasonable number of TPC hits
  if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks) {
    const float pZ = std::abs(momentumAtDca.GetZ());
    const float pT = std::hypot(std::abs(momentumAtDca.GetX()), std::abs(momentumAtDca.GetY()));
    // TODO: Compute
    const float rInnermostHit = 0;
    throw std::runtime_error("DDTrackCreatorILD::PassesQualityCuts: rInnermostHit not implemented");

    if (std::numeric_limits<float>::epsilon() > std::fabs(pT) ||
        std::numeric_limits<float>::epsilon() > std::fabs(pZ) || rInnermostHit == m_tpcOuterR) {
      m_algorithm.error() << "Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit "
                          << rInnermostHit << endmsg;
      return false;
    }

    float nExpectedTpcHits = 0.;

    if (pZ < m_tpcZmax / m_tpcOuterR * pT) {
      const float innerExpectedHitRadius = std::max(m_tpcInnerR, rInnermostHit);
      const float frac = (m_tpcOuterR - innerExpectedHitRadius) / (m_tpcOuterR - m_tpcInnerR);
      nExpectedTpcHits = m_tpcMaxRow * frac;
    }

    if (pZ <= m_tpcZmax / m_tpcInnerR * pT && pZ >= m_tpcZmax / m_tpcOuterR * pT) {
      const float innerExpectedHitRadius = std::max(m_tpcInnerR, rInnermostHit);
      const float frac = (m_tpcZmax * pT / pZ - innerExpectedHitRadius) / (m_tpcOuterR - innerExpectedHitRadius);
      nExpectedTpcHits = frac * m_tpcMaxRow;
    }

    // TODO Get TPC membrane information from GEAR when available
    if (std::fabs(pZ) / momentumAtDca.GetMagnitude() < m_tpcMembraneMaxZ / m_tpcInnerR)
      nExpectedTpcHits = 0;

    const int nTpcHits = GetNTpcHits(pTrack);
    const int nFtdHits = GetNFtdHits(pTrack);

    const int minTpcHits = static_cast<int>(nExpectedTpcHits * m_settings.m_minBarrelTrackerHitFractionOfExpected);

    if (nTpcHits < minTpcHits && nFtdHits < m_settings.m_minFtdHitsForBarrelTrackerHitFraction) {
      m_algorithm.warning() << " Dropping track : " << momentumAtDca.GetMagnitude()
                            << " Number of TPC hits = " << nTpcHits << " < " << minTpcHits << " nftd = " << nFtdHits
                            << endmsg;

      m_algorithm.debug() << " track : " << pTrack << endmsg;
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorILD::DefineTrackPfoUsage(const edm4hep::Track& pTrack,
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

      const float zCutForNonVertexTracks = m_tpcInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks;
      const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) &&
                                   (rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

      const bool isV0 = IsV0(pTrack);
      const bool isDaughter = IsDaughter(pTrack);

      m_algorithm.debug() << " -- track passed quality cuts and has : "
                          << " passRzQualityCuts " << passRzQualityCuts << " isV0 " << isV0 << " isDaughter "
                          << isDaughter << endmsg;

      // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
      if ((firstTrackState.D0 < m_settings.m_d0TrackCut) && (firstTrackState.Z0 < m_settings.m_z0TrackCut) &&
          (rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance)) {
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
        if (firstTrackState.D0 < m_settings.m_d0UnmatchedVertexTrackCut &&
            firstTrackState.Z0 < m_settings.m_z0UnmatchedVertexTrackCut &&
            rInner < m_tpcInnerR + m_settings.m_maxBarrelTrackerInnerRDistance) {
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

void DDTrackCreatorILD::TrackReachesECAL(const edm4hep::Track&, PandoraApi::Track::Parameters& trackParameters) const {
  // fg: return true  for now - there are quality checks in DefineTrackPfoUsage() ...
  trackParameters.m_reachesCalorimeter = true;
  // Everything after the return has been commented out since it is not executed
  return;

  // // at a later stage we could simply check if there is a valid track state at the calorimeter:
  // // ...

  // // Calculate hit position information
  // float hitZMin(std::numeric_limits<float>::max());
  // float hitZMax(-std::numeric_limits<float>::max());
  // float hitOuterR(-std::numeric_limits<float>::max());

  // int maxOccupiedFtdLayer(0);

  // const edm4hep::TrackerHitVec& trackerHitVec(pTrack->getTrackerHits());
  // const unsigned int nTrackHits(trackerHitVec.size());

  // for (const auto& hit : pTrack.getTrackerHits()) {
  //  const float z = static_cast<float>(hit.getPosition()[2]);
  //  const float r = std::hypot(hit.getPosition()[0], hit.getPosition()[1]);

  //   hitZMax = std::max(hitZMax, z);
  //   hitZMin = std::min(hitZMin, z);
  //   hitOuterR = std::max(hitOuterR, r);

  //   if ((r > m_tpcInnerR) && (r < m_tpcOuterR) && (std::fabs(z) <= m_tpcZmax)) {
  //     continue;
  //   }

  //   for (unsigned int j = 0; j < m_nFtdLayers; ++j) {
  //     if ((r > m_ftdInnerRadii[j]) && (r < m_ftdOuterRadii[j]) &&
  //         (std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_ftdZPositions[j]) &&
  //         (std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_ftdZPositions[j])) {
  //       if (static_cast<int>(j) > maxOccupiedFtdLayer)
  //         maxOccupiedFtdLayer = static_cast<int>(j);

  //       break;
  //     }
  //   }
  // }

  // const int nTpcHits(this->GetNTpcHits(pTrack));
  // const int nFtdHits(this->GetNFtdHits(pTrack));

  // // Look to see if there are hits in etd or set, implying track has reached edge of ecal
  // if ((hitOuterR > m_minSetRadius) || (hitZMax > m_minEtdZPosition)) {
  //   trackParameters.m_reachesCalorimeter = true;
  //   return;
  // }

  // // Require sufficient hits in tpc or ftd, then compare extremal hit positions with tracker dimensions
  // if ((nTpcHits >= m_settings.m_reachesECalNBarrelTrackerHits) || (nFtdHits >= m_settings.m_reachesECalNFtdHits)) {
  //   if ((hitOuterR - m_tpcOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance) ||
  //       (std::fabs(hitZMax) - m_tpcZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
  //       (std::fabs(hitZMin) - m_tpcZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
  //       (maxOccupiedFtdLayer >= m_settings.m_reachesECalMinFtdLayer)) {
  //     trackParameters.m_reachesCalorimeter = true;
  //     return;
  //   }
  // }

  // // If track is lowpt, it may curl up and end inside tpc inner radius
  // const pandora::CartesianVector& momentumAtDca(trackParameters.m_momentumAtDca.Get());
  // const float cosAngleAtDca(std::fabs(momentumAtDca.GetZ()) / momentumAtDca.GetMagnitude());
  // const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY());
  // const float pT(std::sqrt(pX * pX + pY * pY));

  // if ((cosAngleAtDca > m_cosTpc) || (pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField *
  // m_tpcOuterR)) {
  //   trackParameters.m_reachesCalorimeter = true;
  //   return;
  // }

  // trackParameters.m_reachesCalorimeter = false;
}

int DDTrackCreatorILD::GetNTpcHits(const edm4hep::Track& pTrack) const {
  // ATTN
  // According to FG: [ 2 * lcio::ILDDetID::TPC - 2 ] is the first number and it is supposed to
  // be the number of hits in the fit and this is what should be used !
  // at least for DD4hep/DDSim

  // ---- use hitsInFit :
  // return pTrack.getSubdetectorHitNumbers()[2 * lcio::ILDDetID::TPC - 2];
  return pTrack.getSubdetectorHitNumbers()[2 * 4 - 2];
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DDTrackCreatorILD::GetNFtdHits(const edm4hep::Track& pTrack) const {
  // ATTN
  // see above
  // ---- use hitsInFit :
  // return pTrack.getSubdetectorHitNumbers()[2 * lcio::ILDDetID::FTD - 2];
  return pTrack.getSubdetectorHitNumbers()[2 * 4 - 2];
}

//
