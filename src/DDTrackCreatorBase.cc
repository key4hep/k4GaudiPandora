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
 *  @file   DDMarlinPandora/src/DDTrackCreatorBase.cc
 *
 *  @brief  Implementation of the track creator class.
 *
 *  $Log: $
 */

#include "edm4hep/ReconstructedParticle.h"


#include <LCObjects/LCTrack.h>
#include "DDTrackCreatorBase.h"
#include "Pandora/PdgTable.h"

#include <MarlinTrk/Factory.h>
#include <MarlinTrk/IMarlinTrack.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include <algorithm>
#include <cmath>
#include <limits>

//forward declaration
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorBase::DDTrackCreatorBase(const Settings& settings, const pandora::Pandora* const pPandora, MsgStream log)
    : m_settings(settings),
      m_pandora(*pPandora),
      m_trackVector(0),
      m_v0TrackList(TrackList()),
      m_parentTrackList(TrackList()),
      m_daughterTrackList(TrackList()),
      m_trackToPidMap(TrackToPidMap()),
      m_minimalTrackStateRadiusSquared(0.f),
      m_log(log) {
  const float ecalInnerR           = settings.m_eCalBarrelInnerR;
  const float tsTolerance          = settings.m_trackStateTolerance;
  m_minimalTrackStateRadiusSquared = (ecalInnerR - tsTolerance) * (ecalInnerR - tsTolerance);
  //wrap in shared_ptr with a dummy destructor
  m_trackingSystem = std::shared_ptr<MarlinTrk::IMarlinTrkSystem>(
      MarlinTrk::Factory::createMarlinTrkSystem(settings.m_trackingSystemName, nullptr, ""),
      [](MarlinTrk::IMarlinTrkSystem*) {});
  m_trackingSystem->init();
  m_encoder        = std::make_shared<BitField64>("encoding string");
  m_lcTrackFactory = std::make_shared<lc_content::LCTrackFactory>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorBase::~DDTrackCreatorBase() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::CreateTrackAssociations(
  const std::vector<const edm4hep::VertexCollection*>& kinkCollections,
  const std::vector<const edm4hep::VertexCollection*>& prongCollections,
  const std::vector<const edm4hep::VertexCollection*>& v0Collections
) {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(kinkCollections));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(prongCollections));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(v0Collections));

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractKinks(const std::vector<const edm4hep::VertexCollection*>& vertexCollections) {
  for (int colIndex = 0; colIndex < vertexCollections.size(); colIndex++) {
    try {
      const edm4hep::VertexCollection* pKinkCollection =vertexCollections[colIndex];

      for (int i = 0, iMax = pKinkCollection->size(); i < iMax; ++i) {
        try {
          edm4hep::Vertex* pVertex = dynamic_cast<Vertex*>(&(pKinkCollection->at(i)));

          if (NULL == pVertex)
            m_log << MSG::ERROR << "Collection type mismatch" << endmsg;

          edm4hep::ReconstructedParticle pReconstructedParticle = pVertex->getParticles(0);

          if (this->IsConflictingRelationship(pReconstructedParticle))
            continue;

          const int vertexPdgCode(pReconstructedParticle.getType());

          // Extract the kink vertex information
          for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack) {
            edm4hep::Track pTrack = pReconstructedParticle.getTracks(iTrack);
            (0 == iTrack) ? m_parentTrackList.insert(&pTrack) : m_daughterTrackList.insert(&pTrack);
            m_log << MSG::DEBUG << "KinkTrack " << iTrack << ", nHits " << pTrack.trackerHits_size()
                                 << endmsg;

            int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

            if (0 == iTrack) {
              trackPdgCode = vertexPdgCode;
            } else {
              switch (vertexPdgCode) {
                case pandora::PI_PLUS:
                case pandora::K_PLUS:
                  trackPdgCode = pandora::MU_PLUS;
                  break;
                case pandora::PI_MINUS:
                case pandora::K_MINUS:
                  trackPdgCode = pandora::MU_MINUS;
                  break;
                case pandora::HYPERON_MINUS_BAR:
                case pandora::SIGMA_PLUS:
                  trackPdgCode = pandora::PI_PLUS;
                  break;
                case pandora::SIGMA_MINUS:
                case pandora::HYPERON_MINUS:
                  trackPdgCode = pandora::PI_PLUS;
                  break;
                default:
                  (pTrack.getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                  break;
              }
            }

            m_trackToPidMap.insert(TrackToPidMap::value_type(&pTrack, trackPdgCode));

            if (0 == m_settings.m_shouldFormTrackRelationships)
              continue;

            // Make track parent-daughter relationships
            if (0 == iTrack) {
              for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
                PANDORA_RETURN_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=,
                    PandoraApi::SetTrackParentDaughterRelationship(m_pandora, pTrack, pReconstructedParticle.getTracks(jTrack)));
              }
            }

            // Make track sibling relationships
            else {
              for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                         PandoraApi::SetTrackSiblingRelationship(m_pandora, pTrack, pReconstructedParticle.getTracks(jTrack)));
              }
            }
          }
        } catch (...) {
          m_log << MSG::WARNING << "Failed to extract kink vertex: " << exception.what() << endmsg;
        }
      }
    } catch (...) {
      m_log << MSG::DEBUG << "Failed to extract kink vertex collection: " << *iter << ", " << exception.what()
                            << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractProngsAndSplits(const std::vector<const edm4hep::VertexCollection*>& vertexCollections) {
  for (int colIndex = 0; colIndex < vertexCollections.size(); colIndex++) {
    try {
      const edm4hep::VertexCollection* pProngOrSplitCollection = vertexCollections[colIndex];

      for (int i = 0, iMax = pProngOrSplitCollection->size(); i < iMax; ++i) {
        try {
          edm4hep::Vertex* pVertex = dynamic_cast<Vertex*>(&(pProngOrSplitCollection->at(i)));

          if (NULL == pVertex)
            m_log << MSG::ERROR << "Collection type mismatch" << endmsg;

          edm4hep::ReconstructedParticle pReconstructedParticle = pVertex->getParticles(0);

          if (this->IsConflictingRelationship(pReconstructedParticle))
            continue;

          // Extract the prong/split vertex information
          for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack) {
            edm4hep::Track pTrack = pReconstructedParticle.getTracks(iTrack);
            (0 == iTrack) ? m_parentTrackList.insert(&pTrack) : m_daughterTrackList.insert(&pTrack);
            m_log << MSG::DEBUG << "Prong or Split Track " << iTrack << ", nHits " << pTrack.trackerHits_size()
                                 << endmsg;

            if (0 == m_settings.m_shouldFormTrackRelationships)
              continue;

            // Make track parent-daughter relationships
            if (0 == iTrack) {
              for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
                PANDORA_RETURN_RESULT_IF(
                    pandora::STATUS_CODE_SUCCESS, !=,
                    PandoraApi::SetTrackParentDaughterRelationship(m_pandora, pTrack, pReconstructedParticle.getTracks(jTrack)));
              }
            }

            // Make track sibling relationships
            else {
              for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                         PandoraApi::SetTrackSiblingRelationship(m_pandora, pTrack, pReconstructedParticle.getTracks(jTrack)));
              }
            }
          }
        } catch (...) {
          m_log << MSG::WARNING << "Failed to extract prong/split vertex: " << exception.what() << endmsg;
        }
      }
    } catch (...) {
      m_log << MSG::DEBUG << "Failed to extract prong/split vertex collection: " << *iter << ", " << exception.what()
                            << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractV0s(const std::vector<const edm4hep::VertexCollection*>& vertexCollections) {
  for (int colIndex = 0; colIndex < vertexCollections.size(); colIndex++) {
    try {
      const edm4hep::VertexCollection* pProngOrSplitCollection = vertexCollections[colIndex];

      for (int i = 0, iMax = pV0Collection->size(); i < iMax; ++i) {
        try {
          edm4hep::Vertex* pVertex = dynamic_cast<Vertex*>(&(pV0Collection->at(i)));

          if (NULL == pVertex)
            m_log << MSG::ERROR << "Collection type mismatch" << endmsg;

          edm4hep::ReconstructedParticle pReconstructedParticle = pVertex->getParticles(0);

          if (this->IsConflictingRelationship(pReconstructedParticle))
            continue;

          // Extract the v0 vertex information
          const int vertexPdgCode(pReconstructedParticle.getType());

          for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack) {
            edm4hep::Track pTrack = pReconstructedParticle.getTracks(iTrack);
            m_v0TrackList.insert(&pTrack);
            m_log << MSG::DEBUG << "V0Track " << iTrack << ", nHits " << pTrack.trackerHits_size() << endmsg;

            int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

            switch (vertexPdgCode) {
              case pandora::PHOTON:
                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::E_PLUS : trackPdgCode = pandora::E_MINUS;
                break;
              case pandora::LAMBDA:
                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PROTON : trackPdgCode = pandora::PI_MINUS;
                break;
              case pandora::LAMBDA_BAR:
                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PROTON_BAR;
                break;
              case pandora::K_SHORT:
                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                break;
              default:
                (pTrack->getOmega() > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                break;
            }

            m_trackToPidMap.insert(TrackToPidMap::value_type(&pTrack, trackPdgCode));

            if (0 == m_settings.m_shouldFormTrackRelationships)
              continue;

            // Make track sibling relationships
            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
              PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                       PandoraApi::SetTrackSiblingRelationship(m_pandora, pTrack, pReconstructedParticle.getTracks(jTrack)));
            }
          }
        } catch (...) {
          m_log << MSG::WARNING << "Failed to extract v0 vertex: " << exception.what() << endmsg;
        }
      }
    } catch (...) {
      m_log << MSG::DEBUG << "Failed to extract v0 vertex collection: " << *iter << ", " << exception.what()
                            << endmsg;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDTrackCreatorBase::IsConflictingRelationship(const edm4hep::ReconstructedParticle pReconstructedParticle) const {
  for (edm4hep::Track pTrack : pReconstructedParticle.getTracks()) {

    if (this->IsDaughter(pTrack) || this->IsParent(pTrack) || this->IsV0(pTrack))
      return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::GetTrackStates(const edm4hep::Track           pTrack,
                                        PandoraApi::Track::Parameters& trackParameters) const {
  const edm4hep::TrackState *pTrackState = &(pTrack.getTrackState(edm4hep::TrackState::AtIP));

  if (!pTrackState)
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState->omega()));
  trackParameters.m_momentumAtDca =
      pandora::CartesianVector(std::cos(pTrackState->phi()), std::sin(pTrackState->phi()),
                               pTrackState->tanLambda()) *
      pt;

  this->CopyTrackState(&(pTrack.getTrackState(TrackState::AtFirstHit)), trackParameters.m_trackStateAtStart);

  //fg: curling TPC tracks have pointers to track segments stored -> need to get track states from last segment!
  const edm4hep::Track pEndTrack = (pTrack.getTracks().empty() ? pTrack : pTrack.getTracks().back());

  this->CopyTrackState(&(pEndTrack.getTrackState(TrackState::AtLastHit)), trackParameters.m_trackStateAtEnd);
  this->CopyTrackState(&(pEndTrack->getTrackState(TrackState::AtCalorimeter)), trackParameters.m_trackStateAtCalorimeter);

  trackParameters.m_isProjectedToEndCap =
      ((std::fabs(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetZ()) < m_settings.m_eCalEndCapInnerZ)
           ? false
           : true);

  // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
  const float minGenericTime(this->CalculateTrackTimeAtCalorimeter(pTrack));
  const float particleMass(trackParameters.m_mass.Get());
  const float particleEnergy(
      std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
  trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 299.792f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::GetTrackStatesAtCalo(edm4hep::Track track, lc_content::LCTrackParameters& trackParameters) {
  if (not trackParameters.m_reachesCalorimeter.Get()) {
    m_log << MSG::DEBUG << "Track does not reach the ECal" << endmsg;
    return;
  }

  const edm4hep::TrackState* trackAtCalo = &(track.getTrackState(TrackState::AtCalorimeter));
  if (!trackAtCalo) {
    m_log << MSG::DEBUG << "Track does not have a trackState at calorimeter" << endmsg;
    m_log << MSG::DEBUG << toString(track) << endmsg;
    return;
  }

  m_log << MSG::DEBUG << "Original" << toString(trackAtCalo) << endmsg;

  const auto* tsPosition = trackAtCalo->getReferencePoint();

  if (std::fabs(tsPosition[2]) < getTrackingRegionExtent()[2]) {
    m_log << MSG::DEBUG << "Original trackState is at Barrel" << endmsg;
    pandora::InputTrackState pandoraTrackState;
    this->CopyTrackState(trackAtCalo, pandoraTrackState);
    trackParameters.m_trackStates.push_back(pandoraTrackState);
  } else {  // if track state is in endcap we do not repeat track state calculation, because the barrel cannot be hit
    m_log << MSG::DEBUG << "Original track state is at Endcap" << endmsg;
    pandora::InputTrackState pandoraTrackState;
    this->CopyTrackState(trackAtCalo, pandoraTrackState);
    trackParameters.m_trackStates.push_back(pandoraTrackState);
    return;
  }

  auto marlintrk  = std::unique_ptr<MarlinTrk::IMarlinTrack>(m_trackingSystem->createTrack());

  for (int iHit = 0; iHit < track.trackerHits_size(); ++iHit) {
    edm4hep::TrackerHit trkHit = track.getTrackerHits(iHit);
    /* TODO this....
    if (UTIL::BitSet32(
            trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {  //it is a composite spacepoint
      //Split it up and add both hits to the MarlinTrk
      const EVENT::LCObjectVec& rawObjects = trkHit->getRawHits();
      for (unsigned k = 0; k < rawObjects.size(); k++) {
        EVENT::TrackerHit* rawHit = static_cast<EVENT::TrackerHit*>(rawObjects[k]);
        if (marlintrk->addHit(rawHit) != MarlinTrk::IMarlinTrack::success) {
          m_log << MSG::DEBUG << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add strip hit " << *rawHit
                                << endmsg;
        }
      }
    } else {*/
      if (marlintrk->addHit(track.getTrackerHits(iHit)) != MarlinTrk::IMarlinTrack::success)
        m_log << MSG::DEBUG << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add tracker hit " << *trkHit
                              << endmsg;
    //}
  }

  bool tanL_is_positive = trackAtCalo->tanLambda() > 0;

  auto trackState = edm4hep::MutableTrackState(*trackAtCalo);

  int return_error = marlintrk->initialise(trackState, m_settings.m_bField, MarlinTrk::IMarlinTrack::modeForward);
  if (return_error != MarlinTrk::IMarlinTrack::success) {
    m_log << MSG::DEBUG << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to initialize track for endcap track : "
                          << endmsg;
    return;
  }

  double chi2 = -DBL_MAX;
  int    ndf  = 0;

  edm4hep::MutableTrackState trackStateAtCaloEndcap;

  unsigned ecal_endcap_face_ID = lcio::ILDDetID::ECAL_ENDCAP;
  int      detElementID        = 0;
  /// TODO this.... 
  m_encoder->reset();  // reset to 0
  (*m_encoder)["subdet"] = ecal_endcap_face_ID;
  (*m_encoder)["side"]   = tanL_is_positive ? lcio::ILDDetID::fwd : lcio::ILDDetID::bwd;
  (*m_encoder)["layer"]  = 0;

  return_error = marlintrk->propagateToLayer(m_encoder->lowWord(), trackStateAtCaloEndcap, chi2, ndf, detElementID,
                                             MarlinTrk::IMarlinTrack::modeForward);
  m_log << MSG::DEBUG << "Found trackState at endcap? Error code: " << return_error << endmsg;

  if (return_error == MarlinTrk::IMarlinTrack::success) {
    m_log << MSG::DEBUG << "Endcap" << toString(&trackStateAtCaloEndcap) << endmsg;
    const auto*  tsEP       = trackStateAtCaloEndcap.getReferencePoint();
    const double radSquared = (tsEP[0] * tsEP[0] + tsEP[1] * tsEP[1]);
    if (radSquared < m_minimalTrackStateRadiusSquared) {
      m_log << MSG::DEBUG << "new track state is below tolerance radius" << endmsg;
      return;
    }
    //for curling tracks the propagated track has the wrong z0 whereas it should be 0. really
    if (std::abs(trackStateAtCaloEndcap.getZ0()) >
        std::abs(2. * M_PI / trackStateAtCaloEndcap.getOmega() * trackStateAtCaloEndcap.getTanLambda())) {
      trackStateAtCaloEndcap.setZ0(0.);
    }
    m_log << MSG::DEBUG << "new track state at endcap accepted" << endmsg;
    pandora::InputTrackState pandoraAtEndcap;
    this->CopyTrackState(&trackStateAtCaloEndcap, pandoraAtEndcap);
    trackParameters.m_trackStates.push_back(pandoraAtEndcap);
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDTrackCreatorBase::CalculateTrackTimeAtCalorimeter(const edm4hep::Track pTrack) const {
  edm4hep::TrackState state = pTrack.getTrackStates(edm4hep::TrackState::AtIP);
  const pandora::Helix            helix(state.phi(), state.d0(), state.z0(), state.omega(),
                                        state.tanLambda(), m_settings.m_bField);
  const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());

  // First project to endcap
  float minGenericTime(std::numeric_limits<float>::max());

  pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
  const int                signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
  (void)helix.GetPointInZ(static_cast<float>(signPz) * m_settings.m_eCalEndCapInnerZ, referencePoint,
                          bestECalProjection, minGenericTime);

  // Then project to barrel surface(s)
  pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
  if (m_settings.m_eCalBarrelInnerSymmetry > 0) {
    // Polygon
    float twopi_n = 2. * M_PI / (static_cast<float>(m_settings.m_eCalBarrelInnerSymmetry));

    for (int i = 0; i < m_settings.m_eCalBarrelInnerSymmetry; ++i) {
      float       genericTime(std::numeric_limits<float>::max());
      const float phi(twopi_n * static_cast<float>(i) + m_settings.m_eCalBarrelInnerPhi0);

      const pandora::StatusCode statusCode(helix.GetPointInXY(
          m_settings.m_eCalBarrelInnerR * std::cos(phi), m_settings.m_eCalBarrelInnerR * std::sin(phi),
          std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

      if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime)) {
        minGenericTime     = genericTime;
        bestECalProjection = barrelProjection;
      }
    }
  } else {
    // Cylinder
    float                     genericTime(std::numeric_limits<float>::max());
    const pandora::StatusCode statusCode(
        helix.GetPointOnCircle(m_settings.m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

    if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime)) {
      minGenericTime     = genericTime;
      bestECalProjection = barrelProjection;
    }
  }

  if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  return minGenericTime;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::CopyTrackState(const edm4hep::TrackState* pTrackState,
                                        pandora::InputTrackState& inputTrackState) const {
  if (!pTrackState)
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState->omega()));

  const double px(pt * std::cos(pTrackState->phi()));
  const double py(pt * std::sin(pTrackState->phi()));
  const double pz(pt * pTrackState->tanLambda());

  const double xs(pTrackState->referencePoint().x);
  const double ys(pTrackState->referencePoint().y);
  const double zs(pTrackState->referencePoint().z);

  inputTrackState = pandora::TrackState(xs, ys, zs, px, py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorBase::Settings::Settings()
    : m_trackCollections(StringVector()),
      m_kinkVertexCollections(StringVector()),
      m_prongVertexCollections(StringVector()),
      m_splitVertexCollections(StringVector()),
      m_v0VertexCollections(StringVector()),
      m_prongSplitVertexCollections(StringVector()),
      m_shouldFormTrackRelationships(1),
      m_minTrackHits(5),
      m_minFtdTrackHits(0),
      m_maxTrackHits(5000.f),
      m_d0TrackCut(50.f),
      m_z0TrackCut(50.f),
      m_usingNonVertexTracks(1),
      m_usingUnmatchedNonVertexTracks(0),
      m_usingUnmatchedVertexTracks(1),
      m_unmatchedVertexTrackMaxEnergy(5.f),
      m_d0UnmatchedVertexTrackCut(5.f),
      m_z0UnmatchedVertexTrackCut(5.f),
      m_zCutForNonVertexTracks(250.f),
      m_reachesECalNBarrelTrackerHits(11),
      m_reachesECalNFtdHits(4),
      m_reachesECalBarrelTrackerOuterDistance(-100.f),
      m_reachesECalMinFtdLayer(9),
      m_reachesECalBarrelTrackerZMaxDistance(-50.f),
      m_reachesECalFtdZMaxDistance(1.f),
      m_curvatureToMomentumFactor(0.3f / 2000.f),
      m_minTrackECalDistanceFromIp(100.f),
      m_maxTrackSigmaPOverP(0.15f),
      m_minMomentumForTrackHitChecks(1.f),
      m_maxBarrelTrackerInnerRDistance(50.f),
      m_minBarrelTrackerHitFractionOfExpected(0.2f),
      m_minFtdHitsForBarrelTrackerHitFraction(2),
      m_trackStateTolerance(0.f),
      m_trackingSystemName("DDKalTest"),
      m_bField(0.f),
      m_eCalBarrelInnerSymmetry(0),
      m_eCalBarrelInnerPhi0(0.f),
      m_eCalBarrelInnerR(0.f),
      m_eCalEndCapInnerZ(0.f) {}
