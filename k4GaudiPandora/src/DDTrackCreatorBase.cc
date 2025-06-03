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

#include "DDTrackCreatorBase.h"
#include "Pandora/PdgTable.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

// From Pandora LCContent
#include <LCObjects/LCTrack.h>

#include <algorithm>
#include <cmath>
#include <limits>

// forward declaration
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorBase::DDTrackCreatorBase(const Settings& settings, const pandora::Pandora* const pPandora)
    : m_settings(settings), m_pandora(*pPandora), m_trackVector(0), m_v0TrackList(TrackList()),
      m_parentTrackList(TrackList()), m_daughterTrackList(TrackList()), m_trackToPidMap(TrackToPidMap()),
      m_minimalTrackStateRadiusSquared(0.f) {
  const float ecalInnerR = settings.m_eCalBarrelInnerR;
  const float tsTolerance = settings.m_trackStateTolerance;
  m_minimalTrackStateRadiusSquared = (ecalInnerR - tsTolerance) * (ecalInnerR - tsTolerance);
  // wrap in shared_ptr with a dummy destructor
  //  m_trackingSystem = std::shared_ptr<MarlinTrk::IMarlinTrkSystem>(
  //      MarlinTrk::Factory::createMarlinTrkSystem(settings.m_trackingSystemName, nullptr, ""),
  //      [](MarlinTrk::IMarlinTrkSystem*) {});
  //  m_trackingSystem->init();
  //  FIXME: get info from metadata, collection, or service
  m_encoder = dd4hep::DDSegmentation::BitFieldCoder("foo:2");
  // m_lcTrackFactory = std::make_shared<lc_content::LCTrackFactory>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorBase::~DDTrackCreatorBase() {}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::CreateTrackAssociations(const edm4hep::VertexCollection& kinkCollection,
                                                                const edm4hep::VertexCollection& prongsCollection,
                                                                const edm4hep::VertexCollection& v0Collection) {
  // FIXME: this should be optional for all collections
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(kinkCollection));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(prongsCollection));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(v0Collection));

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractKinks(const edm4hep::VertexCollection& kinkCollection) {
  for (const auto&& pVertex : kinkCollection) {

    // FIXME: this needs the link collection?  We need the particle that represents the incoming particle of this
    // vertex,
    //  not the particle that comes out of this vertex, which is what getParticles would give me?
    const auto& pReconstructedParticle = edm4hep::ReconstructedParticle(); // pVertex.getAssociactedParticles();
    const auto& trackVec = TrackVector();                                  // pReconstructedParticle.getTracks();

    if (this->IsConflictingRelationship(trackVec))
      continue;

    const int vertexPdgCode(pReconstructedParticle.getPDG());

    // Extract the kink vertex information
    size_t iTrack = -1;
    for (auto const& pTrack : trackVec) {
      iTrack++;
      (iTrack == 0) ? m_parentTrackList.insert(pTrack) : m_daughterTrackList.insert(pTrack);
      // streamlog_out(DEBUG) << "KinkTrack " << iTrack << ", nHits " << pTrack->getTrackerHits().size()
      //                      << std::endl;

      int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

      if (iTrack == 0) {
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
          (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
          break;
        }
      }

      m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

      if (0 == m_settings.m_shouldFormTrackRelationships) {
        continue;
      }

      // Make track parent-daughter relationships
      if (iTrack == 0) {
        for (size_t jTrack = iTrack + 1; jTrack < trackVec.size(); ++jTrack) {
          PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                   PandoraApi::SetTrackParentDaughterRelationship(m_pandora, &pTrack, // (FIXME)
                                                                                  &trackVec[jTrack]));
        }
      } else { // Make track sibling relationships
        for (size_t jTrack = iTrack + 1; jTrack < trackVec.size(); ++jTrack) {
          PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                                   PandoraApi::SetTrackSiblingRelationship(m_pandora, &pTrack, // (FIXME)
                                                                           &trackVec[jTrack]));
        }
      }
    }
  }
  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractProngsAndSplits(const edm4hep::VertexCollection& prongsCollection) {
  for (const auto&& pVertex : prongsCollection) {

    // FIXME: this needs the link collection
    const auto& pReconstructedParticle = edm4hep::ReconstructedParticle(); // pVertex.getAssociactedParticles();
    const auto& trackVec = TrackVector();                                  // pReconstructedParticle.getTracks();

    if (this->IsConflictingRelationship(trackVec))
      continue;

    // Extract the prong/split vertex information
    const size_t nTracks = trackVec.size();
    size_t iTrack = -1;
    for (auto const& pTrack : trackVec) {
      iTrack++;
      (0 == iTrack) ? m_parentTrackList.insert(pTrack) : m_daughterTrackList.insert(pTrack);
      // streamlog_out(DEBUG) << "Prong or Split Track " << iTrack << ", nHits " << pTrack->getTrackerHits().size()
      //                      << std::endl;

      if (0 == m_settings.m_shouldFormTrackRelationships)
        continue;

      // Make track parent-daughter relationships
      if (0 == iTrack) {
        for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
          PANDORA_RETURN_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::SetTrackParentDaughterRelationship(m_pandora, &pTrack, &trackVec[jTrack])); // fixme
        }
      }

      // Make track sibling relationships
      else {
        for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
          PANDORA_RETURN_RESULT_IF(
              pandora::STATUS_CODE_SUCCESS, !=,
              PandoraApi::SetTrackSiblingRelationship(m_pandora, &pTrack, &trackVec[jTrack])); // fixme
        }
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDTrackCreatorBase::ExtractV0s(const edm4hep::VertexCollection& v0Collection) {
  for (const auto&& pVertex : v0Collection) {

    // FIXME: this needs the link collection
    const auto& pReconstructedParticle = edm4hep::ReconstructedParticle(); // pVertex.getAssociactedParticles();
    const auto& trackVec = TrackVector();                                  // pReconstructedParticle.getTracks();

    if (this->IsConflictingRelationship(trackVec))
      continue;

    // Extract the v0 vertex information
    const int vertexPdgCode(pReconstructedParticle.getPDG());
    const size_t nTracks = trackVec.size();
    size_t iTrack = -1;
    for (auto const& pTrack : trackVec) {
      m_v0TrackList.insert(pTrack);
      // streamlog_out(DEBUG) << "V0Track " << iTrack << ", nHits " << pTrack->getTrackerHits().size() << std::endl;

      int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

      switch (vertexPdgCode) {
      case pandora::PHOTON:
        (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::E_PLUS : trackPdgCode = pandora::E_MINUS;
        break;
      case pandora::LAMBDA:
        (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::PROTON : trackPdgCode = pandora::PI_MINUS;
        break;
      case pandora::LAMBDA_BAR:
        (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PROTON_BAR;
        break;
      case pandora::K_SHORT:
        (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
        break;
      default:
        (pTrack.getTrackStates()[0].omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
        break;
      }

      m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

      if (0 == m_settings.m_shouldFormTrackRelationships)
        continue;

      // Make track sibling relationships
      for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack) {
        PANDORA_RETURN_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=,
            PandoraApi::SetTrackSiblingRelationship(m_pandora, &pTrack, &trackVec[jTrack])); // FIXME
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DDTrackCreatorBase::IsConflictingRelationship(TrackVector const& trackVec) const {
  for (unsigned int iTrack = 0, nTracks = trackVec.size(); iTrack < nTracks; ++iTrack) {
    const auto& pTrack = trackVec[iTrack];

    if (this->IsDaughter(pTrack) || this->IsParent(pTrack) || this->IsV0(pTrack))
      return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::GetTrackStates(const edm4hep::Track& pTrack,
                                        PandoraApi::Track::Parameters& trackParameters) const {
  const auto& pTrackState = pTrack.getTrackStates(edm4hep::TrackState::AtIP);

  const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState.omega));
  trackParameters.m_momentumAtDca =
      pandora::CartesianVector(std::cos(pTrackState.phi), std::sin(pTrackState.phi), pTrackState.tanLambda) * pt;

  this->CopyTrackState(pTrack.getTrackStates(edm4hep::TrackState::AtFirstHit), trackParameters.m_trackStateAtStart);

  // fg: curling TPC tracks have pointers to track segments stored -> need to get track states from last segment!
  const auto& pEndTrack = (pTrack.getTracks().empty() ? pTrack : *(pTrack.getTracks().end()--));

  this->CopyTrackState(pEndTrack.getTrackStates(edm4hep::TrackState::AtLastHit), trackParameters.m_trackStateAtEnd);
  this->CopyTrackState(pEndTrack.getTrackStates(edm4hep::TrackState::AtCalorimeter),
                       trackParameters.m_trackStateAtCalorimeter);

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

void DDTrackCreatorBase::GetTrackStatesAtCalo(edm4hep::Track const& track,
                                              lc_content::LCTrackParameters& trackParameters) {
  if (not trackParameters.m_reachesCalorimeter.Get()) {
    // streamlog_out(DEBUG5) << "Track does not reach the ECal" << std::endl;
    return;
  }

  const auto& trackAtCalo = track.getTrackStates(edm4hep::TrackState::AtCalorimeter);
  // FIXME: how to check trackAtCalo exists?
  // if (not trackAtCalo) {
  //   // streamlog_out(DEBUG5) << "Track does not have a trackState at calorimeter" << std::endl;
  //   // streamlog_out(DEBUG3) << toString(track) << std::endl;
  //   return;
  // }

  // streamlog_out(DEBUG3) << "Original" << toString(trackAtCalo) << std::endl;

  const auto& tsPosition = trackAtCalo.referencePoint;

  if (std::fabs(tsPosition[2]) < getTrackingRegionExtent()[2]) {
    // streamlog_out(DEBUG5) << "Original trackState is at Barrel" << std::endl;
    pandora::InputTrackState pandoraTrackState;
    this->CopyTrackState(trackAtCalo, pandoraTrackState);
    trackParameters.m_trackStates.push_back(pandoraTrackState);
  } else { // if track state is in endcap we do not repeat track state calculation, because the barrel cannot be hit
    // streamlog_out(DEBUG5) << "Original track state is at Endcap" << std::endl;
    pandora::InputTrackState pandoraTrackState;
    this->CopyTrackState(trackAtCalo, pandoraTrackState);
    trackParameters.m_trackStates.push_back(pandoraTrackState);
    return;
  }

  // FIXME: use marlintrk ported
  //  auto                        marlintrk  =
  //  std::unique_ptr<MarlinTrk::IMarlinTrack>(m_trackingSystem->createTrack()); const EVENT::TrackerHitVec& trkHits =
  //  track->getTrackerHits(); const int                   nHitsTrack = trkHits.size();

  // for (int iHit = 0; iHit < nHitsTrack; ++iHit) {
  //   EVENT::TrackerHit* trkHit = trkHits[iHit];
  //   if (UTIL::BitSet32(
  //           trkHit->getType())[UTIL::ILDTrkHitTypeBit::COMPOSITE_SPACEPOINT]) {  //it is a composite spacepoint
  //     //Split it up and add both hits to the MarlinTrk
  //     const EVENT::LCObjectVec& rawObjects = trkHit->getRawHits();
  //     for (unsigned k = 0; k < rawObjects.size(); k++) {
  //       EVENT::TrackerHit* rawHit = static_cast<EVENT::TrackerHit*>(rawObjects[k]);
  //       if (marlintrk->addHit(rawHit) != MarlinTrk::IMarlinTrack::success) {
  //         //streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add strip hit " << *rawHit
  //                               << std::endl;
  //       }
  //     }
  //   } else {
  //     if (marlintrk->addHit(trkHits[iHit]) != MarlinTrk::IMarlinTrack::success)
  //       //streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to add tracker hit " << *trkHit
  //                             << std::endl;
  //   }
  // }
  // ENDFIXME: use marlintrk ported

  bool tanL_is_positive = trackAtCalo.tanLambda > 0;

  // FIXME: use ported marlintrk
  //  auto trackState = TrackStateImpl(*trackAtCalo);
  int return_error = -1; // FIXME
  // int return_error = marlintrk->initialise(trackState, m_settings.m_bField, MarlinTrk::IMarlinTrack::modeForward);
  // if (return_error != MarlinTrk::IMarlinTrack::success) {
  //   //streamlog_out(DEBUG4) << "DDTrackCreatorBase::GetTrackStatesAtCalo failed to initialize track for endcap track
  //   : "
  //                         << std::endl;
  //   return;
  // }

  double chi2 = -DBL_MAX;
  int ndf = 0;

  edm4hep::TrackState trackStateAtCaloEndcap;

  unsigned ecal_endcap_face_ID = 999; // FIXMElcio::ILDDetID::ECAL_ENDCAP;
  int detElementID = 0;
  // FIXME m_encoder.reset();  // reset to 0
  // FIXME: use DD4hep subdet IDs?
  // (*m_encoder)[lcio::LCTrackerCellID::subdet()] = ecal_endcap_face_ID;
  // (*m_encoder)[lcio::LCTrackerCellID::side()]   = tanL_is_positive ? lcio::ILDDetID::fwd : lcio::ILDDetID::bwd;
  // (*m_encoder)[lcio::LCTrackerCellID::layer()]  = 0;

  // FIXME use ported marlintrk
  //  return_error = marlintrk->propagateToLayer(m_encoder->lowWord(), trackStateAtCaloEndcap, chi2, ndf, detElementID,
  //                                             MarlinTrk::IMarlinTrack::modeForward);
  // streamlog_out(DEBUG5) << "Found trackState at endcap? Error code: " << return_error << std::endl;

  // FIXME  if (return_error == MarlinTrk::IMarlinTrack::success) {
  if (true) {
    // streamlog_out(DEBUG3) << "Endcap" << toString(&trackStateAtCaloEndcap) << std::endl;
    const auto& tsEP = trackStateAtCaloEndcap.referencePoint;
    const double radSquared = (tsEP[0] * tsEP[0] + tsEP[1] * tsEP[1]);
    if (radSquared < m_minimalTrackStateRadiusSquared) {
      // streamlog_out(DEBUG5) << "new track state is below tolerance radius" << std::endl;
      return;
    }
    // for curling tracks the propagated track has the wrong z0 whereas it should be 0. really
    if (std::abs(trackStateAtCaloEndcap.Z0) >
        std::abs(2. * M_PI / trackStateAtCaloEndcap.omega * trackStateAtCaloEndcap.tanLambda)) {
      trackStateAtCaloEndcap.Z0 = 0.0;
    }
    // Streamlog_out(DEBUG5) << "new track state at endcap accepted" << std::endl;
    pandora::InputTrackState pandoraAtEndcap;
    this->CopyTrackState(trackStateAtCaloEndcap, pandoraAtEndcap);
    trackParameters.m_trackStates.push_back(pandoraAtEndcap);
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DDTrackCreatorBase::CalculateTrackTimeAtCalorimeter(const edm4hep::Track& track) const {
  // APS: get the trackstate at IP?
  auto const& ts = track.getTrackStates(edm4hep::TrackState::AtIP);
  const pandora::Helix helix(ts.phi, ts.D0, ts.Z0, ts.omega, ts.tanLambda, m_settings.m_bField);
  const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());

  // First project to endcap
  float minGenericTime(std::numeric_limits<float>::max());

  pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
  const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
  (void)helix.GetPointInZ(static_cast<float>(signPz) * m_settings.m_eCalEndCapInnerZ, referencePoint,
                          bestECalProjection, minGenericTime);

  // Then project to barrel surface(s)
  pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
  if (m_settings.m_eCalBarrelInnerSymmetry > 0) {
    // Polygon
    float twopi_n = 2. * M_PI / (static_cast<float>(m_settings.m_eCalBarrelInnerSymmetry));

    for (int i = 0; i < m_settings.m_eCalBarrelInnerSymmetry; ++i) {
      float genericTime(std::numeric_limits<float>::max());
      const float phi(twopi_n * static_cast<float>(i) + m_settings.m_eCalBarrelInnerPhi0);

      const pandora::StatusCode statusCode(helix.GetPointInXY(
          m_settings.m_eCalBarrelInnerR * std::cos(phi), m_settings.m_eCalBarrelInnerR * std::sin(phi),
          std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

      if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime)) {
        minGenericTime = genericTime;
        bestECalProjection = barrelProjection;
      }
    }
  } else {
    // Cylinder
    float genericTime(std::numeric_limits<float>::max());
    const pandora::StatusCode statusCode(
        helix.GetPointOnCircle(m_settings.m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

    if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime)) {
      minGenericTime = genericTime;
      bestECalProjection = barrelProjection;
    }
  }

  if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  return minGenericTime;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorBase::CopyTrackState(edm4hep::TrackState const& pTrackState,
                                        pandora::InputTrackState& inputTrackState) const {
  // if (!pTrackState)
  //   throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

  const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(pTrackState.omega));

  const double px(pt * std::cos(pTrackState.phi));
  const double py(pt * std::sin(pTrackState.phi));
  const double pz(pt * pTrackState.tanLambda);

  const double xs(pTrackState.referencePoint[0]);
  const double ys(pTrackState.referencePoint[1]);
  const double zs(pTrackState.referencePoint[2]);

  inputTrackState = pandora::TrackState(xs, ys, zs, px, py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDTrackCreatorBase::Settings::Settings()
    : m_trackCollections(StringVector()), m_kinkVertexCollections(StringVector()),
      m_prongVertexCollections(StringVector()), m_splitVertexCollections(StringVector()),
      m_v0VertexCollections(StringVector()), m_prongSplitVertexCollections(StringVector()),
      m_shouldFormTrackRelationships(1), m_minTrackHits(5), m_minFtdTrackHits(0), m_maxTrackHits(5000.f),
      m_d0TrackCut(50.f), m_z0TrackCut(50.f), m_usingNonVertexTracks(1), m_usingUnmatchedNonVertexTracks(0),
      m_usingUnmatchedVertexTracks(1), m_unmatchedVertexTrackMaxEnergy(5.f), m_d0UnmatchedVertexTrackCut(5.f),
      m_z0UnmatchedVertexTrackCut(5.f), m_zCutForNonVertexTracks(250.f), m_reachesECalNBarrelTrackerHits(11),
      m_reachesECalNFtdHits(4), m_reachesECalBarrelTrackerOuterDistance(-100.f), m_reachesECalMinFtdLayer(9),
      m_reachesECalBarrelTrackerZMaxDistance(-50.f), m_reachesECalFtdZMaxDistance(1.f),
      m_curvatureToMomentumFactor(0.3f / 2000.f), m_minTrackECalDistanceFromIp(100.f), m_maxTrackSigmaPOverP(0.15f),
      m_minMomentumForTrackHitChecks(1.f), m_maxBarrelTrackerInnerRDistance(50.f),
      m_minBarrelTrackerHitFractionOfExpected(0.2f), m_minFtdHitsForBarrelTrackerHitFraction(2),
      m_trackStateTolerance(0.f), m_trackingSystemName("DDKalTest"), m_bField(0.f), m_eCalBarrelInnerSymmetry(0),
      m_eCalBarrelInnerPhi0(0.f), m_eCalBarrelInnerR(0.f), m_eCalEndCapInnerZ(0.f) {}
