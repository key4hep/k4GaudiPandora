/**
 *  @file   DDMarlinPandora/src/DDTrackCreatorALLEGRO.cc
 *
 *  @brief  Implementation of the track creator class for a ALLEGRO all silicon tracker.
 *
 *  $Log: $
 */

#include "DDTrackCreatorALLEGRO.h"
#include "DDTrackCreatorBase.h"

#include <LCObjects/LCTrack.h>
#include <Pandora/PdgTable.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/Detector.h>
#include <DD4hep/DetectorSelector.h>
#include <DDRec/DetectorData.h>

#include <edm4hep/TrackCollection.h>

#include <algorithm>
#include <cmath>
#include <limits>

// forward declarations. See in DDPandoraPFANewProcessor.cc
std::vector<double> getTrackingRegionExtent();

DDTrackCreatorALLEGRO::DDTrackCreatorALLEGRO(const Settings& settings, pandora::Pandora& pandora,
                                             const Gaudi::Algorithm* thisAlg)
    : DDTrackCreatorBase(settings, pandora, thisAlg), m_trackerInnerR(0.f), m_trackerOuterR(0.f), m_trackerZmax(0.f),
      m_cosTracker(0.f), m_endcapDiskInnerRadii(DoubleVector()), m_endcapDiskOuterRadii(DoubleVector()),
      m_endcapDiskZPositions(DoubleVector()), m_nEndcapDiskLayers(0), m_barrelTrackerLayers(0),
      m_tanLambdaEndcapDisk(0.f)

{
  // FIXME! AD: currently ignoring the tracker parameters since we are working with truth tracks...
  /*
  m_trackerInnerR = getTrackingRegionExtent()[0];
  m_trackerOuterR = getTrackingRegionExtent()[1];
  m_trackerZmax = getTrackingRegionExtent()[2];

  ///FIXME: Probably need to be something relating to last disk inner radius
  m_cosTracker = m_trackerZmax / std::sqrt(m_trackerZmax * m_trackerZmax + m_trackerInnerR * m_trackerInnerR);

  dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();

  //Maybe we need to veto the vertex? That was done in the ILD case
  const std::vector< dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors(  (
  dd4hep::DetType::TRACKER | dd4hep::DetType::BARREL )) ;

  m_barrelTrackerLayers = 0;
  for (std::vector< dd4hep::DetElement>::const_iterator iter = barrelDets.begin(), iterEnd = barrelDets.end();iter !=
  iterEnd; ++iter){ try
      {
          dd4hep::rec::ZPlanarData * theExtension = 0;

          const dd4hep::DetElement& theDetector = *iter;
          theExtension = theDetector.extension<dd4hep::rec::ZPlanarData>();

          unsigned int N = theExtension->layers.size();
          m_barrelTrackerLayers=m_barrelTrackerLayers+N;

          streamlog_out( DEBUG2 ) << " Adding layers for barrel tracker from DD4hep for "<< theDetector.name()<< "- n
  layers: " << N<< " sum up to now: "<<m_barrelTrackerLayers<<endmsg; } catch (std::runtime_error &exception){
          m_algorithm.warning() << "DDTrackCreatorALLEGRO exception during Barrel Tracker layer sum for
  "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << endmsg;
      }
  }

  m_nEndcapDiskLayers=0;
  m_endcapDiskInnerRadii.clear();
  m_endcapDiskOuterRadii.clear();

  //Instead of gear, loop over a provided list of forward (read: endcap) tracking detectors. For ILD this would be FTD
   const std::vector< dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors(  (
  dd4hep::DetType::TRACKER | dd4hep::DetType::ENDCAP )) ;

   for (std::vector< dd4hep::DetElement>::const_iterator iter = endcapDets.begin(), iterEnd = endcapDets.end();iter !=
  iterEnd; ++iter){ try
      {
          dd4hep::rec::ZDiskPetalsData * theExtension = 0;

          const dd4hep::DetElement& theDetector = *iter;
          theExtension = theDetector.extension<dd4hep::rec::ZDiskPetalsData>();

          unsigned int N = theExtension->layers.size();
          streamlog_out( DEBUG2 ) << " Filling FTD-like parameters from DD4hep for "<< theDetector.name()<< "- n layers:
  " << N<< endmsg;

          for(unsigned int i = 0; i < N; ++i)
          {
              dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer  = theExtension->layers[i];

              // Create a disk to represent even number petals front side
              //FIXME! VERIFY THAT TIS MAKES SENSE!
              m_endcapDiskInnerRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm);
              m_endcapDiskOuterRadii.push_back(thisLayer.distanceSensitive/dd4hep::mm+thisLayer.lengthSensitive/dd4hep::mm);

              // Take the mean z position of the staggered petals
              const double zpos(thisLayer.zPosition/dd4hep::mm);
              m_endcapDiskZPositions.push_back(zpos);
              streamlog_out( DEBUG2 ) << "     layer " << i << " - mean z position = " << zpos << endmsg;
          }

          m_nEndcapDiskLayers = m_endcapDiskZPositions.size() ;
      } catch (std::runtime_error &exception){
          m_algorithm.warning() << "DDTrackCreatorALLEGRO exception during Forward Tracking Disk parameter construction
  for detector "<<const_cast<dd4hep::DetElement&>(*iter).name()<<" : " << exception.what() << endmsg;
      }
  }

  for (unsigned int iEndcapDiskLayer = 0; iEndcapDiskLayer < m_nEndcapDiskLayers; ++iEndcapDiskLayer)
  {
      if ((std::fabs(m_endcapDiskOuterRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon()) ||
          (std::fabs(m_endcapDiskInnerRadii[iEndcapDiskLayer]) < std::numeric_limits<float>::epsilon()))
      {
          throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
      }
  }

  m_tanLambdaEndcapDisk = m_endcapDiskZPositions[0] / m_endcapDiskOuterRadii[0];
  */
}

pandora::StatusCode DDTrackCreatorALLEGRO::CreateTracks(const std::vector<edm4hep::Track>& tracks) {
  for (const auto& pTrack : tracks) {
    m_algorithm.debug()
        << " Warning! Ignoring expected number of hits and other hit number cuts. Should eventually change!" << endmsg;

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
    auto trackPIDiter = m_trackToPidMap.find(GetTrackID(pTrack));

    if (trackPIDiter != m_trackToPidMap.end()) {
      trackParameters.m_particleId = trackPIDiter->second;
      trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackPIDiter->second);
    }

    if (0.f != signedCurvature)
      trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

    try {
      this->GetTrackStates(pTrack, trackParameters);
      this->TrackReachesECAL(pTrack, trackParameters);
      // FIXME! AD: this is disabled for ALLEGRO since I am not sure what to do here with truth tracks
      // this->GetTrackStatesAtCalo(pTrack, trackParameters);
      this->DefineTrackPfoUsage(pTrack, trackParameters);

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

bool DDTrackCreatorALLEGRO::PassesQualityCuts(const edm4hep::Track& pTrack,
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

  // Check momentum uncertainty is reasonable to use track
  const pandora::CartesianVector& momentumAtDca(trackParameters.m_momentumAtDca.Get());
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

  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::DefineTrackPfoUsage(const edm4hep::Track& pTrack,
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
      const pandora::CartesianVector& momentumAtDca(trackParameters.m_momentumAtDca.Get());
      const float pZ = momentumAtDca.GetZ();
      const float pT = std::hypot(momentumAtDca.GetX(), momentumAtDca.GetY());

      const float zCutForNonVertexTracks = m_trackerInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks;
      const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) &&
                                   (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance));

      const bool isV0 = IsV0(pTrack);
      const bool isDaughter = IsDaughter(pTrack);

      // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
      // if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_trackerInnerR +
      // m_settings.m_maxBarrelTrackerInnerRDistance))
      // FIXME! AD: (rInner < m_trackerInnerR + m_settings.m_maxBarrelTrackerInnerRDistance) check is removed
      if ((firstTrackState.D0 < m_settings.m_d0TrackCut) && (firstTrackState.Z0 < m_settings.m_z0TrackCut)) {
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
  }

  trackParameters.m_canFormPfo = canFormPfo;
  trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DDTrackCreatorALLEGRO::TrackReachesECAL(const edm4hep::Track& pTrack,
                                             PandoraApi::Track::Parameters& trackParameters) const {
  // FIXME! AD: since currently we have do not have track hits -> check the radius of the trackState AtCalo.
  // currently check is done only for ECAL barrel but should check if track reaches the Endcap!
  pandora::CartesianVector posAtCalo(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition());
  double radiusAtCalo = std::sqrt(posAtCalo.GetX() * posAtCalo.GetX() + posAtCalo.GetY() * posAtCalo.GetY());
  if (radiusAtCalo >= m_settings.m_eCalBarrelInnerR)
    trackParameters.m_reachesCalorimeter = true;
  return;

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

    if ((r > m_trackerInnerR) && (r < m_trackerOuterR) && (std::fabs(z) <= m_trackerZmax)) {
      nTrackerHits++;
      continue;
    }

    for (unsigned int j = 0; j < m_nEndcapDiskLayers; ++j) {
      if ((r > m_endcapDiskInnerRadii[j]) && (r < m_endcapDiskOuterRadii[j]) &&
          (std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_endcapDiskZPositions[j]) &&
          (std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_endcapDiskZPositions[j])) {
        if (static_cast<int>(j) > maxOccupiedEndcapDiskLayer)
          maxOccupiedEndcapDiskLayer = static_cast<int>(j);

        nEndcapDiskHits++;
        break;
      }
    }
  }

  // Require sufficient hits in barrel or endcap trackers, then compare extremal hit positions with tracker dimensions
  if ((nTrackerHits >= m_settings.m_reachesECalNBarrelTrackerHits) ||
      (nEndcapDiskHits >= m_settings.m_reachesECalNFtdHits)) {
    if ((hitOuterR - m_trackerOuterR > m_settings.m_reachesECalBarrelTrackerOuterDistance) ||
        (std::fabs(hitZMax) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
        (std::fabs(hitZMin) - m_trackerZmax > m_settings.m_reachesECalBarrelTrackerZMaxDistance) ||
        (maxOccupiedEndcapDiskLayer >= m_settings.m_reachesECalMinFtdLayer)) {
      trackParameters.m_reachesCalorimeter = true;
      return;
    }
  }

  // If track is lowpt, it may curl up and end inside tpc inner radius
  const pandora::CartesianVector& momentumAtDca(trackParameters.m_momentumAtDca.Get());
  const float cosAngleAtDca = std::fabs(momentumAtDca.GetZ()) / momentumAtDca.GetMagnitude();
  const float pT = std::hypot(momentumAtDca.GetX(), momentumAtDca.GetY());

  if ((cosAngleAtDca > m_cosTracker) ||
      (pT < m_settings.m_curvatureToMomentumFactor * m_settings.m_bField * m_trackerOuterR)) {
    trackParameters.m_reachesCalorimeter = true;
    return;
  }

  trackParameters.m_reachesCalorimeter = false;
}
