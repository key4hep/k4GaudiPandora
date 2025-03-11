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
#ifndef DDPFO_CREATOR_H
#define DDPFO_CREATOR_H 1

#include <k4FWCore/Transformer.h>
#include "Api/PandoraApi.h"
#include "ClusterShapes.h"
#include "Objects/CartesianVector.h"

#include <edm4hep/ClusterCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>
#include <edm4hep/VertexCollection.h>
#include <edm4hep/MutableReconstructedParticle.h>
#include <edm4hep/MutableCluster.h>

class CollectionMaps;
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DDPfoCreator class
 */
class DDPfoCreator {
public:
  /**
     *  @brief  Settings class
     */
  class Settings {
  public:
    /**
         *  @brief  Default constructor
         */
    Settings();

    std::string m_startVertexAlgName = "";  ///< The name of the algorithm to fill the start vertex output collection
    float       m_emStochasticTerm   = 0;   ///< The stochastic term for EM shower energy resolution
    float       m_hadStochasticTerm  = 0;   ///< The stochastic term for Hadronic shower energy resolution
    float       m_emConstantTerm     = 0;   ///< The constant term for EM shower energy resolution
    float       m_hadConstantTerm    = 0;   ///< The constant term for Hadronic shower energy resolution
  };

  /**
     *  @brief  Constructor
     *
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
  DDPfoCreator(const Settings& settings, const pandora::Pandora* const pPandora, IMessageSvc* msgSvc);

  /**
     *  @brief  Destructor
     */
  ~DDPfoCreator();

  /**
     *  @brief  Create particle flow objects
     *
     */
  pandora::StatusCode CreateParticleFlowObjects(
      edm4hep::ClusterCollection&               pClusterCollection,
      edm4hep::ReconstructedParticleCollection& pReconstructedParticleCollection,
      edm4hep::VertexCollection&                pStartVertexCollection);


private:
  /**
     *  @brief  index for the subdetector
     */
  enum Index { ECAL_INDEX = 0, HCAL_INDEX = 1, YOKE_INDEX = 2, LCAL_INDEX = 3, LHCAL_INDEX = 4, BCAL_INDEX = 5 };

  /**
     *  @brief  initialise sub detector name strings
     *
     *  @param  subDetectorNames to receive the list of sub detector names
     */
  void InitialiseSubDetectorNames(pandora::StringVector& subDetectorNames) const;

  /**
     *  @brief  Set sub detector energies for a cluster
     *
     *  @param  subDetectorNames the list of sub detector names
     *  @param  pCluster the address of the edm4hep cluster to be set sub detector energies
     *  @param  pandoraCaloHitList the pandora calorimeter hit list
     *  @param  hitE the vector to receive the energy of hits
     *  @param  hitX the vector to receive the x position of hits
     *  @param  hitY the vector to receive the y position of hits
     *  @param  hitZ the vector to receive the z position of hits
     */
  void SetClusterSubDetectorEnergies(const pandora::StringVector& subDetectorNames,
                                     edm4hep::MutableCluster*       pCluster,
                                     const pandora::CaloHitList&  pandoraCaloHitList, pandora::FloatVector& hitE,
                                     pandora::FloatVector& hitX,  pandora::FloatVector& hitY,
                                     pandora::FloatVector& hitZ) const;

  /**
     *  @brief  Set cluster energies and errors
     *
     *  @param  pPandoraPfo the address of the pandora pfo
     *  @param  pPandoraCluster the address of the pandora cluster
     *  @param  pCluster the address of the edm4hep cluster to be set energies and erros
     *  @param  clusterCorrectEnergy a number to receive the cluster correct energy
     */
  void SetClusterEnergyAndError(const pandora::ParticleFlowObject* const pPandoraPfo,
                                const pandora::Cluster* const pPandoraCluster, edm4hep::MutableCluster* pCluster,
                                float& clusterCorrectEnergy) const;

  /**
     *  @brief  Set cluster position, errors and other shape info, by calculating culster shape first
     *
     *  @param  nHitsInCluster number of hits in cluster
     *  @param  hitE the vector of the energy of hits
     *  @param  hitX the vector of the x position of hits
     *  @param  hitY the vector of the y position of hits
     *  @param  hitZ the vector of the z position of hits
     *  @param  pCluster the lcio cluster to be set positions and errors
     *  @param  clusterPosition a CartesianVector to receive the cluster position
     */
  void SetClusterPositionAndError(const unsigned int nHitsInCluster, pandora::FloatVector& hitE,
                                  pandora::FloatVector& hitX, pandora::FloatVector& hitY, pandora::FloatVector& hitZ,
                                  edm4hep::MutableCluster*   pCluster,
                                  pandora::CartesianVector& clusterPositionVec) const;

  /**
     *  @brief  Calculate reference point for pfo with tracks
     *
     *  @param  pPandoraPfo the address of the pandora pfo
     *  @param  referencePoint a CartesianVector to receive the reference point
     */
  pandora::StatusCode CalculateTrackBasedReferencePoint(const pandora::ParticleFlowObject* const pPandoraPfo,
                                                        pandora::CartesianVector&                referencePoint) const;

  /**
     *  @brief  Set reference point of the reconstructed particle
     *
     *  @param  referencePoint a CartesianVector of the reference point
     *  @param  pReconstructedParticle the address of the reconstructed particle to be reference point
     */
  void SetRecoParticleReferencePoint(const pandora::CartesianVector&       referencePoint,
                                     edm4hep::MutableReconstructedParticle* pReconstructedParticle) const;

  /**
     *  @brief  Add tracks to reconstructed particle
     *
     *  @param  pPandoraPfo the address of the pandora pfo
     *  @param  pReconstructedParticle the address of the reconstructed particle to be added tracks
     */
  void AddTracksToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                               edm4hep::MutableReconstructedParticle*    pReconstructedParticle) const;

  /**
     *  @brief  Set properties of reconstructed particle from pandora pfo
     *
     *  @param  pPandoraPfo the address of the pandora pfo
     *  @param  pReconstructedParticle the address of the reconstructed particle to be set properties
     */
  void SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject* const pPandoraPfo,
                                        edm4hep::MutableReconstructedParticle*    pReconstructedParticle) const;

  /**
     *  @brief  Whether parent and daughter tracks are associated with the same pfo
     *
     *  @param  pPandoraTrack the address of the pandora track
     *  @param  allTrackList list of all tracks associated with reconstructed particle
     *
     *  @return boolean
     */
  bool IsValidParentTrack(const pandora::Track* const pPandoraTrack, const pandora::TrackList& allTrackList) const;

  /**
     *  @brief  Whether sibling tracks are associated with the same pfo
     *
     *  @param  pPandoraTrack the address of the pandora track
     *  @param  allTrackList list of all tracks associated with reconstructed particle
     *
     *  @return boolean
     */
  bool HasValidSiblingTrack(const pandora::Track* const pPandoraTrack, const pandora::TrackList& allTrackList) const;

  /**
     *  @brief  Whether the track is the closest (of those associated with the same pfo) to the interaction point
     *
     *  @param  pPandoraTrack the address of the pandora track
     *  @param  allTrackList list of all tracks associated to reconstructed particle
     *
     *  @return boolean
     */
  bool IsClosestTrackToIP(const pandora::Track* const pPandoraTrack, const pandora::TrackList& allTrackList) const;

  /**
     *  @brief  Whether at least one track sibling track is associated to the reconstructed particle
     *
     *  @param  pPandoraTrack the address of the pandora track
     *  @param  allTrackList list of all tracks associated to reconstructed particle
     *
     *  @return boolean
     */
  bool AreAnyOtherSiblingsInList(const pandora::Track* const pPandoraTrack,
                                 const pandora::TrackList&   allTrackList) const;

  const Settings          m_settings;  ///< The pfo creator settings
  const pandora::Pandora& m_pandora;   ///< Reference to the pandora object from which to extract the pfos
  IMessageSvc* m_msgSvc;
};

#endif  // #ifndef DDPFO_CREATOR_H
