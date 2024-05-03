/**
 *  @file   DDMarlinPandora/src/DDMCParticleCreator.cc
 * 
 *  @brief  Implementation of the mc particle creator class.
 * 
 *  $Log: $
 */

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
// Define BaseClass_t
#include "k4FWCore/BaseClass.h"
#include "edm4hep/MCParticle.h" 
#include "edm4hep/MCRecoCaloAssociation.h" 
#include "edm4hep/SimCalorimeterHit.h" 
#include "edm4hep/CaloHitContribution.h" 
#include "edm4hep/Track.h" 
#include "edm4hep/MCRecoTrackerAssociation.h" 
#include "edm4hep/SimTrackerHit.h" 

#include "PandoraPFAlg.h"
#include "DDMCParticleCreator.h"

#include <string>

//forward declarations. See in DDPandoraPFANewProcessor.cc
double getFieldFromCompact();


DDMCParticleCreator::DDMCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pandora(*pPandora),
    m_bField(getFieldFromCompact())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::~DDMCParticleCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateMCParticles(const EVENT::LCEvent *const pLCEvent) const
{
    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCParticleCollection = pLCEvent->getCollection(*iter);

            for (int i = 0, iMax = pMCParticleCollection->getNumberOfElements(); i < iMax; ++i)
            {
                try
                {
                    EVENT::MCParticle *pMcParticle = dynamic_cast<MCParticle*>(pMCParticleCollection->getElementAt(i));

                    if (NULL == pMcParticle)
                        throw EVENT::Exception("Collection type mismatch");

                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy = pMcParticle->getEnergy();
                    mcParticleParameters.m_particleId = pMcParticle->getPDG();
                    mcParticleParameters.m_mcParticleType = pandora::MC_3D;
                    mcParticleParameters.m_pParentAddress = pMcParticle;
                    mcParticleParameters.m_momentum = pandora::CartesianVector(pMcParticle->getMomentum()[0], pMcParticle->getMomentum()[1],
                        pMcParticle->getMomentum()[2]);
                    mcParticleParameters.m_vertex = pandora::CartesianVector(pMcParticle->getVertex()[0], pMcParticle->getVertex()[1],
                        pMcParticle->getVertex()[2]);
                    mcParticleParameters.m_endpoint = pandora::CartesianVector(pMcParticle->getEndpoint()[0], pMcParticle->getEndpoint()[1],
                        pMcParticle->getEndpoint()[2]);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

                    // Create parent-daughter relationships
                    for(MCParticleVec::const_iterator itDaughter = pMcParticle->getDaughters().begin(),
                        itDaughterEnd = pMcParticle->getDaughters().end(); itDaughter != itDaughterEnd; ++itDaughter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora,
                            pMcParticle, *itDaughter));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract MCParticle: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract MCParticles collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateTrackToMCParticleRelationships(const CollectionMaps& collectionMaps,const TrackVector &trackVector) const
{
    for (unsigned ik = 0; ik < trackVector.size(); ik++)
    {
        const edm4hep::Track *pTrack = trackVector.at(ik);
        // Get reconstructed momentum at dca
        const pandora::Helix helixFit(pTrack->getTrackStates(0).phi, pTrack->getTrackStates(0).D0, pTrack->getTrackStates(0).Z0, pTrack->getTrackStates(0).omega, pTrack->getTrackStates(0).tanLambda, m_bField);
        const float recoMomentum(helixFit.GetMomentum().GetMagnitude());

                    // Use momentum magnitude to identify best mc particle
        edm4hep::MCParticle *pBestMCParticle = NULL;
        float bestDeltaMomentum(std::numeric_limits<float>::max());
        try
        {
        
            for (StringVector::const_iterator iter = m_settings.m_TrackRelationCollections.begin(), iterEnd = m_settings.m_TrackRelationCollections.end(); iter != iterEnd; ++iter)
            {       
                if(collectionMaps.collectionMap_TrkRel.find(*iter) == collectionMaps.collectionMap_TrkRel.end()) continue;
                const std::vector<edm4hep::MCRecoTrackerAssociation>& pMCRecoTrackerAssociationCollection = (collectionMaps.collectionMap_TrkRel.find(*iter))->second;
                for(unsigned ith=0 ; ith<pTrack->trackerHits_size(); ith++)
                {
                    for(unsigned ic=0; ic < pMCRecoTrackerAssociationCollection.size(); ic++)
                    {
                        if( pMCRecoTrackerAssociationCollection.at(ic).getRec().id() != pTrack->getTrackerHits(ith).id() ) continue;
                        const edm4hep::ConstSimTrackerHit pSimHit = pMCRecoTrackerAssociationCollection.at(ic).getSim();
                        const edm4hep::ConstMCParticle ipa = pSimHit.getMCParticle();
                        if( m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end() ) continue;
                        const float trueMomentum(pandora::CartesianVector(ipa.getMomentum()[0], ipa.getMomentum()[1], ipa.getMomentum()[2]).GetMagnitude());
                        const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));
                            if (deltaMomentum < bestDeltaMomentum)
                            {
                            pBestMCParticle =const_cast<edm4hep::MCParticle*>((*m_id_pMC_map)[ipa.id()]);
                            bestDeltaMomentum = deltaMomentum;
                             }
                    }
                }
            }

                   
             if (NULL == pBestMCParticle)continue;
             PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack,
                        pBestMCParticle));
        }
        catch (pandora::StatusCodeException &statusCodeException)
        {
        streamlog_out(ERROR) << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
           }
        catch (EVENT::Exception &exception)
         {
            streamlog_out(WARNING) << "Failed to extract track to mc particle relationship: " << exception.what() << std::endl;
         }
            
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract track to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode DDMCParticleCreator::CreateCaloHitToMCParticleRelationships(const CollectionMaps& collectionMaps, const CalorimeterHitVector &calorimeterHitVector) const
{
    typedef std::map<MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(), iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloRel.find(*iter) == collectionMaps.collectionMap_CaloRel.end()) continue;
        try
        {
            const std::vector<edm4hep::MCRecoCaloAssociation>& pMCRecoCaloAssociationCollection = (collectionMaps.collectionMap_CaloRel.find(*iter))->second;
            
            for (unsigned i_calo=0; i_calo < calorimeterHitVector.size(); i_calo++)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();

                    for (unsigned ic=0; ic < pMCRecoCaloAssociationCollection.size(); ic++)
                    {
                        if( pMCRecoCaloAssociationCollection.at(ic).getRec().id() != (*(calorimeterHitVector.at(i_calo))).id() ) continue;
                        const edm4hep::ConstSimCalorimeterHit pSimHit = pMCRecoCaloAssociationCollection.at(ic).getSim();
                        for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont)
                        {
                           edm4hep::ConstCaloHitContribution conb = pSimHit.getContributions(iCont);
                            const edm4hep::ConstMCParticle ipa = conb.getParticle();
                            float  ien = conb.getEnergy();
                            if( m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end() ) continue;
                            const edm4hep::MCParticle * p_tmp = (*m_id_pMC_map)[ipa.id()]; 
                            mcParticleToEnergyWeightMap[p_tmp] += ien;
                        }
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora,
                            *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship: " << exception.what() << std::endl;
                }
            }
        }
        catch (...)
        {
            streamlog_out(DEBUG5) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

DDMCParticleCreator::Settings::Settings()
{
}
