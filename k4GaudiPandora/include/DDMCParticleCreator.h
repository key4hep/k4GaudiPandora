/**
 *  @file   DDMarlinPandora/include/DDMCParticleCreator.h
 * 
 *  @brief  Header file for the mc particle creator class.
 * 
 *  $Log: $
 */

#ifndef DDMCPARTICLECREATOR_H
#define DDMCPARTICLECREATOR_H 1

#include "edm4hep/MCParticle.h"
#include "Api/PandoraApi.h"

#include "DDCaloHitCreator.h"
#include "DDTrackCreatorBase.h"
/**
 *  @brief  DDMCParticleCreator class
 */
class CollectionMaps;
class DDMCParticleCreator
{
public:
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        StringVector    m_mcParticleCollections;                ///< The mc particle collections
        StringVector    m_lcCaloHitRelationCollections;         ///< The SimCaloHit to CaloHit particle relations
        StringVector    m_lcTrackRelationCollections;           ///< The SimTrackerHit to TrackerHit particle relations
        float           m_bField;                               ///< m_bField
    };

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the creator settings
     *  @param  pPandora address of the relevant pandora instance
     */
     DDMCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora);

    /**
     *  @brief  Destructor
     */
     ~DDMCParticleCreator();

    /**
     *  @brief  Create MCParticles
     * 
     *  @param  pLCEvent the lcio event
     */    
    pandora::StatusCode CreateMCParticles(const CollectionMaps& collectionMaps ) const;

    /**
     *  @brief  Create Track to mc particle relationships
     *
     */
    pandora::StatusCode CreateTrackToMCParticleRelationships(const CollectionMaps& collectionMaps, const TrackVector &trackVector) const;

    /**
     *  @brief  Create calo hit to mc particle relationships
     *
     */
    pandora::StatusCode CreateCaloHitToMCParticleRelationships(const CollectionMaps& collectionMaps, const CalorimeterHitVector &calorimeterHitVector) const;

private:
    const Settings          m_settings;                         ///< The mc particle creator settings
    const pandora::Pandora &m_pandora;                          ///< Reference to the pandora object to create the mc particles
    const float             m_bField;                           ///< The bfield
    std::map<unsigned int, const edm4hep::MCParticle*>*  m_id_pMC_map;
};

inline void MCParticleCreator::Reset()
{
    m_id_pMC_map->clear();
}


#endif // #ifndef DDMCPARTICLECREATOR_H
