/**
 *  @file   DDMarlinPandora/include/DDCaloHitCreatorALLEGRO.h
 *
 *  @brief  Header file for the calo hit creator class.
 *
 *  $Log: $
 */

#ifndef DDCALO_HIT_CREATORALLEGRO_H
#define DDCALO_HIT_CREATORALLEGRO_H 1

#include "Api/PandoraApi.h"

#include "DDCaloHitCreator.h"

/**
 *  @brief  DDCaloHitCreator class
 */
class DDCaloHitCreatorALLEGRO : public DDCaloHitCreator {
public:
  /**
   *  @brief  Constructor
   *
   *  @param  settings the creator settings
   *  @param  pPandora address of the relevant pandora instance
   */
  DDCaloHitCreatorALLEGRO(const Settings& settings, pandora::Pandora& pandora, const Gaudi::Algorithm* algorithm);
  ~DDCaloHitCreatorALLEGRO() override = default;

  /**
   *  @brief  Create calo hits
   *
   *  @param  pLCEvent the lcio event
   */
  virtual pandora::StatusCode
  createCaloHits(const std::map<std::string, std::vector<edm4hep::CalorimeterHit>>& eCaloHitsMap,
                 const std::vector<edm4hep::CalorimeterHit>& hCalCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& muonCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& lCalCaloHits,
                 const std::vector<edm4hep::CalorimeterHit>& lhCalCaloHits) const override;

private:
  /**
   *  @brief  Get common calo hit properties: position, parent address, input energy and time
   *
   *  @param  pCaloHit the lcio calorimeter hit
   *  @param  caloHitParameters the calo hit parameters to populate
   */
  void getCommonCaloHitProperties(const edm4hep::CalorimeterHit& hit,
                                  PandoraApi::CaloHit::Parameters& caloHitParameters) const;
};
#endif
