#ifndef DDSimpleMuonDigi_H
#define DDSimpleMuonDigi_H

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCRecoCaloAssociationCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include "CalorimeterHitType.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"
#include "k4FWCore/Transformer.h"
#include "DDRec/SurfaceManager.h"

#include <random>
#include <string>
#include <vector>

struct DDSimpleMuonDigi final
  : k4FWCore::MultiTransformer<
  std::tuple<edm4hep::CalorimeterHitCollection,edm4hep::MCRecoCaloAssociationCollection>(
											  const edm4hep::SimCalorimeterHitCollection&, const edm4hep::EventHeaderCollection&)> {
  DDSimpleMuonDigi(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::MCRecoCaloAssociationCollection> operator()
  (const edm4hep::SimCalorimeterHitCollection& simCaloHits,
   const edm4hep::EventHeaderCollection& headers) const override;

private:
  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", "VXD", "Name of the subdetector"};
  Gaudi::Property<std::vector<int>> m_layersToKeepBarrelVec{
    this, "KeepBarrelLayersVec", {0}, "Vector of Barrel layers to be kept. Layers start at 1!"};
  Gaudi::Property<std::vector<int>> m_layersToKeepEndCapVec{
    this, "KeepEndcapLayersVec", {0}, "Vector of Endcap layers to be kept. Layers start at 1!"};
  //Gaudi::Property<std::vector<bool>> useLayersBarrelVec{this, "useBarrelLayerVector", false, "whether to use the endcap layer vector"};
  //Gaudi::Property<std::vector<bool>> useLayersEndcapVec{this, "useEndCapLayerVector", false, "whether to use the EndCap layer vector"};
  Gaudi::Property<std::vector<std::string>> m_muonCollections{this, "muon_collections", {"muon_collections"}, "Collection of Muons"};
  Gaudi::Property<std::string> outputRelCollection{this, "outputRelCollection", "outputRelCollection",
    "The output collection of relations"};
  Gaudi::Property<std::string> outputMuonCollection{this, "outputMuonCollection","outputMuonCollection",
    "The output collection of muons"};
  Gaudi::Property<std::string> m_encodingStringVariable{
      this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
      "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<std::string> m_cellIDLayerString {this, "CellIDLayerString","Layer", "Name of the part of the cellID that holds the layer"};                                           
  Gaudi::Property<float> m_thresholdMuon{this, "MuonThreshold", {0.025}, "Threshold for muon"};
  Gaudi::Property<float> m_timeThresholdMuon{this, "timethresholdMuon", {0.025}, "time threshold for muons"};
  Gaudi::Property<float> m_calibrCoeffMuon{this, "calibrationCoeffmuon", {120000.0}, "Callibration coefficient of muons"};
  Gaudi::Property<float> m_maxHitEnergyMuon{this, "maxMuonHitEnergy", {2.0}, "Threshold for maximum muon hit energy"};
  Gaudi::Property<std::string> m_detectorNameBarrel{this, "detectornameB", "YokeBarrel", "Name of the subdetector"};
  Gaudi::Property<std::string> m_detectorNameEndcap{this, "detectornameE", "YokeEndcap",
    "Name of the second subdetector"};
  
  std::string   m_collName;
  std::vector<bool>  m_useLayersBarrelVec{}, m_useLayersEndcapVec{};
  SmartIF<IGeoSvc>         m_geoSvc;
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) ;
  float computeHitTime( const edm4hep:: SimCalorimeterHit &h ) const ;

};
DECLARE_COMPONENT(DDSimpleMuonDigi)
#endif

// namespace EVENT {
//   class SimCalorimeterHit ;
// }

// /** === DDSimpleMuonDigi Processor === <br>
//  *  Simple calorimeter digitizer for the muon detectors.
//  *  Converts SimCalorimeterHit collections to one
//  *  CalorimeterHit collection applying a threshold and an calibration constant...
//  *
//  *  @version $Id$
//  */
// class DDSimpleMuonDigi : public Processor {

//  public:

//   virtual Processor*  newProcessor() { return new DDSimpleMuonDigi ; }

//   DDSimpleMuonDigi() ;

//   virtual void init() ;

//   virtual void processRunHeader( LCRunHeader* run ) ;

//   virtual void processEvent( LCEvent * evt ) ;

//   virtual void check( LCEvent * evt ) ;

//   virtual void end() ;

//   bool useLayer(CHT::Layout caloLayout, unsigned int layer) ;
//   float computeHitTime( const EVENT::SimCalorimeterHit *h ) const ;

//  protected:

//   int _nRun = 0;
//   int _nEvt = 0;

//   IntVec _layersToKeepBarrelVec{}, _layersToKeepEndcapVec{};
//   std::vector<bool>  _useLayersBarrelVec{}, _useLayersEndcapVec{};

//   std::vector<std::string> _muonCollections{};

//   std::string _outputMuonCollection = "";
//   std::string _outputRelCollection = "";

//   std::string _cellIDLayerString = "layer";

//   float _thresholdMuon = 0.025;
//   float _timeThresholdMuon = _thresholdMuon ;
//   float _calibrCoeffMuon = 120000;
//   float _maxHitEnergyMuon = 2.0;

//   std::string _detectorNameBarrel = "YokeBarrel";
//   std::string _detectorNameEndcap = "YokeEndcap";

// } ;

// #endif
