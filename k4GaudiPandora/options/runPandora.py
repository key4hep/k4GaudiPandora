#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import DDPandoraPFANewAlgorithm

from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
# from Configurables import RootHistSvc
# from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
import os

iosvc = IOSvc()
iosvc.Input = "output_REC.edm4hep.root"
iosvc.Output = "output_pandora.root"

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

params = {
    "FinalEnergyDensityBin": 110.0,
    "MaxClusterEnergyToApplySoftComp": 200.0,
    "TrackCollections": ["SiTracks_Refitted"],
    "ECalCaloHitCollections": ["ECALBarrel", "ECALEndcap"],
    "HCalCaloHitCollections": ["HCALBarrel", "HCALEndcap", "HCALOther"],
    "LCalCaloHitCollections": [],
    "LHCalCaloHitCollections": [],
    "MuonCaloHitCollections": ["MUON"],
    "MCParticleCollections": ["MCParticles"],
    "RelCaloHitCollections": ["RelationCaloHit", "RelationMuonHit"],
    # "RelTrackCollections": ["SiTracks_Refitted_Relation"],
    "RelTrackCollections": [],
    "KinkVertexCollections": [],
    "ProngVertexCollections": [],
    "SplitVertexCollections": [],
    "V0VertexCollections": [],
    "ClusterCollectionName": ["GaudiPandoraClusters"],
    "PFOCollectionName": ["GaudiPandoraPFOs"],
    "CreateGaps": False,
    "MinBarrelTrackerHitFractionOfExpected": 0,
    "MinFtdHitsForBarrelTrackerHitFraction": 0,
    "MinFtdTrackHits": 0,
    "MinMomentumForTrackHitChecks": 0,
    "MinTrackECalDistanceFromIp": 0,
    "MinTrackHits": 0,
    "ReachesECalBarrelTrackerOuterDistance": -100,
    "ReachesECalBarrelTrackerZMaxDistance": -50,
    "ReachesECalFtdZMaxDistance": 1,
    "ReachesECalMinFtdLayer": 0,
    "ReachesECalNBarrelTrackerHits": 0,
    "ReachesECalNFtdHits": 0,
    "UnmatchedVertexTrackMaxEnergy": 5,
    "UseNonVertexTracks": 1,
    "UseUnmatchedNonVertexTracks": 0,
    "UseUnmatchedVertexTracks": 1,
    "Z0TrackCut": 200,
    "Z0UnmatchedVertexTrackCut": 5,
    "ZCutForNonVertexTracks": 250,
    "MaxTrackHits": 5000,
    "MaxTrackSigmaPOverP": 0.15,
    "CurvatureToMomentumFactor": 0.00015,
    "D0TrackCut": 200,
    "D0UnmatchedVertexTrackCut": 5,
    "StartVertexAlgorithmName": "PandoraPFANew",
    "StartVertexCollectionName": ["GaudiPandoraStartVertices"],
    "YokeBarrelNormalVector": [0, 0, 1],
    "HCalBarrelNormalVector": [0, 0, 1],
    "ECalBarrelNormalVector": [0, 0, 1],
    "MuonBarrelBField": -1.0,
    "MuonEndCapBField": 0.01,
    "EMConstantTerm": 0.01,
    "EMStochasticTerm": 0.17,
    "HadConstantTerm": 0.03,
    "HadStochasticTerm": 0.6,
    "InputEnergyCorrectionPoints": [],
    "LayersFromEdgeMaxRearDistance": 250,
    "NOuterSamplingLayers": 3,
    "TrackStateTolerance": 0,
    "MaxBarrelTrackerInnerRDistance": 200,
    "MinCleanCorrectedHitEnergy": 0.1,
    "MinCleanHitEnergy": 0.5,
    "MinCleanHitEnergyFraction": 0.01,
    "MuonHitEnergy": 0.5,
    "ShouldFormTrackRelationships": 1,
    "TrackCreatorName": "DDTrackCreatorCLIC",
    "TrackSystemName": "DDKalTest",
    "OutputEnergyCorrectionPoints": [],
    "UseEcalScLayers": 0,
    "ECalScMipThreshold": 0,
    "ECalScToEMGeVCalibration": 1,
    "ECalScToHadGeVCalibrationBarrel": 1,
    "ECalScToHadGeVCalibrationEndCap": 1,
    "ECalScToMipCalibration": 1,
    "ECalSiMipThreshold": 0,
    "ECalSiToEMGeVCalibration": 1,
    "ECalSiToHadGeVCalibrationBarrel": 1,
    "ECalSiToHadGeVCalibrationEndCap": 1,
    "ECalSiToMipCalibration": 1,
    "StripSplittingOn": 0,
    # Settings for CalorimeterIntegrationTimeWindow = 10 ns
    "PandoraSettingsXmlFile": "PandoraSettingsCLD/PandoraSettingsDefault.xml",
    "SoftwareCompensationWeights": [
        2.40821,
        -0.0515852,
        0.000711414,
        -0.0254891,
        -0.0121505,
        -1.63084e-05,
        0.062149,
        0.0690735,
        -0.223064,
    ],
    "ECalToMipCalibration": "175.439",
    "HCalToMipCalibration": "45.6621",
    "ECalMipThreshold": "0.5",
    "HCalMipThreshold": "0.3",
    "ECalToEMGeVCalibration": "1.01776966108",
    "HCalToEMGeVCalibration": "1.01776966108",
    "ECalToHadGeVCalibrationBarrel": "1.11490774181",
    "ECalToHadGeVCalibrationEndCap": "1.11490774181",
    "HCalToHadGeVCalibration": "1.00565042407",
    "MuonToMipCalibration": "20703.9",
    "DigitalMuonHits": "0",
    "MaxHCalHitHadronicEnergy": "10000000.",
}


pandora = DDPandoraPFANewAlgorithm("PandoraPFANewProcessor", **params)

ApplicationMgr(
    TopAlg=[pandora],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
