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
from Configurables import PandoraPFAIdeaAlgorithm

import os

iosvc = IOSvc()
iosvc.Input = "input_reco.root"
iosvc.Output = "output_pandora.root"

# Whether the geometry has a dual-readout HCAL endcap.  This defaults to False, i.e. to the
# reduced barrel-only geometry that the CI runs; set it True for production.  The algorithm
# warns loudly when it is False.
hasHcalEndcap = False

# detector geometry
# if K4GEO is empty, this should use relative path to working directory
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/IDEA/compact/IDEA_o2_v01/IDEA_o2_v01.xml'
    if hasHcalEndcap
    else 'FCCee/IDEA/compact/IDEA_o2_v01_CI/IDEA_o2_v01_CI.xml'
]

geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]


# per-calorimeter steering: the five vectors below are indexed together, one entry each
caloHitCollections  = ["SCEPCal_digi_cheren", "SCEPCal_digi_scint", "DRBTScin_digi", "DRBTCher_digi"]
caloSystemIDs       = [4, 5, 28]
caloCollectionTypes = ["ECAL", "ECAL", "HCAL"]
caloLayerFieldNames = ["depth", "depth", ""]
caloEncodingStrings = [
    "system:5,phi:7,theta:11,gamma:4,epsilon:4,depth:1,cherenkov:1",
    "system:5,phi:7,theta:11,gamma:4,epsilon:4,depth:1,cherenkov:1",
    "system:5,stave:10,tower:-8,air:6,col:-16,row:16,clad:1,core:1,cherenkov:1",
]
caloCellSizes       = [10, 10, 2.]

if hasHcalEndcap:
    caloHitCollections  += ["DRETScinLeft_digi", "DRETCherLeft_digi",
                            "DRETScinRight_digi", "DRETCherRight_digi"]
    caloSystemIDs       += [25]
    caloCollectionTypes += ["HCAL"]
    caloLayerFieldNames += [""]
    caloEncodingStrings += ["system:5,stave:10,tower:6,air:1,col:16,row:16,clad:1,core:1,cherenkov:1"]
    caloCellSizes       += [2.]

params = {
    "PandoraSettingsXmlFile": "PandoraSettingsIdea.xml",
    "inputTrackCollection": "TracksFromGenParticles",
    "inputClusterCollections": ["TopoGrownClusters"],
    "outputPfoCollection": "PandoraPfaIdea",
    "outputClusterCollection": "PandoraClusters",
    "CherenkovFieldName": "cherenkov",
    "HasHcalEndcap": hasHcalEndcap,
    "inputCaloHitCollections": caloHitCollections,
    "CaloSystemIDs": caloSystemIDs,
    "CaloCollectionTypes": caloCollectionTypes,
    "CaloLayerFieldNames": caloLayerFieldNames,
    "CaloEncodingStrings": caloEncodingStrings,
    "CaloCellSizes": caloCellSizes,
}

pandoraIdea = PandoraPFAIdeaAlgorithm("PandoraPFAIdeaAlgorithm", **params)

ApplicationMgr(
    TopAlg=[pandoraIdea],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), geoservice],
    OutputLevel=INFO,
)
