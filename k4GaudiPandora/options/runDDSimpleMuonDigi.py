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
from Configurables import DDSimpleMuonDigi
from Configurables import CollectionMerger

from Configurables import GeoSvc
import os

iosvc = IOSvc()
iosvc.Input = "output_REC.edm4hep.root"
iosvc.Output = "output_muon_digi.root"

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

input_collections = ["YokeBarrelCollection", "YokeEndcapCollection"]
output_collections = ["GaudiMuonYokeBarrelCollection", "GaudiMuonYokeEndcapCollection"]
output_relation = ["GaudiRelationMuonYokeBarrelHit", "GaudiRelationMuonYokeEndcapHit"]

digi = [
    DDSimpleMuonDigi("MuonYokeBarrelDigi"),
    DDSimpleMuonDigi("MuonYokeEndcapDigi"),
]

for alg, inputcol, outputcol, outrel in zip(
    digi, input_collections, output_collections, output_relation
):

    alg.SubDetectorName = "VXD"
    alg.KeepBarrelLayersVec = []
    alg.KeepEndcapLayersVec = []

    alg.MUONCollection = [inputcol]
    alg.RelationOutputCollection = [outrel]
    alg.MUONOutputCollection = [outputcol]

    alg.CellIDLayerString = "layer"

    alg.CalibrMUON = 70.1
    alg.MaxHitEnergyMUON = 2.0
    alg.MuonThreshold = 1e-6

# Merge the output collections to make it easy to compare to
# the single one that is obtained by original Marlin processor
merger = CollectionMerger("MuonYokeMerger")
merger.InputCollections = output_collections
merger.OutputCollection = ["GaudiMUON"]

relation_merger = CollectionMerger("MuonRelationMerger")
relation_merger.InputCollections = output_relation
relation_merger.OutputCollection = ["GaudiRelationMUON"]

ApplicationMgr(
    TopAlg=digi + [merger, relation_merger],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
