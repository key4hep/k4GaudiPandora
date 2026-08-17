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

import os

iosvc = IOSvc()
iosvc.Input = "output_pandora_ttbar.root"
iosvc.Output = "output_externalClustering.root"

# GeoSvc used to be configured as a side effect of importing runPandora.py; now that only the
# parameters are imported it has to be set up explicitly here.
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

import sys

# Get the absolute path of the directory where this script resides
# to avoid redefining the pandora input parameters
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

# import the parameters only -- importing runPandora.py itself would re-execute its IOSvc and
# GeoSvc configuration and silently override the input/output files set above
import pandoraParams

# copy, so that the overrides below do not mutate the shared module-level dictionary
params = dict(pandoraParams.params)

# dummy Pandora settings file containing only external clustering algo
params["PandoraSettingsXmlFile"] = current_dir+"/PandoraSettingsExternalClustering.xml"
# rename output collection names to avoid conflicts
params["ClusterCollections"] = ["PandoraClusters"]
params["ClusterCollectionName"] = ["GaudiExternalClusters"]
params["PFOCollectionName"] = ["GaudiExternalClusteringPandoraPFOs"]
params["StartVertexCollectionName"] = ["GaudiExternalClusteringPandoraStartVertices"]

pandora = DDPandoraPFANewAlgorithm("PandoraPFANewAlgorithm", **params)

ApplicationMgr(
    TopAlg=[pandora],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)
