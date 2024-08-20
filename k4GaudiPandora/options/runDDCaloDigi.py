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
from Configurables import DDCaloDigi

from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
import os

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"]+"/FCCee/CLD/compact/CLD_o2_v06/CLD_o2_v06.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

calodigi = DDCaloDigi()

#set properties
calodigi.InputColIsECAL = True                              # True -- ECAL // False -- HCAL
calodigi.InputCaloHitCollection = ["ECalBarrelCollection"]  # "ECalBarrelCollection","ECalEndcapCollection"
                                                            # "HCalBarrelCollection","HCalEndcapCollection","HCalRingCollection"
calodigi.OutputCaloHitCollection = ["ECALBarrel"]
calodigi.RelCollection = ["RelationCaloHit"]

# digitazing parameters for ECAL and HCAL
calodigi.ECALThreshold = 5.0e-5
calodigi.ECALThresholdUnit = "GeV"
calodigi.HCALThreshold = [0.00025]
calodigi.HCALThresholdUnit = "GeV"
calodigi.ECALLayers = [41,100]
calodigi.HCALLayers = [100]
calodigi.CalibrECAL = [37.5227197175, 37.5227197175]
calodigi.CalibrHCALBarrel = [45.9956826061]
calodigi.CalibrHCALEndcap = [46.9252540291]
calodigi.CalibrHCALOther = [57.4588011802]
calodigi.IfDigitalEcal = 0
calodigi.MapsEcalCorrection = 0
calodigi.IfDigitalHcal = 0
calodigi.ECALGapCorrection = 1
calodigi.ECALGapCorrection = 1
calodigi.ECALEndcapCorrectionFactor = 1.03245503522
calodigi.HCALEndcapCorrectionFactor = 1.000
calodigi.ECALGapCorrectionFactor = 1.0
calodigi.ECALModuleGapCorrectionFactor = 0.0
calodigi.HCALModuleGapCorrectionFactor = 0.5

calodigi.Histograms = 0

# timing parameters for ECAL
calodigi.UseEcalTiming = 1
calodigi.ECALCorrectTimesForPropagation = 1
calodigi.ECALTimeWindowMin = -1.0
calodigi.ECALEndcapTimeWindowMax = 10.0
calodigi.ECALBarrelTimeWindowMax = 10.0
calodigi.ECALDeltaTimeHitResolution = 10.0
calodigi.ECALTimeResolution = 10.0
calodigi.ECALSimpleTimingCut = True

# timing parameters for HCAL
calodigi.UseHcalTiming = 1
calodigi.HCALCorrectTimesForPropagation = 1
calodigi.HCALTimeWindowMin = -1.0
calodigi.HCALEndcapTimeWindowMax = 10.0
calodigi.HCALBarrelTimeWindowMax = 10.0
calodigi.HCALDeltaTimeHitResolution = 10.0
calodigi.HCALTimeResolution = 10.0
calodigi.HCALSimpleTimingCut = True

# parameters for extra ECAL digitization effects
calodigi.CalibECALMIP = 1.0e-4
calodigi.ECALApplyRealisticDigi = 0
calodigi.ECAL_PPD_PE_per_MIP = 7.0
calodigi.ECAL_PPD_N_Pixels = 10000
calodigi.ECAL_PPD_N_Pixels_uncertainty = 0.05
calodigi.ECAL_miscalibration_uncorrel = 0.0
calodigi.ECAL_miscalibration_uncorrel_memorise = False
calodigi.ECAL_miscalibration_correl = 0.0
calodigi.ECAL_deadCellRate = 0.0
calodigi.ECAL_deadCell_memorise = False
calodigi.ECAL_strip_absorbtionLength = 1.0e6
calodigi.ECAL_pixel_spread = 0.05
calodigi.ECAL_elec_noise_mips = 0.0
calodigi.energyPerEHpair = 3.6
calodigi.ECAL_maxDynamicRange_MIP = 2500.0
calodigi.StripEcal_default_nVirtualCells = 9
calodigi.ECAL_default_layerConfig = "000000000000000"

# parameters for extra HCAL digitization effects
calodigi.CalibHCALMIP = 1.0e-4
calodigi.HCALApplyRealisticDigi = 0
calodigi.HCAL_PPD_PE_per_MIP = 10.0
calodigi.HCAL_PPD_N_Pixels = 400
calodigi.HCAL_PPD_N_Pixels_uncertainty = 0.05
calodigi.HCAL_miscalibration_uncorrel = 0.0
calodigi.HCAL_miscalibration_uncorrel_memorise = False
calodigi.HCAL_miscalibration_correl = 0.0
calodigi.HCAL_deadCellRate = 0.0
calodigi.HCAL_deadCell_memorise = False
calodigi.HCAL_pixel_spread = 0.0
calodigi.HCAL_elec_noise_mips = 0.0
calodigi.HCAL_maxDynamicRange_MIP = 200.0


iosvc = IOSvc()
iosvc.input = "../simulation/sim_edm4hep.root"
iosvc.output = "../outputfiles/DDCaloDigi/outputCaloDigi_Gaudi.root"

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "../outputfiles/DDCaloDigi/ddcalodigi_hist.root"

ApplicationMgr(TopAlg=[calodigi],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
               OutputLevel=INFO,
               )
