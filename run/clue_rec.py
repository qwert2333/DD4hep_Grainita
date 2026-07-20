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
from Gaudi.Configuration import WARNING, INFO, DEBUG, VERBOSE
from Configurables import ClueGaudiAlgorithmWrapper2D, CLUENtuplizer, THistSvc, EventDataSvc, MetadataSvc, GeoSvc
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input = ["Digi/Digi_photon_2GeV_theta80-100_phi80-100.root"]
iosvc.Output = "Rec/ClueReco_photon_2GeV_theta80-100_phi80-100_dc015_dm025_rho1.root"

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["/afs/cern.ch/user/f/faguo/eos/FCCSW/DD4hep_Grainita/Grainita_ECAL/compact/Grainita_ECAL_Barrel_v02_Long4.xml"]
geoservice.OutputLevel = INFO

dc = 0.015
rho = 0.1
dm = 0.025

MyClueGaudiAlgorithmWrapper = ClueGaudiAlgorithmWrapper2D("ClueGaudiAlgorithmWrapperName",
    CaloHitsCollections = ["GrainitaEcalBarrelDigiHit"],
    CriticalDistance = dc,
    MinLocalDensity = rho,
    FollowerDistance = dm,
    OutputLevel = INFO,
    strategy = "MergeCollections", # "PerDetectorRegion", "PerCollection" , "MergeCollections"
    coordinate = "Polar", # "Cartesian"
    SaveClustersAsHits = True,
    CLUEHitCollName = "CLUECalorimeterHitCollection"
)


MyCLUENtuplizer = CLUENtuplizer("CLUEAnalysis",
    CLUEHitCollName = "CLUECalorimeterHitCollection",
    ReadOutName = "GrainitaEcalBarrelRO",
    LayerFieldName = "rho",
    OutputLevel = INFO
)

iosvc.outputCommands = [
    "drop *",
    "keep CLUEClusters",
]

filename = "rec DATAFILE='Rec/ClueNtuple_photon_2GeV_dc015_dm025_rho1.root' TYP='ROOT' OPT='RECREATE'"
THistSvc().Output = [filename]
THistSvc().OutputLevel = WARNING
THistSvc().PrintAll = False
THistSvc().AutoSave = True
THistSvc().AutoFlush = True

mgr = ApplicationMgr( TopAlg = [MyClueGaudiAlgorithmWrapper, MyCLUENtuplizer],
                  EvtSel = 'NONE',
                  EvtMax   = -1,
                  ExtSvc = [EventDataSvc("EventDataSvc"), MetadataSvc("MetadataSvc"), geoservice],
                  OutputLevel=INFO
                )
