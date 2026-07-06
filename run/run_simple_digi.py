import os

from Gaudi.Configuration import *

# Loading the input SIM file, defining output file
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = os.environ.get(
    "DIGI_INPUT", 
    "CaloSim_gamma_5GeV_theta85-95_phi80-100.root"
)
io_svc.Output = os.environ.get("DIGI_OUTPUT", "CaloDigi_gamma_5GeV_theta85-95_phi80-100.root")


################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("dd4hep_grainita", "")
print(path_to_detector)
detectors_to_use=[ '../Grainita_ECAL/compact/Grainita_ECAL_Barrel_v02_Long4.xml' ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO


################ Dual-readout calorimeter
# SiPM emulation
from Configurables import SimpleSiPMDigiAlg
simpledigi = SimpleSiPMDigiAlg("SimpleSiPMDigiAlg")
simpledigi.InputSimCaloHitCollection = "GrainitaCalorimeterHits"
simpledigi.OutputCaloHitCollection = "GrainitaEcalBarrelDigiHit"
simpledigi.OutputCaloSimLinkCollection = "GrainitaEcalBarrelDigiHit_SimHit_link"
simpledigi.OutputCaloMCPLinkCollection = "GrainitaEcalBarrelDigiHit_MCParticle_link"
simpledigi.ReadOutName = "GrainitaEcalBarrelRO"
simpledigi.UseDigi = True
simpledigi.LightYield = 10  # ph / MeV
simpledigi.OutFileName = "DigiTuple_gamma_5GeV.root"
simpledigi.OutputLevel = INFO


from Configurables import SimpleCalibAlg
simplecali = SimpleCalibAlg("SimpleCalibAlg");
simplecali.InputCaloHitCollection = "GrainitaEcalBarrelDigiHit"
simplecali.OutputCaloHitCollection = "GrainitaEcalBarrelCalibHit"
simplecali.ReadOutName = "GrainitaEcalBarrelRO"
simplecali.InputFormat = "Energy"
simplecali.CalibrationConstant = 1.
simplecali.ApplyRhoPitchCorrection = True
simplecali.AttLength = 0.34   # unit in cm
simplecali.OutputLevel = WARNING

# RNG for sipm emulation (TODO harmonize RNG with other modules)
from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine
rndmEngine = RndmEngine('RndmGenSvc.Engine',
  SetSingleton = True,
  Seeds = [ int(os.environ.get("DIGI_SEED", "1234567")) ] # default seed is 1234567
)

#from Configurables import CreateTruthLinks
#createTruthLinks = CreateTruthLinks("CreateTruthLinks",
#    cell_hit_links=["GrainitaCaloSiPMreadoutDigiHit_link"],
#    clusters=["TopoClusterAll"],
#    mcparticles="MCParticles",
#    cell_mcparticle_links="CaloHitMCParticleLinks",
#    cluster_mcparticle_links="ClusterMCParticleLinks",
#    OutputLevel=INFO
#)

from Configurables import RndmGenSvc
rndmGenSvc = RndmGenSvc("RndmGenSvc",
  Engine = rndmEngine.name()
)

################ Output
io_svc.outputCommands = [
  "keep *",
]

# Profiling
from Configurables import AuditorSvc, ChronoAuditor, UniqueIDGenSvc
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

from k4FWCore import ApplicationMgr
application_mgr = ApplicationMgr(
    TopAlg = [
        simpledigi,
        simplecali
    ],
    EvtSel = 'NONE',
    EvtMax = 100,
    ExtSvc = [
        EventDataSvc("EventDataSvc"),
        geoservice,
        audsvc,
        UniqueIDGenSvc("uidSvc"),
        rndmEngine,
        rndmGenSvc
    ],
    StopOnSignal = True,
)

for algo in application_mgr.TopAlg:
    algo.AuditExecute = True

