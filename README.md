# GRAiNITA EM calorimeter in ALFA detector

This is an inclusive repo for Grainita ECAL simulation. 

It contains the geometry description with dd4hep (Grainita\_ECAL), specific sensitive detector action (GrainitaCaloSDAction) and necessary segmentation (detectorSegmentations). 

Welcome to the DD4hep Tutorials for DRD6 repository!

To build the framework: 
```bash
# Login to an Alma9 machine with CVMFS mounted, e.g. lxplus
source /cvmfs/sw.hsf.org/key4hep/setup.sh
git clone git@github.com:qwert2333/DD4hep_Grainita.git
cd DD4hep_Grainita
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -Wno-dev
make install -j 8
cd ..
k4_local_repo # update the local environment

# Next time when login to lxplus: 
source setup.sh

# To run an example of Grainita full simulation: 
cd run
ddsim --steeringFile fullsim_steering.py

```


