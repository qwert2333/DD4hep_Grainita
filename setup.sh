# DD4hep
#source /cvmfs/sw.hsf.org/key4hep/setup.sh
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh -r 2026-07-03

# Add the sipm lib
export LD_LIBRARY_PATH=$(echo $CMAKE_PREFIX_PATH | tr ':' '\n' | grep simsipm | head -n1)/lib64:$LD_LIBRARY_PATH


# Configures the environment to ensure the system can locate the DD4hep detector builders:
k4_local_repo


