#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Execute one job (works interactively and when executed in lxbtch)
#
# example local command
# ./run.sh selectDelphes.C PhaseII/Configuration4v2 LL-4p-0-100-v1510_14TEV LL-4p-0-100-v1510_14TEV_100005594_PhaseII_Conf4v2_140PileUp.root 1341.36923 /afs/cern.ch/work/k/klawhorn/SnowmassSamples
#
# example lxbtch command
# bsub -q 8nh -W 120 run.sh selectDelphes.C PhaseII/Configuration4v2 LL-4p-0-100-v1510_14TEV LL-4p-0-100-v1510_14TEV_100005594_PhaseII_Conf4v2_140PileUp.root 1341.36923 /afs/cern.ch/work/k/klawhorn/SnowmassSamples
#
# Jay Lawhorn 11/4/13
#---------------------------------------------------------------------------------------------------

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/select_all.C .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/select_all_C.* .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/rootlogon.C .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/big_ttbar_list.txt .

# some basic printing
echo " "; echo "${h}: Show who and where we are";
echo " "
echo " user executing: "`id`;
echo " running on    : "`hostname`;
echo " executing in  : "`pwd`;
echo " submitted from: $HOSTNAME";
echo " ";

# initialize the CMSSW environment
echo " "; echo "${h}: Initialize CMSSW (in $CMSSW_BASE)"; echo " "
workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd -

# get ready to run
echo " "; echo "${h}: Starting root job now"; echo " ";

  root -b -l -q rootlogon.C \
      select_all.C+\(\"big_ttbar_list.txt\"\)

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

exit $status
