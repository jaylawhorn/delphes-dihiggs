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

 root_script=$1
delphes_conf=$2
 sample_name=$3
   file_name=$4
   cross_sec=$5
  output_loc=$6

cd /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

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
echo \
  root -b -l -q rootlogon.C \
      ${root_script}+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${delphes_conf}/140PileUp/${sample_name}/${file_name}\",${cross_sec},\"${output_loc}/${delphes_conf}/${file_name}\"\)

  root -b -l -q rootlogon.C \
      ${root_script}+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${delphes_conf}/140PileUp/${sample_name}/${file_name}\",${cross_sec},\"${output_loc}/${delphes_conf}/${file_name}\"\)

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

exit $status
