B#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Execute one job (works interactively and when executed in lxbtch)
#1;2c
# example local command
# ./run.sh selectDelphes.C PhaseII/Configuration4v2 LL-4p-0-100-v1510_14TEV LL-4p-0-100-v1510_14TEV_100005594_PhaseII_Conf4v2_140PileUp.root 1341.36923 /afs/cern.ch/work/k/klawhorn/SnowmassSamples
#
# example lxbtch command
# bsub -q 8nh -W 120 run.sh selectDelphes.C PhaseII/Configuration4v2 LL-4p-0-100-v1510_14TEV LL-4p-0-100-v1510_14TEV_100005594_PhaseII_Conf4v2_140PileUp.root 1341.36923 /afs/cern.ch/work/k/klawhorn/SnowmassSamples
#
# Jay Lawhorn 11/4/13
#---------------------------------------------------------------------------------------------------

 input_array=("$@")
 root_script=${input_array[0]}
delphes_conf=${input_array[1]}
 sample_name=${input_array[2]}
   cross_sec=${input_array[3]}
  output_loc=${input_array[4]}

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/${root_script}.C .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/${root_script}_C.* .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/rootlogon.C .
cp /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection/mt2.hh .

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

for ((i=5; i<${#input_array[@]}; i++))
do
echo root -b -l -q rootlogon.C ${root_script}.C+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${delphes_conf}/140PileUp/${sample_name}/${input_array[$i]}\",${cross_sec},\"${output_loc}/${input_array[$i]}\"\)
root -b -l -q rootlogon.C ${root_script}.C+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${delphes_conf}/140PileUp/${sample_name}/${input_array[$i]}\",${cross_sec},\"${output_loc}/${input_array[$i]}\"\)
done

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

exit $status