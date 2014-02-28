#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Execute one job (works interactively and when executed in condor)
#---------------------------------------------------------------------------------------------------

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
      ${1}+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${2}/140PileUp/${3}/${4}\",${5},\"/afs/cern.ch/work/k/klawhorn/SnowmassSamples/${2}/${4}\"\)

  root -b -l -q rootlogon.C \
      ${1}+\(\"root://eoscms.cern.ch//store/group/phys_higgs/upgrade/${2}/140PileUp/${3}/${4}\",${5},\"/afs/cern.ch/work/k/klawhorn/SnowmassSamples/${2}/${4}\"\)

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

exit $status
