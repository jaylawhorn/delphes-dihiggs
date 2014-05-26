#!/bin/bash
#------------------------------------------------------------
# Submit a batch of jobs to lxbtch
#
# example command:
#             root_script     conf_file output_location
# ./submit.sh submitDelphes.C xsec.txt  /afs/cern.ch/work/k/klawhorn/SnowmassSamples/
# 
# conf_file has format (no leading "#")
# sample_type cross_section
#
# Note: conf file needs empty last line... clunky I know.
# Jay Lawhorn 11/4/13
#------------------------------------------------------------

cd /afs/cern.ch/user/k/klawhorn/DelphesDiHiggs/Selection

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

in_location=/afs/cern.ch/work/k/klawhorn/test
out_location=/afs/cern.ch/work/k/klawhorn/new_ntuples

while read line #loop over lines in ${conf_file}
do
  array=($line)
  if [ "${array[0]}" != "#" ]; then 
      ./mergeSelFiles.sh ${array[0]} ${in_location} ${out_location}
  fi
done < xsec.txt 