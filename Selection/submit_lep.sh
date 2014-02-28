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

root_script=oldSelect.C
output_loc=/afs/cern.ch/work/k/klawhorn/resTests/

# get list of files in eos for that sample+configuration combination
filelist=`ls /afs/cern.ch/work/k/klawhorn/public/HHToBBTT_14TeV`
for file in $filelist
  do 
  if [[ "${file}" == *root* ]]; then # skip text files in eos
      echo "   bsub -q 8nh -W 120 -J $file run_lep.sh ${root_script} ${file} 2.92 ${output_loc}" 
      bsub -q 8nh -W 120 -J $file run_lep.sh ${root_script} ${file} 2.92 ${output_loc}
  fi
done