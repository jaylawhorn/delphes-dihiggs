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

count=0

while read line #loop over lines in ${conf_file}
do
  file=$line
  echo "   bsub -q 8nh small_run.sh select_all ${file} ${count}"
  bsub -q 8nh small_run.sh select_all ${file} ${count}
  count=$[count+1]
done < big_ttbar_list.txt