#----------------------------------------------------
# script to merge and clean up output from lxbtch submissions
#
# example command
# ./mergeSelFiles.sh LL-4p-0-100-v1510_14TEV /afs/cern.ch/work/k/klawhorn/SnowmassSamples/PhaseII/Configuration4v2 .
#
# Jay Lawhorn 11/4/13
#----------------------------------------------------
#!/bin/bash

  sample_type=$1
file_location=$2
 out_location=$3

count=0

for file in `ls ${file_location}/${sample_type}*`
do
  file_array[$count]=$file
  let "count+=1"
done
 
hadd -f ${out_location}/${sample_type}_temp.root ${file_array[@]}
 
root -l -q cleanUpMergedFiles.C+\(\"${out_location}/${sample_type}_temp.root\",\"${out_location}/${sample_type}.root\"\)
