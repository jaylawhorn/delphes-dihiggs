#!/bin/bash

#file_list=`/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls /store/group/phys_higgs/upgrade/${conf}/140PileUp/${SAMPLE}/`

while read line
do

  array=($line)
  #for conf in PhaseI/Configuration0 PhaseII/Configuration3 PhaseII/Configuration4v2
  for conf in PhaseII/Configuration4v2
  #for conf in PhaseII/Configuration3
  do

    if [ "${array[0]}" != "" ]; then 

	filelist=`/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select ls /store/group/phys_higgs/upgrade/${conf}/140PileUp/${array[0]}/`

	for file in $filelist
	do 

	  if [[ "${file}" == *root* ]]; then
	      echo "   bsub -q 8nh -W 120 -J $file run_cluster.sh ${1} ${conf} ${array[0]} ${file} ${array[1]} "
	      bsub -q 8nh -W 120 -J $file run_cluster.sh ${1} ${conf} ${array[0]} ${file} ${array[1]}
	  fi
	done
    fi
  done
done < $2