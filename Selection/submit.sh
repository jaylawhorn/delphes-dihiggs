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

root_script=$1
  conf_file=$2
 output_loc=$3

while read line #loop over lines in ${conf_file}
do
  array=($line)
  #for conf in PhaseI/Configuration0 PhaseII/Configuration3 PhaseII/Configuration4v2 # to process all three configurations
  for conf in PhaseII/Configuration4v2 # or just one configuration
  do
    if [ "${array[0]}" != "#" ]; then 
	# get list of files in eos for that sample+configuration combination
	filelist=`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/group/phys_higgs/upgrade/${conf}/140PileUp/${array[0]}/`
	for file in $filelist
	do 
	  if [[ "${file}" == *root* ]]; then # skip text files in eos
	      echo "   bsub -q 8nh -W 120 -J $file run.sh ${root_script} ${conf} ${array[0]} ${file} ${array[1]} ${output_loc}" 
	      bsub -q 8nh -W 120 -J $file run.sh ${root_script} ${conf} ${array[0]} ${file} ${array[1]} ${output_loc} # submit to lxbtch
	      #./run.sh ${root_script} ${conf} ${array[0]} ${file} ${array[1]} ${output_loc} # run locally
	  fi
	done
    fi
  done
done < ${conf_file}