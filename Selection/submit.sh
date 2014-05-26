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

conf=PhaseII/Configuration4v2

while read line #loop over lines in ${conf_file}
do
  array=($line)
  if [ "${array[0]}" != "#" ]; then 
      filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/group/phys_higgs/upgrade/${conf}/140PileUp/${array[0]}/`)
      for ((i=0; i<${#filelist[@]}/100+1; i++)) 
      do
	  outarray=()
	  for ((j=0; j<100; j++))
	  do 
	      if [[ ${filelist[100*$i+$j]} == *root ]]; then
		  outarray+=(${filelist[100*$i+$j]})
	      fi
	  done
	  #echo "./run.sh ${root_script} ${conf} ${array[0]} ${array[1]} ${output_loc} ${outarray[@]}"
	  #./run.sh ${root_script} ${conf} ${array[0]} ${array[1]} ${output_loc} ${outarray[@]}
	  echo " bsub -q 8nh -W 120 -J ${array[0]}_$i run.sh ${root_script} ${conf} ${array[0]} ${array[1]} ${output_loc} ${outarray[@]}"
	  bsub -q 8nh -W 120 -J ${array[0]}_$i run.sh ${root_script} ${conf} ${array[0]} ${array[1]} ${output_loc} ${outarray[@]}
      done

      #for file in $filelist
      #do 
	  #if [[ "${file}" == *root* ]]; then # skip text files in eos
	      #echo $file
	      #echo "   bsub -q 8nh -W 120 -J $file run.sh ${root_script} ${conf} ${array[0]} ${file} ${array[1]} ${output_loc}" 
	      #bsub -q 8nh -W 120 -J $file run.sh ${root_script} ${conf} ${array[0]} ${file} ${array[1]} ${output_loc} # submit to lxbtch
	  #fi
      #done
  fi
done < ${conf_file}