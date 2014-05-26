#!/bin/bash

#master run script for pre-processing

conf_file=$1
eos_switch=$2

if [ ${eos_switch} -eq "1" ]; then
    base_path="/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/"
elif [ ${eos_switch} -eq "0" ]; then
    base_path="/afs/cern.ch/work/j/jlawhorn/"
else
    echo "EOS toggle doesn't make sense. please try again"
    exit
fi

while read line #loop over lines in ${conf_file}
do
  array=($line)
  if [ "${array[0]:0:1}" != "#" ]; then 
      echo "Preparing sample" ${array[0]}
      rm ${array[0]}.txt
      filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${base_path}${array[0]}/`)
      outarray=()
      for ((i=0; i<${#filelist[@]}; i++)) 
      do
	  if [[ ${filelist[$i]} == *root ]]; then
	      outarray+=(${filelist[$i]})
	  fi
      done
      echo ${#outarray[@]} "files found in" ${base_path}${array[0]}

      for ((i=0; i<${#outarray[@]}; i++))
      do
	  echo "root://eoscms.cern.ch/"${base_path}${array[0]}${outarray[i]} >> ${array[0]}.txt
      done

      echo `wc -l ${array[0]}.txt | cut -d' ' -f1` "files listed in" ${array[0]}.txt

      root -l -q -b getall.C+\(\"${array[0]}.txt\"\)
  fi
done < ${conf_file}

