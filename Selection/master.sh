#!/bin/bash

#master run script for pre-processing

if [ ${#} -lt "2" ]; then
    echo "Please use the following syntax: " 
    echo "    ./master.sh <conf_file> <submit>"
    echo " "
    echo "<conf_file> is a text file containing the names of all samples"
    echo "<submit> is either 0 if you want to submit the jobs by hand or"
    echo "         1 if you want the jobs to be automatically submitted"
    exit
fi

conf_file=$1
submit=$2

while read line #loop over lines in ${conf_file}
do
  array=($line)
  if [ "${array[0]:0:1}" != "#" ]; then 
      echo "Preparing sample" ${array[0]}
      rm ${array[0]}.txt
      if [ ${array[2]} -eq "1" ]; then
	  filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${array[3]} | grep ${array[0]}`)
      elif [ ${array[2]} -eq "0" ]; then
	  filelist=(`ls ${array[3]} | grep ${array[0]}`)
      fi
      outarray=()
      for ((i=0; i<${#filelist[@]}; i++)) 
      do
	  if [[ ${filelist[$i]} == *root ]]; then
	      outarray+=(${filelist[$i]})
	  fi
      done
      echo ${#outarray[@]} "files found in" ${array[3]}

      for ((i=0; i<${#outarray[@]}; i++))
      do
	  if [ ${array[2]} -eq "1" ]; then
	      echo "root://eoscms.cern.ch/"${array[3]}"/"${outarray[i]} >> ${array[0]}.txt
	  elif [ ${array[2]} -eq "0" ]; then
	      echo ${array[3]}"/"${outarray[i]} >> ${array[0]}.txt
	  fi
      done

      echo `wc -l ${array[0]}.txt | cut -d' ' -f1` "files listed in" ${array[0]}.txt

      root -l -q -b getall.C+\(\"${array[0]}.txt\"\) | tail -1 > ${array[0]}_events.txt

      echo `cat ${array[0]}_events.txt` "events found"

      mkdir -p "/afs/cern.ch/work/j/jlawhorn/public/ntuples/"${array[0]}

      n=0

      rm ${array[0]}_run.sh
      head -7 run_outline.sh > ${array[0]}_run.sh
      while read line2
      do
	  echo root -l -q -b selection.C+\\\(\\\"${line2}\\\",${array[1]},`cat ${array[0]}_events.txt`,\\\""/afs/cern.ch/work/j/jlawhorn/public/ntuples/"${array[0]}/${n}".root"\\\"\\\) >> ${array[0]}_run.sh
	  n=$((n+1))
	      
      done < ${array[0]}.txt

      tail -3 run_outline.sh >> ${array[0]}_run.sh
      chmod u+x ${array[0]}_run.sh

      if [ ${submit} -eq "1" ]; then
	  bsub -q 1nd ${array[0]}_run.sh
      fi
  fi
done < ${conf_file}

