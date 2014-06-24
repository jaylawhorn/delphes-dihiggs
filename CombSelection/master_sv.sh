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
      if [ ${array[3]} -eq "1" ]; then
	  filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${array[5]} | grep ${array[0]}`)
      elif [ ${array[3]} -eq "0" ]; then
	  filelist=(`ls ${array[5]} | grep ${array[0]}`)
      fi
      outarray=()
      for ((i=0; i<${#filelist[@]}; i++)) 
      do
	  if [[ ${filelist[$i]} == *root ]]; then
	      outarray+=(${filelist[$i]})
	  fi
      done
      echo ${#outarray[@]} "files found in" ${array[5]}

      for ((i=0; i<${#outarray[@]}; i++))
      do
	  if [ ${array[3]} -eq "1" ]; then
	      echo "root://eoscms.cern.ch/"${array[5]}"/"${outarray[i]} >> ${array[0]}.txt
	  elif [ ${array[3]} -eq "0" ]; then
	      echo ${array[5]}"/"${outarray[i]} >> ${array[0]}.txt
	  fi
      done

      echo `wc -l ${array[0]}.txt | cut -d' ' -f1` "files listed in" ${array[0]}.txt

      root -l -q -b getall.C+\(\"${array[0]}.txt\"\) | tail -1 > ${array[0]}_events.txt

      nevents=`cat ${array[0]}_events.txt`
      echo $nevents "events found"

      outputDir=/afs/cern.ch/work/a/arapyan/public/comb_ntuples/${array[0]}
      workDir=$CMSSW_BASE #`pwd`  
      runMacro=selection.C
      soFile=`echo $runMacro | sed 's/\./_/'`.so
      script=runjobs.sh
      
      
      # check a few things
      if [ ! "$CMSSW_BASE" ]; then
	  echo "-------> error: define cms environment."
	  exit 1
      fi
      if [ macros/$runMacro -nt macros/$soFile ]; then
	  echo "-------> error: forgot to recompile run macro."
	  exit 1
      fi

      cp setRootEnv.C            $workDir
      cp rootlogon.C             $workDir
      cp $runMacro               $workDir
      cp $soFile                 $workDir

      n=0
      
      while read line2
	do
	if [ ${submit} -eq "1" ]; then
	    echo  $script $workDir $outputDir ${line2} ${array[1]} $nevents ${array[4]} $n $runMacro $soFile
	    bsub  -o out.%J  -q 2nd $script  $workDir $outputDir ${line2}  ${array[1]} $nevents ${array[4]} $n $runMacro $soFile
	fi
	n=$((n+1))
      done < ${array[0]}.txt 
  fi
done < ${conf_file}



  
