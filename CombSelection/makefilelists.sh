#!/bin/bash

#master run script for pre-processing

if [ ${#} -lt "1" ]; then
    echo "Please use the following syntax: " 
    echo "    ./master.sh <conf_file>"
    echo " "
    echo "<conf_file> is a text file containing the names of all samples"
    exit
fi

conf_file=$1

while read line #loop over lines in ${conf_file}
do
    array=($line)
    if [ "${array[0]:0:1}" != "#" ]; then
	echo "**********" ${array[0]} "**********"
	if [ ${array[3]} -eq "1" ]; then
	    filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${array[5]} | grep root`)
	elif [ ${array[3]} -eq "0" ]; then
	    filelist=(`ls ${array[5]} | grep root`)
	fi
	size=${#filelist[@]}
	
	if [ ! -e "${array[0]}.CONF" ]; then
	    echo "Creating filelist for" ${array[0]}
	    if [ ${array[3]} -eq "1" ]; then
		echo "root://eoscms.cern.ch/"${array[5]} ${array[1]} ${array[4]} > ${array[0]}.CONF
	    elif [ ${array[3]} -eq "0" ]; then
		echo ${array[5]} ${array[1]} ${array[4]} > ${array[0]}.CONF
	    fi
	    for ((i=0; i<${#filelist[@]}; i++))
	    do 
		echo ${filelist[i]} >> ${array[0]}.CONF
	    done
	    
	elif [[ -e "${array[0]}.CONF" ]] && [[ `wc -l ${array[0]}.CONF | cut -d' ' -f1` != $((size+1)) ]]; then
	    smallsize=`wc -l ${array[0]}.CONF | cut -d' ' -f1` 
	    echo "Regenerating filelist for" ${array[0]}":"
	    echo ${size} "files found in" ${array[5]} "but only" $((smallsize-1)) "listed in file"
	    if [ ${array[3]} -eq "1" ]; then
		echo "root://eoscms.cern.ch/"${array[5]} ${array[1]} > ${array[0]}.CONF
	    elif [ ${array[3]} -eq "0" ]; then
		echo ${array[5]} ${array[1]} > ${array[0]}.CONF
	    fi
	    for ((i=0; i<${#filelist[@]}; i++))
	    do 
		echo ${filelist[i]} >> ${array[0]}.CONF
	    done
	else
	    echo "Filelist exists and no new files have appeared since generation"
	fi
    fi
done < ${conf_file}