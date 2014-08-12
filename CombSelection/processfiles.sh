#!/bin/bash

indir=$1

for folder in $( ls -d ${indir}*/ )
do
    files=(${folder}/*.root)
    if [[ ! ${#files[@]} -gt 1 ]]; then
	echo "no files in folder... skipping" ${folder}
	continue
    elif [[ -e "${folder%/*}.root" ]]; then
	old="${folder%/*}.root"
	newest=$( ls -rt ${folder}/*.root | tail -n 1 )
	if [[ ${newest} -nt ${old} ]]; then
	    echo "new files have appeared... recombining" ${folder}
	    rm ${old}
	else
	    echo "no new files... leaving" ${folder}
	    continue
	fi
    else
	echo "combining" ${folder}
    fi
    ls ${folder}/*.root > ${folder}temp.txt
    echo root -l -q combinefiles.C+\(\"${folder}temp.txt\",\"${folder%/*}.root\"\)
   root -l -q combinefiles.C+\(\"${folder}temp.txt\",\"${folder%/*}.root\"\)
    #rm ${folder}temp.txt
done