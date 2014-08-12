#!/bin/bash

#master submission script
conf_dir=/afs/cern.ch/work/j/jlawhorn/public/conf-files
#runMacro=ttbar_hack.C
runMacro=selection.C
outputBase=/afs/cern.ch/work/j/jlawhorn/public/ntuples

echo "Checking for new samples and"
echo "remaking config files."

while read line #loop over lines in ${conf_file}                                                                                                         
do
    array=($line)
    if [[ "${array[0]:0:1}" != "#" ]]; then
        if [ ${array[3]} -eq "1" ]; then
            testvar="$( /afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${array[5]} 2>&1 > /dev/null )"
            if [[ ${testvar} != "" ]]; then
                echo "Input folder does not exist... skipping" ${array[0]}
                continue
            fi
            filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${array[5]} | grep root`)
        elif [ ${array[3]} -eq "0" ]; then
            filelist=(`ls ${array[5]} | grep root`)
        fi
        size=${#filelist[@]}

	if [ ! -e "${conf_dir}/${array[0]}.CONF" ]; then
            echo "Creating filelist for" ${array[0]}
            if [ ${array[3]} -eq "1" ]; then
                echo "root://eoscms.cern.ch/"${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            elif [ ${array[3]} -eq "0" ]; then
                echo ${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            fi
            for ((i=0; i<${#filelist[@]}; i++))
            do
                echo ${filelist[i]} >> ${conf_dir}/${array[0]}.CONF
            done

	elif [[ -e "${conf_dir}/${array[0]}.CONF" ]] && [[ `wc -l ${conf_dir}/${array[0]}.CONF | cut -d' ' -f1` != $((size+1)) ]]; then
            smallsize=`wc -l ${conf_dir}/${array[0]}.CONF | cut -d' ' -f1`
            echo "Regenerating filelist for" ${array[0]}":"
            echo ${size} "files found in" ${array[5]} "but only" $((smallsize-1)) "listed in file"
            if [ ${array[3]} -eq "1" ]; then
                echo "root://eoscms.cern.ch/"${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            elif [ ${array[3]} -eq "0" ]; then
                echo ${array[5]} ${array[1]} ${array[4]} > ${conf_dir}/${array[0]}.CONF
            fi
            for ((i=0; i<${#filelist[@]}; i++))
            do
                echo ${filelist[i]} >> ${conf_dir}/${array[0]}.CONF
            done
        fi
    fi
done < "${conf_dir}/conf.txt"

echo " "
echo "Done checking config files."
echo "Submitting jobs."
echo " "

for file in ${conf_dir}/*CONF
do
    info=( $(head -n 1 $file) )
    outname=${file%.*}
    outname=${outname##*/}
    outputDir=${outputBase}/${outname}
    workDir=$CMSSW_BASE #`pwd`                                                                              
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
    
    mkdir -p ${outputDir}
    
    cp setRootEnv.C            $workDir
    cp rootlogon.C             $workDir
    cp $runMacro               $workDir
    cp $soFile                 $workDir
    
    sed 1d ${file} | while read line
    do
	#if [[ ! -e ${outputDir} ]]; then
	    #continue
	    #echo "Output directory doesn't exist. Not submitting."
	if [[ `bjobs -w 2> /dev/null | grep ${line}` ]]; then
	    continue
	    #echo "Job is currently running. Not submitting."
	elif [[ `grep "File broken" ${outputDir}/out.${line%.*}.txt` ]]; then
	    continue
	    #echo "Input file is broken. Not submitting."
	elif [[ -e ${outputDir}/${line} ]] && [[ `grep "Selection complete" ${outputDir}/out.${line%.*}.txt` ]]; then
	    continue
	    #echo "Output file exists and selection completed gracefully. Not submitting."
	#elif [[ -e ${outputDir}/${line} ]] && [[ `grep "Selection complete" \`egrep -lir "${line}" /afs/cern.ch/work/j/jlawhorn/public/ntuples/orphans/\`` ]]; then
	    #echo "Output is in orphaned folder... Not submitting."
	    #continue
	else
	    echo $script $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile 
	    bsub -o ${outputDir}/out.${line%.*}.txt -e ${outputDir}/err.${line%.*}.txt -q 8nh ${script} $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile 
	fi
    done 
done

echo " " 
echo "Done submitting jobs."