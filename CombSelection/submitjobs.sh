#!/bin/bash

#master submission script

directory=${1}
charstr=${2}

for file in ${directory}/${charstr}
do
    info=( $(head -n 1 $file) )
    outname=${file%.*}
    outname=${outname##*/}
    outputDir=/afs/cern.ch/work/j/jlawhorn/public/ntuples/${outname}
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
    
    mkdir -p ${outputDir}
    
    cp setRootEnv.C            $workDir
    cp rootlogon.C             $workDir
    cp $runMacro               $workDir
    cp $soFile                 $workDir
    
    sed 1d ${file} | while read line
    do
	if [[ ! -e ${outputDir} ]]; then
	    continue
	    #echo "Output directory doesn't exist. Not submitting."
	elif [[ `bjobs -w | grep ${line}` ]]; then
	    continue
	    #echo "Job is currently running. Not submitting."
	elif [[ `grep "File broken" ${outputDir}/out.${line%.*}.txt` ]]; then
	    continue
	    #echo "Input file is broken. Not submitting."
	elif [[ -e ${outputDir}/${line} ]] && [[ `grep "Selection complete" ${outputDir}/out.${line%.*}.txt` ]]; then
	    continue
	    #echo "Output file exists and selection completed gracefully. Not submitting."
	elif [[ -e ${outputDir}/${line} ]] && [[ `grep "Selection complete" \`egrep -lir "${line}" /afs/cern.ch/work/j/jlawhorn/public/ntuples/orphans/\`` ]]; then
	    echo "Output is in orphaned folder... Not submitting."
	    continue
	else
	    echo $script $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile 
	    #./${script} $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile 
	    bsub -o ${outputDir}/out.${line%.*}.txt -e ${outputDir}/err.${line%.*}.txt -q 8nh ${script} $workDir $outputDir ${info[0]} ${line} ${info[1]} ${info[2]} $runMacro $soFile 
	fi
    done 
done
