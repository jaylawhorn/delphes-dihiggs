#!/bin/bash

#master submission script
conf_dir=.
runMacro=select.C
#outputBase=/afs/cern.ch/work/j/jlawhorn/public/hh-bacon-feb-11/
outputBase=/afs/cern.ch/work/j/jlawhorn/public/tt-bacon-feb-11/
echo "Checking for new samples"

#filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/group/upgrade/di_higgs_backgrounds/gFHHTobbtt_TuneZ2_14TeV_madgraph/Bacon/ | grep root`)
filelist=(`/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls /store/group/upgrade/di_higgs_backgrounds/TT_jets/Bacon/ | grep root`)

#echo "root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/gFHHTobbtt_TuneZ2_14TeV_madgraph/Bacon/" > ${conf_dir}/hh.CONF
echo "root://eoscms.cern.ch//store/group/upgrade/di_higgs_backgrounds/TT_jets/Bacon/" > ${conf_dir}/tt.CONF

for ((i=0; i<${#filelist[@]/50+1}; i+=50))
do
    echo ${filelist[@]:i:50} >> ${conf_dir}/tt.CONF
done

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

mkdir -p ${outputBase}


cp setRootEnv.C            $workDir
cp rootlogon.C             $workDir
cp $runMacro               $workDir
cp $soFile                 $workDir

info=( $(head -n 1 ${conf_dir}/tt.CONF) )
outname=${file%.*}
outname=${outname##*/}

i=0

sed 1d ${conf_dir}/tt.CONF | while read line
do
    echo $script $workDir ${info[0]} $outputBase $runMacro $soFile $line
	#./$script $workDir ${info[0]} $outputBase $runMacro $soFile $line
    bsub -o ${outputBase}out.${i}.txt -e ${outputBase}err.${i}.txt -q 8nh ${script} $workDir ${info[0]} $outputBase $runMacro $soFile $line
    i=$((i+1))
done