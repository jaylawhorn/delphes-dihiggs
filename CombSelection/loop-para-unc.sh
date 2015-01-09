#!/bin/bash

#master submission script
dir=/afs/cern.ch/work/j/jlawhorn/para-eff
runMacro=flag-events.C

rm ${dir}/*-small.root

for file in ${dir}/*root
do
    outname=${file##*/}
    outname=${outname%.*}
   
    echo root -l -q flag-events.C+\(\"${file}\",\"${dir}/${outname}-small.root\"\)
    root -l -q flag-events.C+\(\"${file}\",\"${dir}/${outname}-small.root\"\)
done
