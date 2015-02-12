#!/bin/bash

inputArray=("$@")

if [ ${#inputArray[@]} -lt 6 ]
then
    echo "not enough arguments! try again"
    exit -1
fi

scramdir=${inputArray[0]}
inputDir=${inputArray[1]}
outputDir=${inputArray[2]}
runMacro=${inputArray[3]}
soFile=${inputArray[4]}

workDir=`pwd`
echo `hostname`
echo "args:    $*"

cd ${scramdir}/src
eval `scramv1 runtime -sh`
cd $workDir

cp ${scramdir}/setRootEnv.C . 
cp ${scramdir}/rootlogon.C .
cp ${scramdir}/$runMacro  .
cp ${scramdir}/$soFile  .

for ((i=5; i<${#inputArray[@]}; i++)) 
do
    echo root -l -b -q ${runMacro}+\(\"${inputDir}${inputArray[i]}\",\"${outputDir}${inputArray[i]}\"\)
    root -l -b -q ${runMacro}+\(\"${inputDir}${inputArray[i]}\",\"${outputDir}${inputArray[i]}\"\)
done

status=`echo $?`
echo "Status - $status"

exit $status