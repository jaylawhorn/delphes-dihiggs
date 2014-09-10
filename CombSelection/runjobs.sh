#!/bin/bash
      
scramdir=$1
outputDir=$2
inputDir=$3
inputFile=$4
xsec=$5
id=$6
runMacro=$7
soFile=$8


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

echo root -l -b -q ${runMacro}+\(\"${inputDir}${inputFile}\",${xsec},${id},\"${inputFile}\"\)
root -l -b -q ${runMacro}+\(\"${inputDir}${inputFile}\",${xsec},${id},\"${inputFile}\"\)

cp ${inputFile} $outputDir/${inputFile}

status=`echo $?`
echo "Status - $status"

exit $status
