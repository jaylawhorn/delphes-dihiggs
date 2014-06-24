#!/bin/bash
      
scramdir=$1
outputDir=$2
infile=$3
xsec=$4
events=$5
id=$6
x=$7
runMacro=$8
soFile=$9

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

mkdir ${outputDir}
root -l -b -q ${runMacro}+\(\"${infile}\",${xsec},${events},${id},\"${outputDir}/${x}.root\"\)

status=`echo $?`
echo "Status - $status"

exit $status
