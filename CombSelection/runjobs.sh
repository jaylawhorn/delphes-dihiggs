#!/bin/bash
      
scramdir=$1
outputDir=$2
xsec=$3
events=$4
id=$5
runMacro=$6
soFile=$7
inarray=("$@")

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

loop=`expr ${#inarray[@]} - 7`

if [ $((loop%2)) -ne "0" ]; then
    echo "I/O broken, please try again"
    exit 1
fi

for ((i=0; i<$((loop/2)); i++)) 
do 
    x=`expr 7 + $i + $loop / 2`
    echo root -l -b -q ${runMacro}+\(\"${inarray[`expr $i + 7`]}\",${xsec},${events},${id},\"${outputDir}/${inarray[$x]}.root\"\)
    root -l -b -q ${runMacro}+\(\"${infile}\",${xsec},${events},${id},\"${outputDir}/${x}.root\"\)

done

status=`echo $?`
echo "Status - $status"

exit $status
