#!/bin/bash

IN_FILE=$1
IN_DIR=$2
OUT_DIR=$3
DET_CARD=$4

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

# some basic printing                                                                                                                                                      
echo " "; echo "${h}: Show who and where we are";
echo " "
echo " user executing: "`id`;
echo " running on    : "`hostname`;
echo " executing in  : "`pwd`;
echo " submitted from: $HOSTNAME";
echo " ";

#initialize the CMSSW environment
echo " "; echo "${h}: Initialize CMSSW (in $CMSSW_BASE)"; echo " "
WORK_DIR=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd -

#setup delphes
cp /afs/cern.ch/work/j/jlawhorn/public/for-prod/meden-agan.tar .
tar -xvf $WORK_DIR/meden-agan.tar
cd Delphes
chmod -R 777 *
./configure
make

#copy input file
xrdcp root://eoscms.cern.ch/${IN_DIR}/${IN_FILE} .
gunzip ${IN_FILE}

cp scripts/qcd14.py .
cmsRun qcd14.py ${IN_FILE%.*.*} $RANDOM

#produce minbias file
sed -i 's/50000/50/g' Cards/MinBias_TuneZ2star_14TeV_pythia6_cff.py
cmsRun Cards/MinBias_TuneZ2star_14TeV_pythia6_cff.py
mv GenEvent_ASCII.dat GenEvent_ASCII.hepmc

./hepmc2pileup MinBias.pileup GenEvent_ASCII.hepmc
if [ $? -ne 0 ]; then
    exit 15
fi

rm -rf GenEvent_ASCII.hepmc

echo Cards/$DET_CARD

#sed -i 's/#Set MaxEvents/Set MaxEvents/' Cards/$DET_CARD

OUT_FILE=${IN_FILE%.*.*}_ma.root

./DelphesCMSFWLite Cards/$DET_CARD ${OUT_FILE} ${IN_FILE%.*.*}.root
if [ $? -ne 0 ]; then
    exit 16
fi

cmsStage $OUT_FILE $OUT_DIR/$OUT_FILE
