#!/bin/bash

#for sample in BB-4p-0-300-v1510_14TEV BB-4p-300-700-v1510_14TEV BB-4p-700-1300-v1510_14TEV BB-4p-1300-2100-v1510_14TEV BB-4p-2100-100000_14TEV BBB-4p-0-600-v1510_14TEV BBB-4p-600-1300-v1510_14TEV BBB-4p-1300-100000-v1510_14TEV LL-4p-0-100-v1510_14TEV LL-4p-100-200-v1510_14TEV LL-4p-200-500-v1510_14TEV LL-4p-500-900-v1510_14TEV LL-4p-900-1400-v1510_14TEV LL-4p-1400-100000-v1510_14TEV tt-4p-0-600-v1510_14TEV tt-4p-600-1100-v1510_14TEV tt-4p-1100-1700-v1510_14TEV tt-4p-1700-2500-v1510_14TEV tt-4p-2500-100000-v1510_14TEV
for sample in HHToTTBB_14TeV
#for sample in tt-4p-0-600-v1510_14TEV tt-4p-600-1100-v1510_14TEV tt-4p-1100-1700-v1510_14TEV tt-4p-1700-2500-v1510_14TEV tt-4p-2500-100000-v1510_14TEV
do
  ./mergeSelFiles.sh ${sample} /afs/cern.ch/work/k/klawhorn/resTestsOff/PhaseII/Configuration4v2 /afs/cern.ch/work/k/klawhorn
  #./mergeSelFiles.sh ${sample} /afs/cern.ch/work/k/klawhorn/jetPerformance/PhaseII/Configuration4v2 /afs/cern.ch/work/k/klawhorn
done