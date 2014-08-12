#!/bin/bash

input_loc=/store/group/phys_higgs/future/sixie/Madgraph/BBJJ_M60To200_14TeV_v2/LHE/
output_loc=/store/user/jlawhorn/qcd_bbjj_meden_agan/

for file in `/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls ${input_loc}`
do
    echo    bsub -q 2nd -o /afs/cern.ch/work/j/jlawhorn/logfile_qcd_${file%.*}.txt -e /afs/cern.ch/work/j/jlawhorn/errfile_qcd_${file%.*}.txt -J ${file} run.sh ${file} $input_loc $output_loc JetStudies_Phase_II_140PileUp_conf4.tcl
    #bsub -q 2nd -o /afs/cern.ch/work/j/jlawhorn/logfile_qcd_${line%.*}.txt -e /afs/cern.ch/work/j/jlawhorn/errfile_qcd_${line%.*}.txt -J ${file} run.sh ${line} $input_loc $output_loc JetStudies_Phase_II_140PileUp_conf4.tcl
    ./run.sh ${file} $input_loc $output_loc QcdStudies_Phase_II_140PileUp_conf4.tcl
    exit
done 