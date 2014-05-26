#!/bin/bash
#------------------------------------------------------------
# Submit a batch of jobs to lxbtch
#
# example command:
#             root_script     conf_file output_location
# ./submit.sh submitDelphes.C xsec.txt  /afs/cern.ch/work/k/klawhorn/SnowmassSamples/
# 
# conf_file has format (no leading "#")
# sample_type cross_section
#
# Note: conf file needs empty last line... clunky I know.
# Jay Lawhorn 11/4/13
#------------------------------------------------------------

root_script=$1
 output_loc=$2

#xsec="225.2"
 xsec="2.92"

for file in `ls /afs/cern.ch/work/k/klawhorn/hh_param/*root`
do
    echo "   bsub -q 8nh -W 120 -J $file run.sh ${root_script} ${file} ${xsec} ${output_loc}" 
    bsub -q 8nh -W 120 -J $file run_mine.sh ${root_script} ${file} ${xsec} ${output_loc} # submit to lxbtch
    #./run_mine.sh ${root_script} ${file} ${xsec} ${output_loc}
done
