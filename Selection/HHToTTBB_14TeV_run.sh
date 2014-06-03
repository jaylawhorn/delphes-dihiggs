cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/selection* .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/rootlogon.C .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/mt2.hh .
workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd - 
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_0.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/0.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_1.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/1.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_10.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/2.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_11.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/3.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_12.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/4.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_13.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/5.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_14.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/6.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_15.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/7.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_16.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/8.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_17.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/9.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_18.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/10.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_19.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/11.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_2.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/12.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_3.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/13.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_4.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/14.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_5.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/15.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_6.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/16.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_7.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/17.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_8.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/18.root\"\)
root -l -q -b selection.C+\(\"/afs/cern.ch/work/j/jlawhorn/public/HHToTTBB_14TeV//HHToTTBB_14TeV_9.root\",2.92,97968,\"/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV/19.root\"\)
status=`echo $?`
echo "${h}: Status - $status"
exit $status