cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/selection* .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/rootlogon.C .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/mt2.hh .
workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd - 
status=`echo $?`
echo "${h}: Status - $status"
exit $status