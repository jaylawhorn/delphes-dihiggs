cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/vbf* .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/Selection/rootlogon.C .

workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd - 
status=`echo $?`
echo "${h}: Status - $status"
exit $status