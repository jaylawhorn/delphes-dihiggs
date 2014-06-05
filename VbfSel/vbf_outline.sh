cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/VbfSel/selection* .
cp /afs/cern.ch/user/j/jlawhorn/delphes-dihiggs/VbfSel/rootlogon.C .

workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`
cd - 
status=`echo $?`
echo "${h}: Status - $status"
exit $status