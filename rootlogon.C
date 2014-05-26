{
  bool foundIt=true;
  cout << " ROOT Version: " << gROOT->GetVersion() << "." << endl;
  char *cmsbase=gSystem->Getenv("CMSSW_BASE");
  if (cmsbase==NULL) {
    cout << " CMSSW environment has not been set up." << endl;
    foundIt=false;
  } else {
    TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());

    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
    cout << " CMSSW environment is set up." << endl;
    gSystem->Load("/afs/cern.ch/user/j/jlawhorn/Delphes-3.0.10/libDelphes.so");
    gROOT->ProcessLine(".include /afs/cern.ch/user/j/jlawhorn/Delphes-3.0.10");
    gROOT->ProcessLine(".include /afs/cern.ch/user/j/jlawhorn/Delphes-3.0.10/external");
    gROOT->Macro("../Utils/MitStyleRemix.cc++");
  }
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
