{
  bool foundIt=true;
  cout << " ROOT Version: " << gROOT->GetVersion() << "." << endl;
  char *cmsbase=gSystem->Getenv("CMSSW_BASE");
  if (cmsbase==NULL) {
    cout << " CMSSW environment has not been set up." << endl;
    foundIt=false;
  } else {
    TString rfitpath("/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
    //TString rfitpath("/afs/cern.ch/cms/slc6_amd64_gcc472/lcg/roofit/5.32.03/include/");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
    
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
    str.ReplaceAll(" -c ", " -std=c++0x -c ");
    gSystem->SetMakeSharedLib(str);
    
    str = gSystem->GetMakeExe();
    str.ReplaceAll(" -c ", " -std=c++0x -c ");
    gSystem->SetMakeExe(str);

    cout << " CMSSW environment is set up." << endl;
    TString addedLibs(gSystem->GetLibraries());
    if(!addedLibs.Contains("setRootEnv_C.so")) {
      gROOT->Macro("setRootEnv.C+");
    }
    gSystem->Load("$CMSSW_BASE/src/Delphes/libDelphes.so");
    loadLibraries("libTauAnalysisSVFitHelper.so");
    loadLibraries("libTauAnalysisCandidateTools.so");
    loadLibraries("libJECJECHelper.so");
    loadLibraries("libPhysicsToolsKinFitter.so");
    gROOT->ProcessLine(".include $CMSSW_BASE/src/Delphes");
    gROOT->ProcessLine(".include $CMSSW_BASE/src/Delphes/external");
    gROOT->Macro("../Utils/hhMVA.c+");
    //gROOT->Macro("DelWeight.cc+");
    //gROOT->Macro("../Utils/HttStyles.cc+");
    //gROOT->Macro("../Utils/CPlot.cc+");
  }
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
