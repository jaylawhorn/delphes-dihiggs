void make_training() {

  TFile *f1 = new TFile("/afs/cern.ch/work/j/jlawhorn/public/ntuples/HHToTTBB_14TeV.root","READ");
  TTree *tree1 = (TTree*) f1->Get("Events");
  TFile *f2 = new TFile("/afs/cern.ch/work/j/jlawhorn/public/ntuples/tt.root","READ");
  TTree *tree2 = (TTree*) f2->Get("Events");

  //TFile *newfile = new TFile("/afs/cern.ch/work/j/jlawhorn/public/mva-training/mt.root", "RECREATE");
  TFile *newfile = new TFile("mt_test.root", "RECREATE");

  TTree *treeS = (TTree*) tree1->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==3 && ptTau1>30 && ptTau2>20 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeS->SetName("treeS");
  treeS->SetTitle("treeS");
  TTree *treeB = (TTree*) tree2->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==3 && ptTau1>30 && ptTau2>20 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeB->SetName("treeB");
  treeB->SetTitle("treeB");

  newfile->Write();
  newfile->Close();
  /*
  newfile=0; treeS=0; treeB=0;

  newfile = new TFile("/afs/cern.ch/work/j/jlawhorn/public/mva-training/et.root", "RECREATE");

  treeS = (TTree*) tree1->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==2 && ptTau1>30 && ptTau2>20 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeS->SetName("treeS");
  treeS->SetTitle("treeS");
  treeB = (TTree*) tree2->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==2 && ptTau1>30 && ptTau2>20 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeB->SetName("treeB");
  treeB->SetTitle("treeB");

  newfile->Write();
  newfile->Close();

  newfile=0; treeS=0; treeB=0;

  newfile = new TFile("/afs/cern.ch/work/j/jlawhorn/public/mva-training/em.root", "RECREATE");

  treeS = (TTree*) tree1->CopyTree("isBBTT==1 && tauCat1==3 && tauCat2==2 && ptTau1>20 && ptTau2>20 && tauIso1<0.4 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeS->SetName("treeS");
  treeS->SetTitle("treeS");
  treeB = (TTree*) tree2->CopyTree("isBBTT==1 && tauCat1==3 && tauCat2==2 && ptTau1>20 && ptTau2>20 && tauIso1<0.4 && tauIso2<0.4 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  treeB->SetName("treeB");
  treeB->SetTitle("treeB");

  newfile->Write();
  newfile->Close();

  newfile=0; treeS=0; treeB=0;

  newfile = new TFile("/afs/cern.ch/work/j/jlawhorn/public/mva-training/tt.root", "RECREATE");

  treeS = (TTree*) tree1->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && mTT>60 && mTT<140 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  //treeS = (TTree*) tree1->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && bTag1>0 && bTag2>0 && ptB1>20 && ptB2>20");
  treeS->SetName("treeS");
  treeS->SetTitle("treeS");
  treeB = (TTree*) tree2->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && (bTag1==2||bTag1==3||bTag1==6||bTag1==7) && (bTag2==1||bTag2==3||bTag2==6||bTag2==7) && ptB1>20 && ptB2>20 && mTT>60 && mTT<140 && sqrt( (etaTau1-etaTau2)**2 + (phiTau1-phiTau2)**2 )>0.4");
  //treeB = (TTree*) tree2->CopyTree("isBBTT==1 && tauCat1==1 && tauCat2==1 && ptTau1>45 && ptTau2>45 && bTag1>0 && bTag2>0 && ptB1>20 && ptB2>20");
  treeB->SetName("treeB");
  treeB->SetTitle("treeB");

  newfile->Write();
  newfile->Close();
  */
}
