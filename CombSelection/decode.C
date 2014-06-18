{

  Int_t test=903362060;

  Int_t nH=0, nW=0, nZ=0, nT=0, nB=0, nQ=0, nG=0, nP=0, nL=0;

  nH=test%10; test=(test-nH)/10;
  nW=test%10; test=(test-nW)/10;
  nZ=test%10; test=(test-nZ)/10;
  nT=test%10; test=(test-nT)/10;
  nB=test%10; test=(test-nB)/10;
  nQ=test%10; test=(test-nQ)/10;
  nG=test%10; test=(test-nG)/10;
  nP=test%10; test=(test-nP)/10;
  nL=test%10; test=(test-nL)/10;

  cout << nH << " Higgs bosons" << endl;
  cout << nW << " W bosons" << endl;
  cout << nZ << " Z bosons" << endl;
  cout << nT << " top quarks" << endl;
  cout << nB << " b quarks" << endl;
  cout << nQ << " light quarks" << endl;
  cout << nG << " gluons" << endl;
  cout << nP << " photons" << endl;
  cout << nL << " leptons" << endl;

}
