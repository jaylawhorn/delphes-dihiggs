#!/bin/bash

mkdir -p plots

#pt ================================================

root -l -q -b bbtt_upg_mt.C+\(\"ptTau1\",15,0,300,\"\#tau_{h}\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"ptTau2\",15,0,300,\"\#tau_{\#mu}\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"ptB1\",15,0,300,\"leading\ b\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"ptB2\",15,0,300,\"second\ b\ p_{T}\",\"Events\",10000\)

root -l -q -b bbtt_upg_et.C+\(\"ptTau1\",15,0,300,\"\#tau_{h}\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptTau2\",15,0,300,\"\#tau_{e}\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptB1\",15,0,300,\"leading\ b\ p_{T}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptB2\",15,0,300,\"second\ b\ p_{T}\",\"Events\",10000\)

root -l -q -b bbtt_upg_em.C+\(\"ptTau1\",15,0,300,\"\#tau_{\#mu}\ p_{T}\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"ptTau2\",15,0,300,\"\#tau_{e}\ p_{T}\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"ptB1\",15,0,300,\"leading\ b\ p_{T}\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"ptB2\",15,0,300,\"second\ b\ p_{T}\",\"Events\",100000\)

root -l -q -b bbtt_upg_tt.C+\(\"ptTau1\",15,0,300,\"leading\ \#tau_{h}\ p_{T}\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptTau2\",15,0,300,\"second\ \#tau_{h}\ p_{T}\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptB1\",15,0,300,\"leading\ b\ p_{T}\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptB2\",15,0,300,\"second\ b\ p_{T}\",\"Events\",1000\)

#eta ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"etaTau1\",16,-4,4,\"\#tau_{h}\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"etaTau2\",16,-4,4,\"\#tau_{\#mu}\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"etaB1\",16,-4,4,\"leading\ b\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_mt.C+\(\"etaB2\",16,-4,4,\"second\ b\ \#eta\",\"Events\",10000\)

root -l -q -b bbtt_upg_et.C+\(\"etaTau1\",16,-4,4,\"\#tau_{h}\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"etaTau2\",16,-4,4,\"\#tau_{e}\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"etaB1\",16,-4,4,\"leading\ b\ \#eta\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"etaB2\",16,-4,4,\"second\ b\ \#eta\",\"Events\",10000\)

root -l -q -b bbtt_upg_em.C+\(\"etaTau1\",16,-4,4,\"\#tau_{\#mu}\ \#eta\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"etaTau2\",16,-4,4,\"\#tau_{e}\ \#eta\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"etaB1\",16,-4,4,\"leading\ b\ \#eta\",\"Events\",100000\)
root -l -q -b bbtt_upg_em.C+\(\"etaB2\",16,-4,4,\"second\ b\ \#eta\",\"Events\",100000\)

root -l -q -b bbtt_upg_tt.C+\(\"etaTau1\",16,-4,4,\"leading\ \#tau_{h}\ \#eta\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"etaTau2\",16,-4,4,\"second\ \#tau_{h}\ \#eta\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"etaB1\",16,-4,4,\"leading\ b\ \#eta\",\"Events\",1000\)
root -l -q -b bbtt_upg_tt.C+\(\"etaB2\",16,-4,4,\"second\ b\ \#eta\",\"Events\",1000\)

#mTT ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"mTT\",25,0,500,\"m_{\#tau\#tau}\ \(vis\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"mTT\",25,0,500,\"m_{\#tau\#tau}\ \(vis\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"mTT\",25,0,500,\"m_{\#tau\#tau}\ \(vis\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"mTT\",10,0,500,\"m_{\#tau\#tau}\ \(vis\)\",\"Events\",1000\)

#mBB ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"mBB1\",25,0,500,\"m_{bb}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"mBB1\",25,0,500,\"m_{bb}\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"mBB1\",25,0,500,\"m_{bb}\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"mBB1\",10,0,500,\"m_{bb}\",\"Events\",1000\)

#mHH ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"mHH\",30,100,700,\"m_{HH}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"mHH\",30,100,500,\"m_{HH}\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"mHH\",30,100,700,\"m_{HH}\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"mHH\",10,100,700,\"m_{HH}\",\"Events\",1000\)

#mT2 ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"mt2pileup\",25,0,500,\"m_{T2}\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"mt2pileup\",25,0,500,\"m_{T2}\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"mt2pileup\",25,0,500,\"m_{T2}\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"mt2pileup\",10,0,500,\"m_{T2}\",\"Events\",1000\)

#ptTT ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"ptTT\",25,0,500,\"p_{T}\ \(\#tau\#tau\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptTT\",25,0,500,\"p_{T}\ \(\#tau\#tau\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"ptTT\",25,0,500,\"p_{T}\ \(\#tau\#tau\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptTT\",25,0,500,\"p_{T}\ \(\#tau\#tau\)\",\"Events\",1000\)

#ptBB ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"ptBB1\",25,0,500,\"p_{T}\ \(bb\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptBB1\",25,0,500,\"p_{T}\ \(bb\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"ptBB1\",25,0,500,\"p_{T}\ \(bb\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptBB1\",10,0,500,\"p_{T}\ \(bb\)\",\"Events\",1000\)

#ptHH ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"ptHH\",25,0,500,\"p_{T}\ \(HH\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"ptHH\",25,0,500,\"p_{T}\ \(HH\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"ptHH\",25,0,500,\"p_{T}\ \(HH\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"ptHH\",10,0,500,\"p_{T}\ \(HH\)\",\"Events\",1000\)

#drBB ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"sqrt\(\(etaB1-etaB2\)**2+\(phiB1-phiB2\)**2\)\",16,0,8,\"\#DeltaR\(bb\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"sqrt\(\(etaB1-etaB2\)**2+\(phiB1-phiB2\)**2\)\",16,0,8,\"\#DeltaR\(bb\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"sqrt\(\(etaB1-etaB2\)**2+\(phiB1-phiB2\)**2\)\",16,0,8,\"\#DeltaR\(bb\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"sqrt\(\(etaB1-etaB2\)**2+\(phiB1-phiB2\)**2\)\",16,0,8,\"\#DeltaR\(bb\)\",\"Events\",1000\)

#drBB ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"sqrt\(\(etaTau1-etaTau2\)**2+\(phiTau1-phiTau2\)**2\)\",16,0,8,\"\#DeltaR\(\#tau\#tau\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"sqrt\(\(etaTau1-etaTau2\)**2+\(phiTau1-phiTau2\)**2\)\",16,0,8,\"\#DeltaR\(\#tau\#tau\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"sqrt\(\(etaTau1-etaTau2\)**2+\(phiTau1-phiTau2\)**2\)\",16,0,8,\"\#DeltaR\(\#tau\#tau\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"sqrt\(\(etaTau1-etaTau2\)**2+\(phiTau1-phiTau2\)**2\)\",16,0,8,\"\#DeltaR\(\#tau\#tau\)\",\"Events\",1000\)

#drBB ==============================================

root -l -q -b bbtt_upg_mt.C+\(\"sqrt\(\(etaBB1-etaTT\)**2+\(phiBB1-phiTT\)**2\)\",16,0,8,\"\#DeltaR\(HH\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"sqrt\(\(etaBB1-etaTT\)**2+\(phiBB1-phiTT\)**2\)\",16,0,8,\"\#DeltaR\(HH\)\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"sqrt\(\(etaBB1-etaTT\)**2+\(phiBB1-phiTT\)**2\)\",16,0,8,\"\#DeltaR\(HH\)\",\"Events\",100000\)
root -l -q -b bbtt_upg_tt.C+\(\"sqrt\(\(etaBB1-etaTT\)**2+\(phiBB1-phiTT\)**2\)\",16,0,8,\"\#DeltaR\(HH\)\",\"Events\",1000\)

#bdt ===============================================

root -l -q -b bbtt_upg_mt.C+\(\"bdtVal\",40,-1,1,\"BDT\ Value\",\"Events\",10000\)
root -l -q -b bbtt_upg_et.C+\(\"bdtVal\",40,-1,1,\"BDT\ Value\",\"Events\",10000\)
root -l -q -b bbtt_upg_em.C+\(\"bdtVal\",40,-1,1,\"BDT\ Value\",\"Events\",100000\)

mv *png plots/

#root -l -q bbtt_upg_tt.C+\(\"mt2pileup\",10,0,500,\"m_{T2}\",\"Events\",1000\)
#root -l -q bbtt_upg_mt.C+\(\"bdtVal\",20,-0.5,0.5,\"BDT\ Value\",\"Events\",1000\)
#root -l -q bbtt_upg_et.C+\(\"bdtVal\",20,-0.5,0.5,\"BDT\ Value\",\"Events\",1000\)
#root -l -q bbtt_upg_em.C+\(\"bdtVal\",20,-0.5,0.5,\"BDT\ Value\",\"Events\",10000\)