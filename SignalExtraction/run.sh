#!/bin/bash

#pt ================================================

#root -l -q -b draw_obj.C+\(\"ptTau1\",\"\(ptTau2\>45\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",0,0,1,25,0,500,\"leading\ \#tau\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptTau2\",\"\(ptTau1\>45\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",0,0,1,15,0,300,\"second\ \#tau\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB1\",\"\(ptB2\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",0,1,0,25,0,500,\"leading\ b\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB2\",\"\(ptB1\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",0,1,0,15,0,300,\"second\ b\ p_{T}\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"ptTau1\",\"\(ptTau2\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",1,0,1,25,0,500,\"\#tau_{h}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptTau2\",\"\(ptTau1\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",1,0,1,25,0,500,\"\#tau_{\#mu}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB1\",\"\(ptB2\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",1,1,0,25,0,500,\"leading\ b\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB2\",\"\(ptB1\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",1,1,0,15,0,300,\"second\ b\ p_{T}\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"ptTau1\",\"\(ptTau2\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",2,0,1,25,0,500,\"\#tau_{h}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptTau2\",\"\(ptTau1\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",2,0,1,25,0,500,\"\#tau_{e}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB1\",\"\(ptB2\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",2,1,0,25,0,500,\"leading\ b\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB2\",\"\(ptB1\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",2,1,0,15,0,300,\"second\ b\ p_{T}\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"ptTau1\",\"\(ptTau2\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",3,0,1,25,0,500,\"\#tau_{\#mu}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptTau2\",\"\(ptTau1\>30\&\&abs\(etaTau1\)\<4.0\&\&abs\(etaTau2\)\<4.0\)\",3,0,1,25,0,500,\"\#tau_{e}\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB1\",\"\(ptB2\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",3,1,0,25,0,500,\"leading\ b\ p_{T}\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"ptB2\",\"\(ptB1\>30\&\&abs\(etaB1\)\<4.0\&\&abs\(etaB2\)\<4.0\)\",3,1,0,15,0,300,\"second\ b\ p_{T}\",\"Events\"\)

#eta ===============================================

#root -l -q -b draw_obj.C+\(\"etaTau1\",\"\(ptTau1\>45\&\&ptTau2\>45\&\&abs\(etaTau2\)\<4.0\)\",0,0,1,10,-5,5,\"leading\ \#tau\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaTau2\",\"\(ptTau1\>45\&\&ptTau2\>45\&\&abs\(etaTau1\)\<4.0\)\",0,0,1,10,-5,5,\"second\ \#tau\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB1\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB2\)\<4.0\)\",0,1,0,10,-5,5,\"leading\ b\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB2\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB1\)\<4.0\)\",0,1,0,10,-5,5,\"second\ b\ \#eta\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"etaTau1\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau2\)\<4.0\)\",1,0,1,10,-5,5,\"\#tau_{h}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaTau2\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau1\)\<4.0\)\",1,0,1,10,-5,5,\"\#tau_{\#mu}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB1\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB2\)\<4.0\)\",1,1,0,10,-5,5,\"leading\ b\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB2\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB1\)\<4.0\)\",1,1,0,10,-5,5,\"second\ b\ \#eta\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"etaTau1\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau2\)\<4.0\)\",2,0,1,10,-5,5,\"\#tau_{h}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaTau2\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau1\)\<4.0\)\",2,0,1,10,-5,5,\"\#tau_{e}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB1\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB2\)\<4.0\)\",2,1,0,10,-5,5,\"leading\ b\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB2\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB1\)\<4.0\)\",2,1,0,10,-5,5,\"second\ b\ \#eta\",\"Events\"\)

#root -l -q -b draw_obj.C+\(\"etaTau1\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau1\)\<4.0\)\",3,0,1,10,-5,5,\"\#tau_{\#mu}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaTau2\",\"\(ptTau1\>30\&\&ptTau2\>30\&\&abs\(etaTau2\)\<4.0\)\",3,0,1,10,-5,5,\"\#tau_{e}\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB1\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB1\)\<4.0\)\",3,1,0,10,-5,5,\"leading\ b\ \#eta\",\"Events\"\)
#root -l -q -b draw_obj.C+\(\"etaB2\",\"\(ptB1\>30\&\&ptB2\>30\&\&abs\(etaB2\)\<4.0\)\",3,1,0,10,-5,5,\"second\ b\ \#eta\",\"Events\"\)

#m(tt) =============================================
#root -l -q -b draw.C+\(\"mTT\",\"\(1\)\",0,8,0,200,\"m_{vis}\(\#tau\#tau\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mTT\",\"\(1\)\",1,8,0,200,\"m_{vis}\(\#tau\#tau\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mTT\",\"\(1\)\",2,8,0,200,\"m_{vis}\(\#tau\#tau\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mTT\",\"\(1\)\",3,8,0,200,\"m_{vis}\(\#tau\#tau\)\",\"Events\"\)

#m(bb) =============================================
#root -l -q -b draw.C+\(\"mBB1\",\"\(1\)\",0,10,0,250,\"m\(bb\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mBB1\",\"\(1\)\",1,10,0,250,\"m\(bb\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mBB1\",\"\(1\)\",2,10,0,250,\"m\(bb\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mBB1\",\"\(1\)\",3,10,0,250,\"m\(bb\)\",\"Events\"\)

#mT2 ===============================================
#root -l -q -b draw.C+\(\"mt2\",\"\(1\)\",0,20,0,500,\"m_{T2}\",\"Events\"\)
#root -l -q -b draw.C+\(\"mt2\",\"\(1\)\",1,20,0,500,\"m_{T2}\",\"Events\"\)
#root -l -q -b draw.C+\(\"mt2\",\"\(1\)\",2,20,0,500,\"m_{T2}\",\"Events\"\)
#root -l -q -b draw.C+\(\"mt2\",\"\(1\)\",3,20,0,500,\"m_{T2}\",\"Events\"\)

#m(HH) =============================================
#root -l -q -b draw.C+\(\"mHH\",\"\(1\)\",0,20,0,1500,\"m\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mHH\",\"\(1\)\",1,60,0,1500,\"m\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mHH\",\"\(1\)\",2,60,0,1500,\"m\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"mHH\",\"\(1\)\",3,60,0,1500,\"m\(HH\)\",\"Events\"\)

#pT(H) ============================================
root -l draw.C+\(\"ptHH\",\"\(mBB1\>100\&\&mBB1\<150\&\&mTT\>50\&\&mTT\<150\&\&mHH\>300\&\&mt2\>130\)\",0,20,0,500,\"p_{T}\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"ptHH\",\"\(1\)\",1,20,0,500,\"p_{T}\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"ptHH\",\"\(1\)\",2,20,0,500,\"p_{T}\(HH\)\",\"Events\"\)
#root -l -q -b draw.C+\(\"ptHH\",\"\(1\)\",3,20,0,500,\"p_{T}\(HH\)\",\"Events\"\)