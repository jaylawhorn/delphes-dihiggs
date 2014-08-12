#ifndef TDRSTYLE_H
#define TDRSTYLE_H

#include <TH1.h>
#include <TF1.h>
#include <TPad.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>

/// define common styles                                                                                                                                       
void SetStyle();
/// create a canveas in common style (default is square)                                                                                                       
TCanvas* MakeCanvas(const char* name, const char* title, int dX=600, int dY=600);
/// define a histogram file (not yet including poly markers, which are used fro data)                                                                          
void InitHist(TH1* hist, const char* xtit, const char* ytit="Number of Entries", int color=kBlack, int style=0);
/// define line style, width and color for signal sample                                                                                                       
void InitSignal(TH1* hist);
/// define common legend style
void SetLegendStyle(TLegend* leg);
/// add CMSPreliminary logo, lumi and channel                                                                                                                  
void CMSPrelim(const char* dataset, const char* channel, double lowX=0.4, double lowY=0.74);

#endif
