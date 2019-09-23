#include "TROOT.h"
#include "ReadDreamFile.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include "TPad.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TApplication.h"

void PlotsCorrelationBBar (const char* STflag = "") {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);


TString add1="1";
TString add2="2";
TString add3="3";
TString add4="4";
TString add5="5";
TString addon;
//Accessing all the directory in the root file
  TFile* file0data[4];
  TFile* file0MC[4];

  bool whichfolder = true;// true ->CorrectBin folder

  if(!whichfolder){
  if(strcmp(STflag, add1)==0)
  {
    addon=add1;
  std::cout << "Sphericity [0.,0.3] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st1.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st1.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st1.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st1.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st1.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st1.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st1.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st1.root");
  }
  else if(strcmp(STflag, add2)==0)
  {
    addon=add2;

  std::cout << "Sphericity [0.3,0.7] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st2.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st2.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st2.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st2.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st2.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st2.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st2.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st2.root");
  }
  else if(strcmp(STflag, add3)==0)
  {
    addon=add3;

  std::cout << "Sphericity [0.7,1.] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st3.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st3.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st3.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st3.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st3.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st3.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st3.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st3.root");
  }
  else if(strcmp(STflag, add4)==0)
  {
    addon=add4;

  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_full.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_full.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_full.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_full.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_full.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_full.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_full.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_full.root");
  }
}
else if(whichfolder){
  if(strcmp(STflag, add1)==0)
  {
    addon=add1;

  std::cout << "Sphericity [0.,0.3] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st1.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st1.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st1.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st1.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st1.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st1.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st1.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st1.root");
  }
  else if(strcmp(STflag, add2)==0)
  {
    addon=add2;

  std::cout << "Sphericity [0.3,0.7] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st2.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st2.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st2.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st2.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st2.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st2.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st2.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st2.root");
  }
  else if(strcmp(STflag, add3)==0)
  {
    addon=add3;

  std::cout << "Sphericity [0.7,1.] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st3.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st3.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st3.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st3.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st3.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st3.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st3.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st3.root");
  }
  else if(strcmp(STflag, add4)==0)
  {
    addon=add4;

  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_full.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_full.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_full.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_full.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_full.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_full.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_full.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_full.root");
  }
}

//REMEMBER THAT WHEN ACCESSING DIRECTLY FROM FILE USE FindObjectAny!!!!
   TH1F* hCFdata[3];
   TH1F* hCFmc[3];

   hCFdata[0] = (TH1F*)(file0data[0]->FindObjectAny("hCk_ReweightedpApData_0"));
   hCFdata[0]->Sumw2();
   hCFmc[0] = (TH1F*)(file0MC[0]->FindObjectAny("hCk_ReweightedpApMC_0"));
   hCFmc[0]->Sumw2();

   hCFdata[1] = (TH1F*)(file0data[1]->FindObjectAny("hCk_ReweightedpALData_0"));
   hCFdata[1]->Sumw2();
   hCFmc[1] = (TH1F*)(file0MC[1]->FindObjectAny("hCk_ReweightedpALMC_0"));
   hCFmc[1]->Sumw2();
   hCFdata[2] = (TH1F*)(file0data[2]->FindObjectAny("hCk_ReweightedLALData_0"));
   hCFdata[2]->Sumw2();
   hCFmc[2] = (TH1F*)(file0MC[2]->FindObjectAny("hCk_ReweightedLALMC_0"));
   hCFmc[2]->Sumw2();

   hCFdata[0]->SetLineColor(kBlack);
   hCFdata[1]->SetLineColor(kRed);
   hCFdata[2]->SetLineColor(kBlue);

   hCFdata[0]->SetLineWidth(2);
   hCFdata[1]->SetLineWidth(2);
   hCFdata[2]->SetLineWidth(2);

   // hCFdata[0]->SetMarkerStyle(kFullCircle);
   // hCFdata[1]->SetMarkerStyle(kFullSquare);
   // hCFdata[2]->SetMarkerStyle(kFullDiamond);
   // hCFdata[0]->SetMarkerSize(0.5);
   // hCFdata[1]->SetMarkerSize(0.5);
   // hCFdata[2]->SetMarkerSize(0.5);

   hCFmc[0]->SetLineColor(kOrange+1);
   hCFmc[1]->SetLineColor(kOrange+1);
   hCFmc[2]->SetLineColor(kOrange+1);

   hCFmc[0]->SetLineWidth(2);
   hCFmc[1]->SetLineWidth(2);
   hCFmc[2]->SetLineWidth(2);

   hCFmc[0]->SetMarkerStyle(kFullCircle);
   hCFmc[1]->SetMarkerStyle(kFullSquare);
   hCFmc[2]->SetMarkerStyle(kFullDiamond);
   hCFmc[0]->SetMarkerSize(0.4);
   hCFmc[1]->SetMarkerSize(0.4);
   hCFmc[2]->SetMarkerSize(0.4);

   TString foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/CFPlots/";

   auto* cPlotpAp_data = new TCanvas("cPlotpAp_data","cPlotpAp_data",0,0,5000,5000);
   hCFdata[0]->GetYaxis()->SetTitle("C(k *)");
   hCFdata[0]->GetXaxis()->SetTitle("k* [GeV/c]");
   hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.5);
   hCFdata[0]->Draw();
   TLegend *legpAp = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
 legpAp->SetFillStyle(0);
 gStyle->SetLegendTextSize(0.02);
 legpAp->AddEntry(hCFdata[0],"p + #bar{p}");
 legpAp->SetLineColor(0);
 legpAp->Draw("same");
 cPlotpAp_data->SaveAs(TString::Format(foldernameplot+"pAp_%s.pdf",addon.Data()));


 auto* cPlotpAp_datazoom = new TCanvas("cPlotpAp_datazoom","cPlotpAp_datazoom",0,0,5000,5000);
 hCFdata[0]->GetYaxis()->SetTitle("C(k *)");
 hCFdata[0]->GetXaxis()->SetTitle("k * [GeV/c]");
 hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.25);
 hCFdata[0]->Draw();
 cPlotpAp_datazoom->SaveAs(TString::Format(foldernameplot+"pApzoom_%s.pdf",addon.Data()));


 auto* cPlotpAp_dataMC = new TCanvas("cPlotpAp_dataMC","cPlotpAp_dataMC",0,0,5000,5000);
 hCFdata[0]->GetYaxis()->SetTitle("C(k^{*})");
 hCFdata[0]->GetXaxis()->SetTitle("k * [GeV/c]");
 hCFmc[0]->GetXaxis()->SetRangeUser(0.,0.5);
 hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.5);
 hCFmc[0]->Draw();
 hCFdata[0]->Draw("same");
 TLegend *legpApMC = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legpApMC->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legpApMC->AddEntry(hCFdata[0]," data p + #bar{p}");
legpApMC->AddEntry(hCFmc[0]," MC p + #bar{p}");
legpApMC->SetLineColor(0);
legpApMC->Draw("same");
cPlotpAp_dataMC->SaveAs(TString::Format(foldernameplot+"pApMCdata_%s.pdf",addon.Data()));


 auto* cPlotpAL_data = new TCanvas("cPlotpAL_data","cPlotpAL_data",0,0,5000,5000);
 hCFdata[1]->GetYaxis()->SetTitle("C(k *)");
 hCFdata[1]->GetXaxis()->SetTitle("k * [GeV/c]");
 hCFdata[1]->GetXaxis()->SetRangeUser(0.,0.5);
 hCFdata[1]->Draw();
 TLegend *legpAL = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legpAL->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legpAL->AddEntry(hCFdata[1],"p + #bar{#Lambda}");
legpAL->SetLineColor(0);
legpAL->Draw("same");
cPlotpAL_data->SaveAs(TString::Format(foldernameplot+"pAL_%s.pdf",addon.Data()));

auto* cPlotpAL_datazoom = new TCanvas("cPlotpAL_datazoom","cPlotpAL_datazoom",0,0,5000,5000);
hCFdata[1]->GetYaxis()->SetTitle("C(k *)");
hCFdata[1]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[1]->GetXaxis()->SetRangeUser(0.,0.25);
hCFdata[1]->Draw();
cPlotpAL_datazoom->SaveAs(TString::Format(foldernameplot+"pALzoom_%s.pdf",addon.Data()));

auto* cPlotpAL_dataMC = new TCanvas("cPlotpAL_dataMC","cPlotpAL_dataMC",0,0,5000,5000);
hCFdata[1]->GetYaxis()->SetTitle("C(k *)");
hCFdata[1]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFmc[1]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[1]->GetXaxis()->SetRangeUser(0.,0.5);
hCFmc[1]->Draw();
hCFdata[1]->Draw("same");
TLegend *legpALMC = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legpALMC->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legpALMC->AddEntry(hCFdata[1]," data p + #bar{#Lambda}");
legpALMC->AddEntry(hCFmc[1]," MC p + #bar{#Lambda}");
legpALMC->SetLineColor(0);
legpALMC->Draw("same");
cPlotpAL_dataMC->SaveAs(TString::Format(foldernameplot+"pALMCdata_%s.pdf",addon.Data()));


auto* cPlotLAL_data = new TCanvas("cPlotLAL_data","cPlotLAL_data",0,0,5000,5000);
hCFdata[2]->GetYaxis()->SetTitle("C(k *)");
hCFdata[2]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[2]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[2]->Draw();
TLegend *legLAL = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legLAL->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legLAL->AddEntry(hCFdata[2],"#Lambda + #bar{#Lambda}");
legLAL->SetLineColor(0);
legLAL->Draw("same");
cPlotLAL_data->SaveAs(TString::Format(foldernameplot+"LAL_%s.pdf",addon.Data()));


auto* cPlotLAL_datazoom = new TCanvas("cPlotLAL_datazoom","cPlotLAL_datazoom",0,0,5000,5000);
hCFdata[2]->GetYaxis()->SetTitle("C(k *)");
hCFdata[2]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[2]->GetXaxis()->SetRangeUser(0.,0.25);
hCFdata[2]->Draw();
cPlotLAL_datazoom->SaveAs(TString::Format(foldernameplot+"LALzoom_%s.pdf",addon.Data()));

auto* cPlotLAL_dataMC = new TCanvas("cPlotLAL_dataMC","cPlotLAL_dataMC",0,0,5000,5000);
hCFdata[2]->GetYaxis()->SetTitle("C(k *)");
hCFdata[2]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFmc[2]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[2]->GetXaxis()->SetRangeUser(0.,0.5);
hCFmc[2]->Draw();
hCFdata[2]->Draw("same");
TLegend *legLALMC = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legLALMC->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legLALMC->AddEntry(hCFdata[2]," data #Lambda + #bar{#Lambda}");
legLALMC->AddEntry(hCFmc[2]," MC #Lambda + #bar{#Lambda}");
legLALMC->SetLineColor(0);
legLALMC->Draw("same");
cPlotLAL_dataMC->SaveAs(TString::Format(foldernameplot+"LALMCdata_%s.pdf",addon.Data()));


auto* cPlotAll_data = new TCanvas("cPlotAll_data","cPlotAll_data",0,0,5000,5000);
hCFdata[0]->GetYaxis()->SetTitle("C(k *)");
hCFdata[0]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[1]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[2]->GetXaxis()->SetRangeUser(0.,0.5);
hCFdata[0]->GetYaxis()->SetRangeUser(0.1,3.42);
hCFdata[0]->Draw();
hCFdata[1]->Draw("same");
hCFdata[2]->Draw("same");
TLegend *legall = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legall->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legall->AddEntry(hCFdata[0]," p + #bar{p}");
legall->AddEntry(hCFdata[1]," p + #bar{#Lambda}");
legall->AddEntry(hCFdata[2]," #Lambda + #bar{#Lambda}");
legall->SetLineColor(0);
legall->Draw("same");
cPlotAll_data->SaveAs(TString::Format(foldernameplot+"All_data_%s.pdf",addon.Data()));


auto* cPlotAll_datazoom = new TCanvas("cPlotAll_datazoom","cPlotAll_datazoom",0,0,5000,5000);
hCFdata[0]->GetYaxis()->SetTitle("C(k *)");
hCFdata[0]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.25);
hCFdata[0]->GetYaxis()->SetRangeUser(0.1,3.42);

hCFdata[0]->Draw();
hCFdata[1]->Draw("same");
hCFdata[2]->Draw("same");
TLegend *legallz = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legallz->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legallz->AddEntry(hCFdata[0]," p + #bar{p}");
legallz->AddEntry(hCFdata[1]," p + #bar{#Lambda}");
legallz->AddEntry(hCFdata[2]," #Lambda + #bar{#Lambda}");
legallz->SetLineColor(0);
legallz->Draw("same");
cPlotAll_datazoom->SaveAs(TString::Format(foldernameplot+"Allzoom_data_%s.pdf",addon.Data()));


auto* cPlotAll_dataLong = new TCanvas("cPlotAll_dataLong","cPlotAll_dataLong",0,0,5000,5000);
hCFdata[0]->GetYaxis()->SetTitle("C(k *)");
hCFdata[0]->GetXaxis()->SetTitle("k * [GeV/c]");
hCFdata[0]->GetXaxis()->SetRangeUser(0.,1.5);
hCFdata[1]->GetXaxis()->SetRangeUser(0.,1.5);
hCFdata[2]->GetXaxis()->SetRangeUser(0.,1.5);
hCFdata[0]->GetYaxis()->SetRangeUser(0.1,3.42);
hCFdata[0]->Draw();
hCFdata[1]->Draw("same");
hCFdata[2]->Draw("same");
TLegend *legallL = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
legallL->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
legallL->AddEntry(hCFdata[0]," p + #bar{p}");
legallL->AddEntry(hCFdata[1]," p + #bar{#Lambda}");
legallL->AddEntry(hCFdata[2]," #Lambda + #bar{#Lambda}");
legallL->SetLineColor(0);
legallL->Draw("same");
cPlotAll_dataLong->SaveAs(TString::Format(foldernameplot+"All_dataLong_%s.pdf",addon.Data()));







}
