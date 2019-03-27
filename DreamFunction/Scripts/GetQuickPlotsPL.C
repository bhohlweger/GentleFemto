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

void GetQuickPlotsPL() {


DreamPlot::SetStyle();


  TFile* file0data=TFile::Open("/Users/Valentina/HMpp13TeV/AnalysisData/ClosePairRej/SelectedPairs/CFOutput_pL.root");
  TFile* file0mc=TFile::Open("/Users/Valentina/HMpp13TeV/Sim_MB_5/CFOutput_pL.root");

  TList* listPairDistdata[2];
  TList* listDistdata[2];
  TH1F* hCFdata[2];
  TH1F* hCFd[2];
  TH1F* hCFm[2];


  hCFd[0] = (TH1F*)(file0data->FindObjectAny("hCk_ReweightedMeV_0"));
  hCFd[0]->Sumw2();

  hCFm[0] = (TH1F*)(file0mc->FindObjectAny("hCk_ReweightedMeV_0"));
  hCFm[0]->Sumw2();

  listPairDistdata[0]=(TList*)(file0data->FindObjectAny("PairDist"));
  listPairDistdata[1]=(TList*)(file0data->FindObjectAny("AntiPairDist"));

  listDistdata[0]=(TList*)listPairDistdata[0]->FindObject("PairReweighted");
  listDistdata[1]=(TList*)listPairDistdata[1]->FindObject("PairReweighted");

  hCFdata[0]=(TH1F*)listDistdata[0]->FindObject("CFDist_Particle0_Particle2_clone_Shifted_FixShifted_Rebinned_4_Reweighted");
  hCFdata[1]=(TH1F*)listDistdata[1]->FindObject("CFDist_Particle1_Particle3_clone_Shifted_FixShifted_Rebinned_4_Reweighted");

  TH1F *ratio;


  ratio = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_ratioCFData0",hCFdata[0]->GetName()));
  ratio->Divide(hCFdata[1]);

  TCanvas* b = new TCanvas();
  ratio->SetTitle("; k* (GeV/#it{c}); #it{C}_{p-#Lambda}(k*)/#it{C}_{#bar{p}-#bar{#Lambda}}(k*)");
  ratio->GetYaxis()->SetTitleOffset(1.5);
  ratio->GetXaxis()->SetRangeUser(0.,0.5);
  ratio->GetYaxis()->SetRangeUser(0.9,1.03);

  ratio->SetMarkerStyle(21);
  ratio->SetMarkerColor(1);
  ratio->SetMarkerSize(0.6);
  ratio->SetLineColor(1);
  ratio->SetLineWidth(2);
  ratio->Draw();



  TCanvas* a = new TCanvas();
  hCFdata[0]->SetTitle("; k* (GeV/#it{c}); #it{C}(k*) ");
  hCFdata[0]->GetYaxis()->SetTitleOffset(1.5);
  hCFdata[0]->GetXaxis()->SetRangeUser(0.,0.5);
  hCFdata[0]->GetYaxis()->SetRangeUser(0.9,2.);
  hCFdata[0]->SetMarkerStyle(8);
  hCFdata[0]->SetMarkerColor(1);
  hCFdata[0]->SetLineColor(1);
  hCFdata[0]->SetLineWidth(2);
  hCFdata[0]->Draw();

  hCFdata[1]->SetMarkerStyle(21);
  hCFdata[1]->SetMarkerSize(0.6);
  hCFdata[1]->SetMarkerColor(6);
  hCFdata[1]->SetLineColor(6);
  hCFdata[1]->SetLineWidth(2);
  hCFdata[1]->Draw("same");

  TCanvas* c = new TCanvas();
  hCFd[0]->SetTitle("; k* (MeV/#it{c}); #it{C}(k*) ");
  hCFd[0]->GetYaxis()->SetTitleOffset(1.5);
  hCFd[0]->GetXaxis()->SetRangeUser(0.,1200);
  hCFd[0]->GetYaxis()->SetRangeUser(0.5,2.);
  hCFd[0]->SetMarkerStyle(8);
  hCFd[0]->SetMarkerColor(1);
  hCFd[0]->SetLineColor(1);
  hCFd[0]->SetLineWidth(2);
  hCFd[0]->Draw();

  hCFm[0]->SetMarkerStyle(21);
  hCFm[0]->SetMarkerSize(0.6);
  hCFm[0]->SetMarkerColor(4);
  hCFm[0]->SetLineColor(4);
  hCFm[0]->SetLineWidth(2);
  hCFm[0]->Draw("same");

  auto* leg0= new TLegend(0.65,0.65,0.8,0.8);
  leg0->AddEntry(hCFd[0], "Data", "l");
  leg0->AddEntry(hCFm[0], "Pythia 8", "l");
  leg0->Draw("same");


  TString foldernameplot = "/Users/Valentina/cernbox/AN/PLOTS_AN/";

  a->SaveAs(foldernameplot+"Comparison_pL.pdf");
  b->SaveAs(foldernameplot+"Ratio_pL.pdf");
  c->SaveAs(foldernameplot+"pLDataVsMC.pdf");



}
