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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void BbarB_PurityPlots(const char* filename, const char* prefix,
                     const char* addon = "") {

DreamPlot::SetStyle();

TString add1="1";
TString add2="2";
TString add3="3";
TString add4="4";
TString add5="5";
//Accessing all the Lists in the root file

TList *listTrackCuts=0;
TList *listv0Cuts=0;
TList *listAntiTrackCuts=0;
TList *listAntiv0Cuts=0;
//Accessing all the directory in the root file

  TFile* _file0=TFile::Open(filename);

  TDirectoryFile *dirTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sTrkCutsMC%s", prefix, addon)));
  dirTrackCuts->GetObject(Form("%sTrkCutsMC%s",prefix,addon),listTrackCuts);

  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sv0CutsMC%s", prefix, addon)));
  dirv0Cuts->GetObject(Form("%sv0CutsMC%s",prefix,addon),listv0Cuts);
  TList* listv0CutsMC=(TList*)listv0Cuts->FindObject("v0MonteCarlo");


  TDirectoryFile *dirAntiTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiTrkCutsMC%s", prefix, addon)));
  dirAntiTrackCuts->GetObject(Form("%sAntiTrkCutsMC%s",prefix,addon),listAntiTrackCuts);

  TDirectoryFile *dirAntiv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiv0CutsMC%s", prefix, addon)));
  dirAntiv0Cuts->GetObject(Form("%sAntiv0CutsMC%s",prefix,addon),listAntiv0Cuts);
  TList* listAntiv0CutsMC=(TList*)listAntiv0Cuts->FindObject("v0MonteCarlo");


// Track Cuts plots
std::cout << "------------------------Purity Plots------------------------" << std::endl;
std::cout << std::endl;

  auto* hMCCorrectlyIdentifiedProtons = (TH1F*)listTrackCuts->FindObject("CorrParPt");
  auto* hMCGeneratedProtons = (TH1F*)listTrackCuts->FindObject("GenPartPt");
  auto* hMCIdentifiedProtons = (TH1F*)listTrackCuts->FindObject("IdentPartPt");

  auto* hMCCorrectlyIdentifiedAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("CorrParPt");
  auto* hMCGeneratedAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("GenPartPt");
  auto* hMCIdentifiedAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("IdentPartPt");

  auto* hMCCorrectlyIdentifiedLambda = (TH1F*)listv0CutsMC->FindObject("CorrParPt");
  auto* hMCGeneratedLambda = (TH1F*)listv0CutsMC->FindObject("GenPartPt");
  auto* hMCIdentifiedLambda = (TH1F*)listv0CutsMC->FindObject("IdentPartPt");

  auto* hMCCorrectlyIdentifiedAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("CorrParPt");
  auto* hMCGeneratedAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("GenPartPt");
  auto* hMCIdentifiedAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("IdentPartPt");

  TH1F *PurityProtons = (TH1F*)hMCCorrectlyIdentifiedProtons->Clone();
  PurityProtons->Divide(hMCIdentifiedProtons);
  TH1F *EfficiencyProtons = (TH1F*)hMCCorrectlyIdentifiedProtons->Clone();
  EfficiencyProtons->Divide(hMCGeneratedProtons);

  TH1F *PurityAntiProtons = (TH1F*)hMCCorrectlyIdentifiedAntiProtons->Clone();
  PurityAntiProtons->Divide(hMCIdentifiedAntiProtons);
  TH1F *EfficiencyAntiProtons = (TH1F*)hMCCorrectlyIdentifiedAntiProtons->Clone();
  EfficiencyAntiProtons->Divide(hMCGeneratedAntiProtons);

  TH1F *PurityLambda = (TH1F*)hMCCorrectlyIdentifiedLambda->Clone();
  PurityLambda->Divide(hMCIdentifiedLambda);
  TH1F *EfficiencyLambda = (TH1F*)hMCCorrectlyIdentifiedLambda->Clone();
  EfficiencyLambda->Divide(hMCGeneratedLambda);

  TH1F *PurityAntiLambda = (TH1F*)hMCCorrectlyIdentifiedAntiLambda->Clone();
  PurityAntiLambda->Divide(hMCIdentifiedAntiLambda);
  TH1F *EfficiencyAntiLambda = (TH1F*)hMCCorrectlyIdentifiedAntiLambda->Clone();
  EfficiencyAntiLambda->Divide(hMCGeneratedAntiLambda);

  auto* cProtons = new TCanvas("Can_p","Can_p",0,0,5000,3100);
  cProtons->Divide(2,2);
  PurityProtons->SetTitle(" ; p_{T} [GeV/c]; Purity protons (%)");
  PurityProtons->GetXaxis()->SetRangeUser(0.5,4.05);
  PurityProtons->GetYaxis()->SetRangeUser(0.,1.05);
  EfficiencyProtons->SetTitle("; p_{T} [GeV/c]; Efficiency protons (%)");
  EfficiencyProtons->GetXaxis()->SetRangeUser(0.5,4.05);
  PurityAntiProtons->SetTitle(" ; p_{T} [GeV/c]; Purity Anti-protons (%)");
  PurityAntiProtons->GetXaxis()->SetRangeUser(0.5,4.05);
  PurityAntiProtons->GetYaxis()->SetRangeUser(0.,1.05);
  EfficiencyAntiProtons->SetTitle("; p_{T} [GeV/c]; Efficiency Anti-protons (%)");
  EfficiencyAntiProtons->GetXaxis()->SetRangeUser(0.5,4.05);

  cProtons->cd(1);
  PurityProtons->Draw();
  cProtons->cd(2);
  EfficiencyProtons->Draw();
  cProtons->cd(3);
  PurityAntiProtons->Draw();
  cProtons->cd(4);
  EfficiencyAntiProtons->Draw();

  auto* cLambda = new TCanvas("Can_l","Can_l",0,0,5000,3100);
  cLambda->Divide(2,2);
  PurityLambda->SetTitle(" ; p_{T} [GeV/c]; Purity #Lambda (%)");
  PurityLambda->GetXaxis()->SetRangeUser(0.5,8);
  PurityLambda->GetYaxis()->SetRangeUser(0.,1.05);

  EfficiencyLambda->SetTitle("; p_{T} [GeV/c]; Efficiency #Lambda (%)");
  EfficiencyLambda->GetXaxis()->SetRangeUser(0.5,10.);
  PurityAntiLambda->SetTitle(" ; p_{T} [GeV/c]; Purity Anti-#Lambda (%)");
  PurityAntiLambda->GetXaxis()->SetRangeUser(0.5,10.);
  PurityAntiLambda->GetYaxis()->SetRangeUser(0.,1.05);

  EfficiencyAntiLambda->SetTitle("; p_{T} [GeV/c]; Efficiency Anti-#Lambda (%)");
  EfficiencyAntiLambda->GetXaxis()->SetRangeUser(0.5,10.);


   cLambda->cd(1);
   PurityLambda->Draw();
   cLambda->cd(2);
   EfficiencyLambda->Draw();
   cLambda->cd(3);
   PurityAntiLambda->Draw();
   cLambda->cd(4);
   EfficiencyAntiLambda->Draw();


  std::cout << "--------Primary, Secondary, Material------------------------" << std::endl;
  std::cout << std::endl;

  auto* hMCPrimaryProtons = (TH1F*)listTrackCuts->FindObject("PrimaryPt");
  auto* hMCContProtons = (TH1F*)listTrackCuts->FindObject("ContPt");
  auto* hMCMatProtons = (TH1F*)listTrackCuts->FindObject("MatPt");

  auto* hMCPrimaryAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("PrimaryPt");
  auto* hMCContAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("ContPt");
  auto* hMCMatAntiProtons = (TH1F*)listAntiTrackCuts->FindObject("MatPt");

  auto* hMCPrimaryLambda = (TH1F*)listv0CutsMC->FindObject("PrimaryPt");
  auto* hMCContLambda = (TH1F*)listv0CutsMC->FindObject("ContPt");
  auto* hMCMatLambda = (TH1F*)listv0CutsMC->FindObject("MatPt");

  auto* hMCPrimaryAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("PrimaryPt");
  auto* hMCContAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("ContPt");
  auto* hMCMatAntiLambda = (TH1F*)listAntiv0CutsMC->FindObject("MatPt");
  //
  //ratioCFDataMC[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_ratioCFDataMC0",hCFdata[0]->GetName()));

   TH1F *SumProtonsA = (TH1F*)hMCPrimaryProtons->Clone(Form("%s_clone_SumProtonsA",hMCPrimaryProtons->GetName()));
   SumProtonsA->Add(hMCContProtons);
   TH1F *SumProtonsTot = (TH1F*)SumProtonsA->Clone(Form("%s_clone_SumProtonsTot",SumProtonsA->GetName()));
   SumProtonsTot->Add(hMCMatProtons);

   TH1F *PrimaryProtons = (TH1F*)hMCPrimaryProtons->Clone(Form("%s_clone_PrimaryProtons",hMCPrimaryProtons->GetName()));
   PrimaryProtons->Divide(SumProtonsTot);
   TH1F *ContProtons = (TH1F*)hMCContProtons->Clone(Form("%s_clone_ContProtons",hMCContProtons->GetName()));
   ContProtons->Divide(SumProtonsTot);
   TH1F *MatProtons = (TH1F*)hMCMatProtons->Clone(Form("%s_clone_MatProtons",hMCMatProtons->GetName()));
   MatProtons->Divide(SumProtonsTot);

  // TH1F *EfficiencyProtons = (TH1F*)hMCCorrectlyIdentifiedProtons->Clone();
  // EfficiencyProtons->Divide(hMCGeneratedProtons);
  //
  // TH1F *PurityLambda = (TH1F*)hMCCorrectlyIdentifiedLambda->Clone();
  // PurityLambda->Divide(hMCIdentifiedLambda);
  // TH1F *EfficiencyLambda = (TH1F*)hMCCorrectlyIdentifiedLambda->Clone();
  // EfficiencyLambda->Divide(hMCGeneratedLambda);


  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"The system only dreams in total darkness"<<std::endl;

  // Making PDF file
  TString foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/QA/";

  if(strcmp(addon, add1)==0){
    cProtons->SaveAs(foldername+"PurityEfficiencyProtons_st1.pdf");
    cLambda->SaveAs(foldername+"PurityEfficiencyLambda_st1.pdf");
}
else if(strcmp(addon, add2)==0){
  cProtons->SaveAs(foldername+"PurityEfficiencyProtons_st2.pdf");
  cLambda->SaveAs(foldername+"PurityEfficiencyLambda_st2.pdf");
}
else if(strcmp(addon, add3)==0){
  cProtons->SaveAs(foldername+"PurityEfficiencyProtons_st3.pdf");
  cLambda->SaveAs(foldername+"PurityEfficiencyLambda_st3.pdf");
}
else if(strcmp(addon, add4)==0){
  cProtons->SaveAs(foldername+"PurityEfficiencyProtons_full.pdf");
  cLambda->SaveAs(foldername+"PurityEfficiencyLambda_full.pdf");
}

}
