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

void GetComparisonMC(const char* filename1, const char* filename2, const char* prefix,
                     const char* addon = "") {


DreamPlot::SetStyle();

TString add1="1";
TString add2="2";
TString add3="3";
TString add4="4";
TString add5="5";

  TFile* _file0mc1=TFile::Open(filename1);
  TFile* _file0mc2=TFile::Open(filename2);

  TDirectoryFile *dirResultsmc1=(TDirectoryFile*)(_file0mc1->FindObjectAny(Form("%sResults%s", prefix, addon)));
  TDirectoryFile *dirResultsmc2=(TDirectoryFile*)(_file0mc2->FindObjectAny(Form("%sResults%s", prefix, addon)));

  TList *Resultsmc1;
  TList *Resultsmc2;
  dirResultsmc1->GetObject(Form("%sResults%s", prefix, addon),Resultsmc1);
  dirResultsmc2->GetObject(Form("%sResults%s", prefix, addon),Resultsmc2);
//Proton-AntiProton----------------------------------------------------
  TList* tmpFoldermc1=(TList*)Resultsmc1->FindObject("Particle0_Particle1");
  TList* tmpFoldermc2=(TList*)Resultsmc2->FindObject("Particle0_Particle1");

  TH1F* SEMC1[4];
  TH1F* SEMC2[4];
  TH1F* MEMC1[4];
  TH1F* MEMC2[4];

  SEMC1[0] = (TH1F*)tmpFoldermc1->FindObject("SEDist_Particle0_Particle1");
  SEMC1[0]->Sumw2();
  MEMC1[0]= (TH1F*)tmpFoldermc1->FindObject("MEDist_Particle0_Particle1");
  MEMC1[0]->Sumw2();
  SEMC2[0] = (TH1F*)tmpFoldermc2->FindObject("SEDist_Particle0_Particle1");
  SEMC2[0]->Sumw2();
  MEMC2[0]= (TH1F*)tmpFoldermc2->FindObject("MEDist_Particle0_Particle1");
  MEMC2[0]->Sumw2();

  double entriesSEMC1 = SEMC1[0]->GetEntries();
  double entriesMEMC1 = MEMC1[0]->GetEntries();

  double entriesSEMC2 = SEMC2[0]->GetEntries();
  double entriesMEMC2 = MEMC2[0]->GetEntries();


  TH1F *scaleSEMC1[4];
  TH1F *scaleSEMC2[4];

  TH1F *scaleMEMC1[4];
  TH1F *scaleMEMC2[4];

  std::cout<<"SE Distributions for Proton-Antiproton"<<std::endl;

  scaleSEMC1[0] = (TH1F*)SEMC1[0]->Clone(Form("%s_clone_scaleSEMC10",SEMC1[0]->GetName()));
  scaleSEMC1[0]->Scale(1/entriesSEMC1);
  scaleSEMC2[0] = (TH1F*)SEMC2[0]->Clone(Form("%s_clone_scaleSEMC20",SEMC2[0]->GetName()));
  scaleSEMC2[0]->Scale(1/entriesSEMC2);


  scaleMEMC1[0] = (TH1F*)MEMC1[0]->Clone(Form("%s_clone_scaleMEMC10",MEMC1[0]->GetName()));
  scaleMEMC1[0]->Scale(1/entriesMEMC1);
  scaleMEMC2[0] = (TH1F*)MEMC2[0]->Clone(Form("%s_clone_scaleMEMC20",MEMC2[0]->GetName()));
  scaleMEMC2[0]->Scale(1/entriesMEMC2);

  TH1F *TmpscaleSEMC1[4];
  TH1F *TmpscaleSEMC2[4];

  TH1F *TmpscaleMEMC1[4];
  TH1F *TmpscaleMEMC2[4];

  TmpscaleSEMC1[0] = (TH1F*)scaleSEMC1[0]->Clone(Form("%s_clone_TmpscaleSEMC10",scaleSEMC1[0]->GetName()));
  TmpscaleSEMC2[0] = (TH1F*)scaleSEMC2[0]->Clone(Form("%s_clone_TmpscaleSEMC20",scaleSEMC2[0]->GetName()));

  TmpscaleMEMC1[0] = (TH1F*)scaleMEMC1[0]->Clone(Form("%s_clone_TmpscaleMEMC10",scaleMEMC1[0]->GetName()));
  TmpscaleMEMC2[0] = (TH1F*)scaleMEMC2[0]->Clone(Form("%s_clone_TmpscaleMEMC20",scaleMEMC2[0]->GetName()));


  TCanvas* Can1 = new TCanvas("Can1");
  Can1->Divide(2, 2);
  Can1->cd(1);
  scaleSEMC1[0]->SetTitle("; k* (GeV/#it{c}); #it{SE}(k*) ");
  scaleSEMC1[0]->GetYaxis()->SetTitleOffset(1.5);

  scaleSEMC1[0]->SetMarkerStyle(1);
  scaleSEMC1[0]->SetMarkerColor(2);
  scaleSEMC1[0]->SetLineColor(2);

  scaleSEMC1[0]->Draw();
  scaleSEMC2[0]->SetMarkerStyle(1);
  scaleSEMC2[0]->SetMarkerColor(1);
  scaleSEMC2[0]->SetLineColor(1);

  scaleSEMC2[0]->Draw("same");
  auto* leg1= new TLegend(0.55,0.65,0.75,0.8);
  leg1->AddEntry(scaleSEMC1[0], "Pythia 8", "l");//"l" sets the legend as lines
  leg1->AddEntry(scaleSEMC2[0], "EPOS LHC", "l");
  leg1->Draw("same");

  Can1->cd(2);
  scaleMEMC1[0]->SetTitle("; k* (GeV/#it{c}); #it{ME}(k*)");
  scaleMEMC1[0]->GetYaxis()->SetTitleOffset(1.5);

  scaleMEMC1[0]->SetMarkerStyle(1);
  scaleMEMC1[0]->SetMarkerColor(2);
  scaleMEMC1[0]->SetLineColor(2);

  scaleMEMC2[0]->SetMarkerStyle(1);
  scaleMEMC2[0]->SetMarkerColor(1);
  scaleMEMC2[0]->SetLineColor(1);

  scaleMEMC1[0]->Draw();
  // DreamPlot::SetStyleHisto(scaleMEMC1[0], 1, 2);
  // DreamPlot::SetStyleHisto(scaleMEMC2[0], 2, 1);
  scaleMEMC2[0]->Draw("same");

//Comparison SE ME within same MC production
   Can1->cd(3);
   TmpscaleSEMC1[0]->SetTitle("; k* (GeV/#it{c}); ");
   TmpscaleSEMC1[0]->SetMarkerStyle(1);
   TmpscaleSEMC1[0]->SetMarkerColor(3);
   TmpscaleSEMC1[0]->SetLineColor(3);

   TmpscaleMEMC1[0]->SetMarkerStyle(1);
   TmpscaleMEMC1[0]->SetMarkerColor(4);
   TmpscaleMEMC1[0]->SetLineColor(4);

   TmpscaleSEMC1[0]->Draw();
   TmpscaleMEMC1[0]->Draw("same");


  // DreamPlot::SetStyleHisto(scaleSEMC1[0], 1, 2);
  // DreamPlot::SetStyleHisto(scaleMEMC1[0], 1, 1);
   auto* leg= new TLegend(0.55,0.65,0.75,0.8);
   leg->AddEntry(TmpscaleSEMC1[0], "Same event", "l");//"l" sets the legend as lines
   leg->AddEntry(TmpscaleMEMC1[0], "Mixed event", "l");
   TLatex label1;
   label1.SetNDC(kTRUE);
   label1.DrawLatex(0.55, 0.4, TString::Format("PYTHIA"));
   leg->Draw("same");
  //
  //
   Can1->cd(4);
   TmpscaleMEMC2[0]->SetTitle("; k* (GeV/#it{c}); ");

   TmpscaleMEMC2[0]->SetMarkerStyle(1);
   TmpscaleMEMC2[0]->SetMarkerColor(4);
   TmpscaleMEMC2[0]->SetLineColor(4);

   TmpscaleSEMC2[0]->SetMarkerStyle(1);
   TmpscaleSEMC2[0]->SetMarkerColor(3);
   TmpscaleSEMC2[0]->SetLineColor(3);
   TmpscaleMEMC2[0]->Draw();
   TmpscaleSEMC2[0]->Draw("same");

  TLatex label2;
  label2.SetNDC(kTRUE);
  label2.DrawLatex(0.55, 0.4, TString::Format("EPOS"));
  leg->Draw("same");

  TString foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/ComparisonMC/Plot/";

  if(strcmp(addon, add4)==0){
  Can1->SaveAs(foldernameplot+"SEME_PYTHIA_EPOS_full.pdf");
  // Can2->SaveAs(foldernameplot+"SEME_singleComp_full.pdf");

  }
   // //Lambda-AntiLambda----------------------------------------------------
   // DreamDist* LAL_pairdata = new DreamDist();
   // tmpFolderdata=(TList*)Resultsdata->FindObject("Particle2_Particle3");
   // TH1F* SEDataLAL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle2_Particle3");
   // TH1F* MEDataLAL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle2_Particle3");
   // TH2F* SEDataMultLAL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle2_Particle3");
   // TH2F* MEDataMultLAL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle2_Particle3");
   //  LAL_pairdata->SetSEDist(SEDataLAL_pair , "_Shifted");
   //  LAL_pairdata->SetSEMultDist(SEDataMultLAL_pair, "_Shifted");
   //  LAL_pairdata->SetMEDist(MEDataLAL_pair, "_Shifted");
   //  LAL_pairdata->SetMEMultDist(MEDataMultLAL_pair, "_Shifted");
   //
   // DreamDist* LAL_pairmc = new DreamDist();
   // tmpFoldermc=(TList*)Resultsmc->FindObject("Particle2_Particle3");
   // TH1F* SEMCLAL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle2_Particle3");
   // TH1F* MEMCLAL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle2_Particle3");
   // TH2F* SEMCMultLAL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle2_Particle3");
   // TH2F* MEMCMultLAL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle2_Particle3");
   //  LAL_pairmc->SetSEDist(SEMCLAL_pair , "_Shifted");
   //  LAL_pairmc->SetSEMultDist(SEMCMultLAL_pair, "_Shifted");
   //  LAL_pairmc->SetMEDist(MEMCLAL_pair, "_Shifted");
   //  LAL_pairmc->SetMEMultDist(MEMCMultLAL_pair, "_Shifted");
   //
   //  //Proton-AntiLambda----------------------------------------------------
   //  DreamDist* pAL_pairdata = new DreamDist();
   //  tmpFolderdata=(TList*)Resultsdata->FindObject("Particle0_Particle3");
   //  TH1F* SEDatapAL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle0_Particle3");
   //  TH1F* MEDatapAL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle0_Particle3");
   //  TH2F* SEDataMultpAL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle0_Particle3");
   //  TH2F* MEDataMultpAL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle0_Particle3");
   //   pAL_pairdata->SetSEDist(SEDatapAL_pair , "_Shifted");
   //   pAL_pairdata->SetSEMultDist(SEDataMultpAL_pair, "_Shifted");
   //   pAL_pairdata->SetMEDist(MEDatapAL_pair, "_Shifted");
   //   pAL_pairdata->SetMEMultDist(MEDataMultpAL_pair, "_Shifted");
   //
   //   DreamDist* pAL_pairmc = new DreamDist();
   //   tmpFoldermc=(TList*)Resultsmc->FindObject("Particle0_Particle3");
   //   TH1F* SEMCpAL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle0_Particle3");
   //   TH1F* MEMCpAL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle0_Particle3");
   //   TH2F* SEMCMultpAL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle0_Particle3");
   //   TH2F* MEMCMultpAL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle0_Particle3");
   //    pAL_pairmc->SetSEDist(SEMCpAL_pair , "_Shifted");
   //    pAL_pairmc->SetSEMultDist(SEMCMultpAL_pair, "_Shifted");
   //    pAL_pairmc->SetMEDist(MEMCpAL_pair, "_Shifted");
   //    pAL_pairmc->SetMEMultDist(MEMCMultpAL_pair, "_Shifted");
   //
   //   //AntiProton-Lambda----------------------------------------------------
   //   DreamDist* ApL_pairdata = new DreamDist();
   //   tmpFolderdata=(TList*)Resultsdata->FindObject("Particle1_Particle2");
   //   TH1F* SEDataApL_pair = (TH1F*)tmpFolderdata->FindObject("SEDist_Particle1_Particle2");
   //   TH1F* MEDataApL_pair = (TH1F*)tmpFolderdata->FindObject("MEDist_Particle1_Particle2");
   //   TH2F* SEDataMultApL_pair = (TH2F*)tmpFolderdata->FindObject("SEMultDist_Particle1_Particle2");
   //   TH2F* MEDataMultApL_pair = (TH2F*)tmpFolderdata->FindObject("MEMultDist_Particle1_Particle2");
   //    ApL_pairdata->SetSEDist(SEDataApL_pair , "_Shifted");
   //    ApL_pairdata->SetSEMultDist(SEDataMultApL_pair, "_Shifted");
   //    ApL_pairdata->SetMEDist(MEDataApL_pair, "_Shifted");
   //    ApL_pairdata->SetMEMultDist(MEDataMultApL_pair, "_Shifted");
   //
   //    DreamDist* ApL_pairmc = new DreamDist();
   //    tmpFoldermc=(TList*)Resultsmc->FindObject("Particle1_Particle2");
   //    TH1F* SEMCApL_pair = (TH1F*)tmpFoldermc->FindObject("SEDist_Particle1_Particle2");
   //    TH1F* MEMCApL_pair = (TH1F*)tmpFoldermc->FindObject("MEDist_Particle1_Particle2");
   //    TH2F* SEMCMultApL_pair = (TH2F*)tmpFoldermc->FindObject("SEMultDist_Particle1_Particle2");
   //    TH2F* MEMCMultApL_pair = (TH2F*)tmpFoldermc->FindObject("MEMultDist_Particle1_Particle2");
   //     ApL_pairmc->SetSEDist(SEMCApL_pair , "_Shifted");
   //     ApL_pairmc->SetSEMultDist(SEMCMultApL_pair, "_Shifted");
   //     ApL_pairmc->SetMEDist(MEMCApL_pair, "_Shifted");
   //     ApL_pairmc->SetMEMultDist(MEMCMultApL_pair, "_Shifted");



}
