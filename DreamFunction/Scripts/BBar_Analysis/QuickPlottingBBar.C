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

int main(int argc, char* argv[]) {
  const char* date = argv[1];
  const char* maxkstar = (argv[2]) ? argv[2] : "";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  // gStyle->SetLegendTextSize(0.8);


  DreamPlot::SetStyle();

  TString add1="1";
  TString add2="2";
  TString add3="3";
  TString add4="4";
  TString add5="5";
  TString add6="6";
  TString add8="8";


  TString add01="01";
  TString add02="02";
  TString add03="03";
  TString add04="04";
  TString add05="05";
  TString add06="06";
  TString add08="08";


  double charmask = atof (maxkstar);

  std::cout<<"charmask = "<<charmask<<std::endl;

  gSystem->Exec(TString::Format("mkdir ../../../BBbar/GentleFemto_Output/CFPlots/%s",date));
  TString OutputFolder = TString::Format("../../../BBbar/GentleFemto_Output/CFPlots/%s/",date);
  TFile* fOut = new TFile(TString::Format(OutputFolder+"ComparisonHM_StS0.root",date),"update");

//Accessing all the directory in the root file
  TFile* fileHMdataS0_01[4];
  TFile* fileHMdataS0_02[4];
  TFile* fileHMdataS0_03[4];
  TFile* fileHMdataS0_04[4];
  TFile* fileHMdataS0_05[4];
  TFile* fileHMdataS0_06[4];
  TFile* fileHMdataS0_08[4];
  TString s0_label_01;TString s0_label_02;TString s0_label_03;TString s0_label_04;
  TString s0_label_05;TString s0_label_06;TString s0_label_08;

  TFile* fileHMdataST_01[4];
  TFile* fileHMdataST_02[4];
  TFile* fileHMdataST_03[4];
  TFile* fileHMdataST_04[4];
  TFile* fileHMdataST_05[4];
  TFile* fileHMdataST_06[4];
  TFile* fileHMdataST_08[4];
  TString st_label_01;TString st_label_02;TString st_label_03;TString st_label_04;
  TString st_label_05;TString st_label_06;TString st_label_08;

  std::cout << "Reading SpherIcity/SpherOcity [0.,0.3] \n";
   st_label_01 = "0 < S_{T} < 0.3";
   s0_label_01 = "0 < S_{0} < 0.3";
  fileHMdataS0_01[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st1.root");
  fileHMdataST_01[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st1.root");
  fileHMdataS0_01[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st1.root");
  fileHMdataST_01[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st1.root");
  fileHMdataS0_01[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st1.root");
  fileHMdataST_01[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st1.root");

  std::cout << "Reading SpherIcity/SpherOcity [0.3,0.7] \n";
  st_label_02 = "0.3 < S_{T} < 0.7";
  s0_label_02 = "0.3 < S_{0} < 0.7";
  fileHMdataS0_02[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st2.root");
  fileHMdataST_02[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st2.root");
  fileHMdataS0_02[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st2.root");
  fileHMdataST_02[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st2.root");
  fileHMdataS0_02[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st2.root");
  fileHMdataST_02[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st2.root");

  std::cout << "Reading SpherIcity/SpherOcity [0.7,1] \n";
  st_label_03 = "0.7 < S_{T} < 1";
  s0_label_03 = "0.7 < S_{0} < 1";
  fileHMdataS0_03[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st3.root");
  fileHMdataST_03[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st3.root");
  fileHMdataS0_03[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st3.root");
  fileHMdataST_03[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st3.root");
  fileHMdataS0_03[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st3.root");
  fileHMdataST_03[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st3.root");

  std::cout << "Reading SpherIcity/SpherOcity [0,1] \n";
  st_label_04 = "0 < S_{T} < 1";
  s0_label_04 = "0 < S_{0} < 1";
  fileHMdataS0_04[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_full.root");
  fileHMdataST_04[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_full.root");
  fileHMdataS0_04[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_full.root");
  fileHMdataST_04[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_full.root");
  fileHMdataS0_04[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_full.root");
  fileHMdataST_04[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_full.root");

  std::cout << "Reading SpherIcity/SpherOcity [0.8,1] \n";
  st_label_05 = "0.8 < S_{T} < 1";
  s0_label_05 = "0.8 < S_{0} < 1";
  fileHMdataS0_05[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st5.root");
  fileHMdataST_05[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAp_App_st5.root");
  fileHMdataS0_05[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st5.root");
  fileHMdataST_05[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAL_ApL_st5.root");
  fileHMdataS0_05[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st5.root");
  fileHMdataST_05[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_LAL_ALL_st5.root");

  std::cout << "Reading SpherIcity/SpherOcity [0.9,1] \n";
  st_label_06 = "0.9 < S_{T} < 1";
  s0_label_06 = "0.9 < S_{0} < 1";
  fileHMdataS0_06[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st6.root");
  fileHMdataST_06[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAp_App_st6.root");
  fileHMdataS0_06[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st6.root");
  fileHMdataST_06[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAL_ApL_st6.root");
  fileHMdataS0_06[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st6.root");
  fileHMdataST_06[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_LAL_ALL_st6.root");

  std::cout << "Reading NO SpherIcity/SpherOcity Cuts\n";
  st_label_08 = "No S_{T} cuts";
  s0_label_08 = "No S_{0} cuts";
  fileHMdataS0_08[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAp_App_st8.root");
  fileHMdataST_08[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAp_App_st8.root");
  fileHMdataS0_08[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_pAL_ApL_st8.root");
  fileHMdataST_08[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_pAL_ApL_st8.root");
  fileHMdataS0_08[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/S0_16082019/CFOutput_LAL_ALL_st8.root");
  fileHMdataST_08[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/ST_16082019/CFOutput_LAL_ALL_st8.root");

//REMEMBER THAT WHEN ACCESSING DIRECTLY FROM FILE USE FindObjectAny!!!!
   TH1F* hCFHMdataST_01[3];
   TH1F* hCFHMdataST_02[3];
   TH1F* hCFHMdataST_03[3];
   TH1F* hCFHMdataST_04[3];
   TH1F* hCFHMdataST_05[3];
   TH1F* hCFHMdataST_06[3];
   TH1F* hCFHMdataST_08[3];

   TH1F* hCFHMdataS0_01[3];
   TH1F* hCFHMdataS0_02[3];
   TH1F* hCFHMdataS0_03[3];
   TH1F* hCFHMdataS0_04[3];
   TH1F* hCFHMdataS0_05[3];
   TH1F* hCFHMdataS0_06[3];
   TH1F* hCFHMdataS0_08[3];
//pAp flag 1 = bin 4 MeV
//pAL/LAL flag 1 = bin 16 MeV
   for(int i=0;i<3;i++){
     hCFHMdataST_01[i] = (TH1F*)(fileHMdataST_01[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_01[i]->Sumw2();
     hCFHMdataST_02[i] = (TH1F*)(fileHMdataST_02[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_02[i]->Sumw2();
     hCFHMdataST_03[i] = (TH1F*)(fileHMdataST_03[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_03[i]->Sumw2();
     hCFHMdataST_04[i] = (TH1F*)(fileHMdataST_04[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_04[i]->Sumw2();
     hCFHMdataST_05[i] = (TH1F*)(fileHMdataST_05[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_05[i]->Sumw2();
     hCFHMdataST_06[i] = (TH1F*)(fileHMdataST_06[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_06[i]->Sumw2();
     hCFHMdataST_08[i] = (TH1F*)(fileHMdataST_08[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataST_08[i]->Sumw2();
     DreamPlot::SetStyleHistoCF(hCFHMdataST_01[i], 8, kMagenta+2, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_02[i], 8, kMagenta-4, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_03[i], 8, kGreen+2, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_04[i], 8, kRed+1, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_05[i], 8, kOrange+7, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_06[i], 8, kAzure-2, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataST_08[i], 8, kCyan+1, 20);

     hCFHMdataS0_01[i] = (TH1F*)(fileHMdataS0_01[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_01[i]->Sumw2();
     hCFHMdataS0_02[i] = (TH1F*)(fileHMdataS0_02[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_02[i]->Sumw2();
     hCFHMdataS0_03[i] = (TH1F*)(fileHMdataS0_03[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_03[i]->Sumw2();
     hCFHMdataS0_04[i] = (TH1F*)(fileHMdataS0_04[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_04[i]->Sumw2();
     hCFHMdataS0_05[i] = (TH1F*)(fileHMdataS0_05[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_05[i]->Sumw2();
     hCFHMdataS0_06[i] = (TH1F*)(fileHMdataS0_06[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_06[i]->Sumw2();
     hCFHMdataS0_08[i] = (TH1F*)(fileHMdataS0_08[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdataS0_08[i]->Sumw2();
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_01[i], 22, kMagenta-2, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_02[i], 22, kViolet+1, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_03[i], 22, kGreen+4, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_04[i], 22, kRed+3, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_05[i], 22, kOrange-2, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_06[i], 22, kAzure+6, 20);
     DreamPlot::SetStyleHistoCF(hCFHMdataS0_08[i], 22, kBlue, 20);

   }


   fOut->cd();
   for(int i=0;i<3;i++)
   {
     hCFHMdataST_01[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add1.Data()),TObject::kOverwrite);
     hCFHMdataST_02[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add2.Data()),TObject::kOverwrite);
     hCFHMdataST_03[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add3.Data()),TObject::kOverwrite);
     hCFHMdataST_04[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add4.Data()),TObject::kOverwrite);
     hCFHMdataST_05[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add5.Data()),TObject::kOverwrite);
     hCFHMdataST_06[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add6.Data()),TObject::kOverwrite);
     hCFHMdataST_08[i]->Write(TString::Format("hCFHMdataST_%i_%s",i,add8.Data()),TObject::kOverwrite);

     hCFHMdataS0_01[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add1.Data()),TObject::kOverwrite);
     hCFHMdataS0_02[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add2.Data()),TObject::kOverwrite);
     hCFHMdataS0_03[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add3.Data()),TObject::kOverwrite);
     hCFHMdataS0_04[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add4.Data()),TObject::kOverwrite);
     hCFHMdataS0_05[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add5.Data()),TObject::kOverwrite);
     hCFHMdataS0_06[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add6.Data()),TObject::kOverwrite);
     hCFHMdataS0_08[i]->Write(TString::Format("hCFHMdataS0_%i_%s",i,add8.Data()),TObject::kOverwrite);

   }
   fOut->Close();
TCanvas *can_3GeV_ST[3]; TLegend *leg_3GeV_ST[3];
TCanvas *can_1GeV_ST[3]; TLegend *leg_1GeV_ST[3];
TCanvas *can_800MeV_ST[3];
TCanvas *can_500MeV_ST[3]; TLegend *leg_500MeV_ST[3];
TCanvas *can_250MeV_ST[3]; TLegend *leg_250MeV_ST[3];

TCanvas *can_3GeV_S0[3]; TLegend *leg_3GeV_S0[3];
TCanvas *can_1GeV_S0[3]; TLegend *leg_1GeV_S0[3];
TCanvas *can_800MeV_S0[3];
TCanvas *can_500MeV_S0[3]; TLegend *leg_500MeV_S0[3];
TCanvas *can_250MeV_S0[3]; TLegend *leg_250MeV_S0[3];

TCanvas *can_3GeV_03[3]; TLegend *leg_3GeV_03[3];
TCanvas *can_1GeV_03[3]; TLegend *leg_1GeV_03[3];
TCanvas *can_800MeV_03[3];
TCanvas *can_500MeV_03[3]; TLegend *leg_500MeV_03[3];
TCanvas *can_250MeV_03[3]; TLegend *leg_250MeV_03[3];

TCanvas *can_3GeV_05[3]; TLegend *leg_3GeV_05[3];
TCanvas *can_1GeV_05[3]; TLegend *leg_1GeV_05[3];
TCanvas *can_800MeV_05[3];
TCanvas *can_500MeV_05[3]; TLegend *leg_500MeV_05[3];
TCanvas *can_250MeV_05[3]; TLegend *leg_250MeV_05[3];

TCanvas *can_3GeV_06[3]; TLegend *leg_3GeV_06[3];
TCanvas *can_1GeV_06[3]; TLegend *leg_1GeV_06[3];
TCanvas *can_800MeV_06[3];
TCanvas *can_500MeV_06[3]; TLegend *leg_500MeV_06[3];
TCanvas *can_250MeV_06[3]; TLegend *leg_250MeV_06[3];

// Sphericity plots

    //----------1GeV-------------------------
    for(int i=0;i<3;i++){

      can_1GeV_ST[i] = new TCanvas("can_1GeV_ST[i]", "can_1GeV_ST[i]");
      hCFHMdataST_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
      hCFHMdataST_03[i] ->GetYaxis()->SetTitle("C(k*)");
      hCFHMdataST_03[i] ->GetXaxis()->SetRangeUser(0.,1000);
      hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
      hCFHMdataST_03[i] ->Draw("ep");
      hCFHMdataST_05[i] ->Draw("same ep");
      hCFHMdataST_06[i] ->Draw("same ep");
      hCFHMdataST_08[i] ->Draw("same ep");
       leg_1GeV_ST[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
       leg_1GeV_ST[i]->SetFillStyle(0);
       leg_1GeV_ST[i]->SetLineColor(0);
       // gStyle->SetLegendTextSize(0.05);
       leg_1GeV_ST[i]->AddEntry((TObject*)0, "", "");
       if(i==0){
         leg_1GeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
       }else if(i==1){
         leg_1GeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
       }else if(i==2){
         leg_1GeV_ST[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
       }
       leg_1GeV_ST[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_ST[i]->AddEntry(hCFHMdataST_03[i],st_label_03);
       leg_1GeV_ST[i]->AddEntry(hCFHMdataST_05[i],st_label_05);
       leg_1GeV_ST[i]->AddEntry(hCFHMdataST_06[i],st_label_06);
       leg_1GeV_ST[i]->AddEntry(hCFHMdataST_08[i],st_label_08);
       leg_1GeV_ST[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_ST[i]->Draw();
       can_1GeV_ST[i]->SaveAs(OutputFolder+TString::Format("CFSphericity_1GeV_sys%i.pdf",i));
// 500 MeV
       can_500MeV_ST[i] = new TCanvas("can_500MeV_ST[i]", "can_500MeV_ST[i]");
       hCFHMdataST_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
       hCFHMdataST_03[i] ->GetYaxis()->SetTitle("C(k*)");
       hCFHMdataST_03[i] ->GetXaxis()->SetRangeUser(0.,500);
       if(i==0){
       hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
     }else if(i==1){
       hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }else if(i==2){
       hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }
       hCFHMdataST_03[i] ->Draw("ep");
       hCFHMdataST_05[i] ->Draw("same ep");
       hCFHMdataST_06[i] ->Draw("same ep");
       hCFHMdataST_08[i] ->Draw("same ep");
       if(i==0){
        leg_500MeV_ST[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_500MeV_ST[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_500MeV_ST[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_500MeV_ST[i]->SetFillStyle(0);
        leg_500MeV_ST[i]->SetLineColor(0);
        // gStyle->SetLegendTextSize(0.05);
        leg_500MeV_ST[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_500MeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_500MeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_500MeV_ST[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_500MeV_ST[i]->AddEntry((TObject*)0, "", "");
        leg_500MeV_ST[i]->AddEntry(hCFHMdataST_03[i],st_label_03);
        leg_500MeV_ST[i]->AddEntry(hCFHMdataST_05[i],st_label_05);
        leg_500MeV_ST[i]->AddEntry(hCFHMdataST_06[i],st_label_06);
        leg_500MeV_ST[i]->AddEntry(hCFHMdataST_08[i],st_label_08);
        leg_500MeV_ST[i]->AddEntry((TObject*)0, "", "");
        leg_500MeV_ST[i]->Draw();
        can_500MeV_ST[i]->SaveAs(OutputFolder+TString::Format("CFSphericity_500MeV_sys%i.pdf",i));
// 250 MeV
        can_250MeV_ST[i] = new TCanvas("can_250MeV_ST[i]", "can_250MeV_ST[i]");
        hCFHMdataST_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
        hCFHMdataST_03[i] ->GetYaxis()->SetTitle("C(k*)");
        hCFHMdataST_03[i] ->GetXaxis()->SetRangeUser(0.,250);
        if(i==0){
        hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
      }else if(i==1){
        hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
      }else if(i==2){
        hCFHMdataST_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
      }
        hCFHMdataST_03[i] ->Draw("ep");
       hCFHMdataST_05[i] ->Draw("same ep");
       hCFHMdataST_06[i] ->Draw("same ep");
       hCFHMdataST_08[i] ->Draw("same ep");
       if(i==0){
        leg_250MeV_ST[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_250MeV_ST[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_250MeV_ST[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_250MeV_ST[i]->SetFillStyle(0);
        leg_250MeV_ST[i]->SetLineColor(0);
        leg_250MeV_ST[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_250MeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_250MeV_ST[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_250MeV_ST[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_250MeV_ST[i]->AddEntry((TObject*)0, "", "");
        leg_250MeV_ST[i]->AddEntry(hCFHMdataST_03[i],st_label_03);
        leg_250MeV_ST[i]->AddEntry(hCFHMdataST_05[i],st_label_05);
        leg_250MeV_ST[i]->AddEntry(hCFHMdataST_06[i],st_label_06);
        leg_250MeV_ST[i]->AddEntry(hCFHMdataST_08[i],st_label_08);
        leg_250MeV_ST[i]->AddEntry((TObject*)0, "", "");
        leg_250MeV_ST[i]->Draw();
        can_250MeV_ST[i]->SaveAs(OutputFolder+TString::Format("CFSphericity_250MeV_sys%i.pdf",i));
    }

// Spherocity plots
    //----------1GeV-------------------------
    for(int i=0;i<3;i++){
      can_1GeV_S0[i] = new TCanvas("can_1GeV_S0[i]", "can_1GeV_S0[i]");
      hCFHMdataS0_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
      hCFHMdataS0_03[i] ->GetYaxis()->SetTitle("C(k*)");
      hCFHMdataS0_03[i] ->GetXaxis()->SetRangeUser(0.,1000);
      hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
      hCFHMdataS0_03[i] ->Draw("ep");
      hCFHMdataS0_05[i] ->Draw("same ep");
      hCFHMdataS0_06[i] ->Draw("same ep");
      hCFHMdataS0_08[i] ->Draw("same ep");
       leg_1GeV_S0[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
       leg_1GeV_S0[i]->SetFillStyle(0);
       leg_1GeV_S0[i]->SetLineColor(0);
       // gStyle->SetLegendTextSize(0.05);
       leg_1GeV_S0[i]->AddEntry((TObject*)0, "", "");
       if(i==0){
         leg_1GeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
       }else if(i==1){
         leg_1GeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
       }else if(i==2){
         leg_1GeV_S0[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
       }
       leg_1GeV_S0[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_S0[i]->AddEntry(hCFHMdataS0_03[i],s0_label_03);
       leg_1GeV_S0[i]->AddEntry(hCFHMdataS0_05[i],s0_label_05);
       leg_1GeV_S0[i]->AddEntry(hCFHMdataS0_06[i],s0_label_06);
       leg_1GeV_S0[i]->AddEntry(hCFHMdataS0_08[i],s0_label_08);
       leg_1GeV_S0[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_S0[i]->Draw();
       can_1GeV_S0[i]->SaveAs(OutputFolder+TString::Format("CFSpherocity_1GeV_sys%i.pdf",i));
// 500 MeV
       can_500MeV_S0[i] = new TCanvas("can_500MeV_S0[i]", "can_500MeV_S0[i]");
       hCFHMdataS0_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
       hCFHMdataS0_03[i] ->GetYaxis()->SetTitle("C(k*)");
       hCFHMdataS0_03[i] ->GetXaxis()->SetRangeUser(0.,500);
       if(i==0){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
     }else if(i==1){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }else if(i==2){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }
       hCFHMdataS0_03[i] ->Draw("ep");
       hCFHMdataS0_05[i] ->Draw("same ep");
       hCFHMdataS0_06[i] ->Draw("same ep");
       hCFHMdataS0_08[i] ->Draw("same ep");
       if(i==0){
        leg_500MeV_S0[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_500MeV_S0[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_500MeV_S0[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_500MeV_S0[i]->SetFillStyle(0);
        leg_500MeV_S0[i]->SetLineColor(0);
        // gStyle->SetLegendTextSize(0.05);
        leg_500MeV_S0[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_500MeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_500MeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_500MeV_S0[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_500MeV_S0[i]->AddEntry((TObject*)0, "", "");
        leg_500MeV_S0[i]->AddEntry(hCFHMdataS0_03[i],s0_label_03);
        leg_500MeV_S0[i]->AddEntry(hCFHMdataS0_05[i],s0_label_05);
        leg_500MeV_S0[i]->AddEntry(hCFHMdataS0_06[i],s0_label_06);
        leg_500MeV_S0[i]->AddEntry(hCFHMdataS0_08[i],s0_label_08);
        leg_500MeV_S0[i]->AddEntry((TObject*)0, "", "");
        leg_500MeV_S0[i]->Draw();
        can_500MeV_S0[i]->SaveAs(OutputFolder+TString::Format("CFSpherocity_500MeV_sys%i.pdf",i));
// 250 MeV
       can_250MeV_S0[i] = new TCanvas("can_250MeV_S0[i]", "can_250MeV_S0[i]");
       hCFHMdataS0_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
       hCFHMdataS0_03[i] ->GetYaxis()->SetTitle("C(k*)");
       hCFHMdataS0_03[i] ->GetXaxis()->SetRangeUser(0.,250);
       if(i==0){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
     }else if(i==1){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }else if(i==2){
       hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }
       hCFHMdataS0_03[i] ->Draw("ep");
       hCFHMdataS0_05[i] ->Draw("same ep");
       hCFHMdataS0_06[i] ->Draw("same ep");
       hCFHMdataS0_08[i] ->Draw("same ep");
       if(i==0){
        leg_250MeV_S0[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_250MeV_S0[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_250MeV_S0[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_250MeV_S0[i]->SetFillStyle(0);
        leg_250MeV_S0[i]->SetLineColor(0);
        // gStyle->SetLegendTextSize(0.05);
        leg_250MeV_S0[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_250MeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_250MeV_S0[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_250MeV_S0[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_250MeV_S0[i]->AddEntry((TObject*)0, "", "");
        leg_250MeV_S0[i]->AddEntry(hCFHMdataS0_03[i],s0_label_03);
        leg_250MeV_S0[i]->AddEntry(hCFHMdataS0_05[i],s0_label_05);
        leg_250MeV_S0[i]->AddEntry(hCFHMdataS0_06[i],s0_label_06);
        leg_250MeV_S0[i]->AddEntry(hCFHMdataS0_08[i],s0_label_08);
        leg_250MeV_S0[i]->AddEntry((TObject*)0, "", "");
        leg_250MeV_S0[i]->Draw();
        can_250MeV_S0[i]->SaveAs(OutputFolder+TString::Format("CFSpherocity_250MeV_sys%i.pdf",i));
    }

// Spherocity vs Sphericity plots
    //----------1GeV-------------------------
    for(int i=0;i<3;i++){
      can_1GeV_03[i] = new TCanvas("can_1GeV_03[i]", "can_1GeV_03[i]");
      hCFHMdataS0_03[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
      hCFHMdataS0_03[i] ->GetYaxis()->SetTitle("C(k*)");
      hCFHMdataS0_03[i] ->GetXaxis()->SetRangeUser(0.,1000);
      hCFHMdataST_03[i] ->GetXaxis()->SetRangeUser(0.,1000);
      if(i==0){
      hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
    }else if(i==1){
      hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
    }else if(i==2){
      hCFHMdataS0_03[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
    }
      hCFHMdataS0_03[i] ->Draw("ep");
      hCFHMdataST_03[i] ->Draw("same ep");
      if(i==0){
       leg_1GeV_03[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
     }else if(i==1){
       leg_1GeV_03[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
     }else if(i==2){
       leg_1GeV_03[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
     }
       leg_1GeV_03[i]->SetFillStyle(0);
       leg_1GeV_03[i]->SetLineColor(0);
       gStyle->SetLegendTextSize(0.03);
       leg_1GeV_03[i]->AddEntry((TObject*)0, "", "");
       if(i==0){
         leg_1GeV_03[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
       }else if(i==1){
         leg_1GeV_03[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
       }else if(i==2){
         leg_1GeV_03[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
       }
       leg_1GeV_03[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_03[i]->AddEntry(hCFHMdataS0_03[i],s0_label_03);
       leg_1GeV_03[i]->AddEntry(hCFHMdataST_03[i],st_label_03);
       leg_1GeV_03[i]->AddEntry((TObject*)0, "", "");
       leg_1GeV_03[i]->Draw();
       can_1GeV_03[i]->SaveAs(OutputFolder+TString::Format("CFS0ST03_1GeV_sys%i_.pdf",i));


       can_1GeV_05[i] = new TCanvas("can_1GeV_05[i]", "can_1GeV_05[i]");
       hCFHMdataS0_05[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
       hCFHMdataS0_05[i] ->GetYaxis()->SetTitle("C(k*)");
       hCFHMdataS0_05[i] ->GetXaxis()->SetRangeUser(0.,1000);
       hCFHMdataST_05[i] ->GetXaxis()->SetRangeUser(0.,1000);
       if(i==0){
       hCFHMdataS0_05[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
     }else if(i==1){
       hCFHMdataS0_05[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }else if(i==2){
       hCFHMdataS0_05[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }
       hCFHMdataS0_05[i] ->Draw("ep");
       hCFHMdataST_05[i] ->Draw("same ep");
       if(i==0){
        leg_1GeV_05[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_1GeV_05[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_1GeV_05[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_1GeV_05[i]->SetFillStyle(0);
        leg_1GeV_05[i]->SetLineColor(0);
        gStyle->SetLegendTextSize(0.02);
        leg_1GeV_05[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_1GeV_05[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_1GeV_05[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_1GeV_05[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_1GeV_05[i]->AddEntry((TObject*)0, "", "");
        leg_1GeV_05[i]->AddEntry(hCFHMdataS0_05[i],s0_label_05);
        leg_1GeV_05[i]->AddEntry(hCFHMdataST_05[i],st_label_05);
        leg_1GeV_05[i]->AddEntry((TObject*)0, "", "");
        leg_1GeV_05[i]->Draw();
        can_1GeV_05[i]->SaveAs(OutputFolder+TString::Format("CFS0ST05_1GeV_sys%i_.pdf",i));


       can_1GeV_06[i] = new TCanvas("can_1GeV_06[i]", "can_1GeV_06[i]");
       hCFHMdataS0_06[i] ->GetXaxis()->SetTitle("k* [MeV/c]");
       hCFHMdataS0_06[i] ->GetYaxis()->SetTitle("C(k*)");
       hCFHMdataS0_06[i] ->GetXaxis()->SetRangeUser(0.,1000);
       hCFHMdataST_06[i] ->GetXaxis()->SetRangeUser(0.,1000);
       if(i==0){
       hCFHMdataS0_06[i] ->GetYaxis()->SetRangeUser(0.3,2.5);
     }else if(i==1){
       hCFHMdataS0_06[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }else if(i==2){
       hCFHMdataS0_06[i] ->GetYaxis()->SetRangeUser(0.3,1.2);
     }
       hCFHMdataS0_06[i] ->Draw("ep");
       hCFHMdataST_06[i] ->Draw("same ep");
       if(i==0){
        leg_1GeV_06[i] = new TLegend(0.6, 0.4, 0.9, 0.9);
      }else if(i==1){
        leg_1GeV_06[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }else if(i==2){
        leg_1GeV_06[i] = new TLegend(0.6, 0.2, 0.9, 0.7);
      }
        leg_1GeV_06[i]->SetFillStyle(0);
        leg_1GeV_06[i]->SetLineColor(0);
        gStyle->SetLegendTextSize(0.02);
        leg_1GeV_06[i]->AddEntry((TObject*)0, "", "");
        if(i==0){
          leg_1GeV_06[i]->AddEntry((TObject*)0, "HM p-#bar{p}", "");
        }else if(i==1){
          leg_1GeV_06[i]->AddEntry((TObject*)0, "HM p-#bar{#Lambda}", "");
        }else if(i==2){
          leg_1GeV_06[i]->AddEntry((TObject*)0, "HM #Lambda-#bar{#Lambda}", "");
        }
        leg_1GeV_06[i]->AddEntry((TObject*)0, "", "");
        leg_1GeV_06[i]->AddEntry(hCFHMdataS0_06[i],s0_label_06);
        leg_1GeV_06[i]->AddEntry(hCFHMdataST_06[i],st_label_06);
        leg_1GeV_06[i]->AddEntry((TObject*)0, "", "");
        leg_1GeV_06[i]->Draw();
        can_1GeV_06[i]->SaveAs(OutputFolder+TString::Format("CFS0ST06_1GeV_sys%i_.pdf",i));
    }

return 0;
}
