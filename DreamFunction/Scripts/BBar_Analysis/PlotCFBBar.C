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
  const char* STflag = argv[1];
  const char* date = (argv[2]) ? argv[2] : "";
  const char* maxkstar = (argv[3]) ? argv[3] : "";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetLegendTextSize(0.5);


  DreamPlot::SetStyle();

  TString add1="1";
  TString add2="2";
  TString add3="3";
  TString add4="4";
  TString add5="5";
  TString add6="6";
  TString add8="8";
  TString addon;
  TString st_label;



  double charmask = atof (maxkstar);

  std::cout<<"charmask = "<<charmask<<std::endl;

  gSystem->Exec(TString::Format("mkdir ../../../BBbar/GentleFemto_Output/CFPlots/%s",date));
  TString OutputFolder = TString::Format("../../../BBbar/GentleFemto_Output/CFPlots/%s/",date);
  TFile* fOut = new TFile(TString::Format(OutputFolder+"ComparisonHMMB.root",date),"update");

//Accessing all the directory in the root file
  TFile* fileMBdata[4];
  TFile* fileHMdata[4];


  if(strcmp(STflag, add1)==0)
  {
    addon=add1;
    st_label = "0 < s_{T} < 0.3";
  std::cout << "Sphericity [0.,0.3] \n";
  fileMBdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAp_App_st1.root");
  fileHMdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st1.root");
  fileMBdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAL_ApL_st1.root");
  fileHMdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st1.root");
  fileMBdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_LAL_ALL_st1.root");
  fileHMdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st1.root");
  }
  else if(strcmp(STflag, add2)==0)
  {
    addon=add2;
    st_label = "0.3 < s_{T} < 0.7";
  std::cout << "Sphericity [0.3,0.7] \n";
  fileMBdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAp_App_st2.root");
  fileHMdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st2.root");
  fileMBdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAL_ApL_st2.root");
  fileHMdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st2.root");
  fileMBdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_LAL_ALL_st2.root");
  fileHMdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st2.root");
  }
  else if(strcmp(STflag, add3)==0)
  {
    addon=add3;
    st_label = "0.7 < s_{T} < 1";
  std::cout << "Sphericity [0.7,1.] \n";
  fileMBdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAp_App_st3.root");
  fileHMdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_st3.root");
  fileMBdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAL_ApL_st3.root");
  fileHMdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_st3.root");
  fileMBdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_LAL_ALL_st3.root");
  fileHMdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_st3.root");
  }
  else if(strcmp(STflag, add4)==0)
  {
    addon=add4;
    st_label = "0 < s_{T} < 1";
    fileMBdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAp_App_full.root");
    fileHMdata[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAp_App_full.root");
    fileMBdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_pAL_ApL_full.root");
    fileHMdata[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_pAL_ApL_full.root");
    fileMBdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/MB/CFOutput_LAL_ALL_full.root");
    fileHMdata[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/HM/22July2019/CFOutput_LAL_ALL_full.root");
  }

//REMEMBER THAT WHEN ACCESSING DIRECTLY FROM FILE USE FindObjectAny!!!!
   TH1F* hCFMBdata[3];
   TH1F* hCFHMdata[3];

   for(int i=0;i<3;i++){
     hCFMBdata[i] = (TH1F*)(fileMBdata[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFMBdata[i]->Sumw2();
     hCFHMdata[i] = (TH1F*)(fileHMdata[i]->FindObjectAny("hCk_ReweightedMeV_1"));
     hCFHMdata[i]->Sumw2();
   }

   DreamPlot::SetStyleHistoCF(hCFMBdata[0], 8, kAzure+8, 20);
   DreamPlot::SetStyleHistoCF(hCFMBdata[1], 8, kOrange-3, 20);
   DreamPlot::SetStyleHistoCF(hCFMBdata[2], 8, kGreen+4, 20);

   DreamPlot::SetStyleHistoCF(hCFHMdata[0], 8, kBlue, 20);
   DreamPlot::SetStyleHistoCF(hCFHMdata[1], 8, kRed, 20);
   DreamPlot::SetStyleHistoCF(hCFHMdata[2], 8, kGreen, 20);

   fOut->cd();
   for(int i=0;i<3;i++)
   {
     hCFMBdata[i]->Write(TString::Format("hCFMBdata_%i_%s",i,addon.Data()),TObject::kOverwrite);
     hCFHMdata[i]->Write(TString::Format("hCFHMdata_%i_%s",i,addon.Data()),TObject::kOverwrite);
   }
   fOut->Close();

   TCanvas *can_MB_3GeV[3];   TCanvas *can_MB_1GeV[3];   TCanvas *can_MB_500MeV[3];   TCanvas *can_MB_300MeV[3];
   TLegend *leg_MB_3GeV[3];   TLegend *leg_MB_1GeV[3];   TLegend *leg_MB_500MeV[3];   TLegend *leg_MB_300MeV[3];
   TCanvas *can_HM_3GeV[3];   TCanvas *can_HM_1GeV[3];   TCanvas *can_HM_500MeV[3];   TCanvas *can_HM_300MeV[3];
   TLegend *leg_HM_3GeV[3];   TLegend *leg_HM_1GeV[3];   TLegend *leg_HM_500MeV[3];   TLegend *leg_HM_300MeV[3];
   //----------3GeV-------------------------
   can_MB_3GeV[0] = new TCanvas("can_MB_3GeV[0]", "can_MB_3GeV[0]");
   hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
   hCFMBdata[0] ->Draw();
   leg_MB_3GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_MB_3GeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
   leg_MB_3GeV[0]->AddEntry((TObject*)0, "", "");
   leg_MB_3GeV[0]->AddEntry((TObject*)0, st_label, "");
   leg_MB_3GeV[0]->Draw();
   can_MB_3GeV[0]->SaveAs(OutputFolder+TString::Format("MBCF_3GeV_pantip_%s.pdf",addon.Data()));
   can_MB_3GeV[1] = new TCanvas("can_MB_3GeV[1]", "can_MB_3GeV[1]");
   hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
   hCFMBdata[1] ->Draw();
   leg_MB_3GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_MB_3GeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
   leg_MB_3GeV[1]->AddEntry((TObject*)0, "", "");
   leg_MB_3GeV[1]->AddEntry((TObject*)0, st_label, "");
   leg_MB_3GeV[1]->Draw();
   can_MB_3GeV[1]->SaveAs(OutputFolder+TString::Format("MBCF_3GeV_pantiL_%s.pdf",addon.Data()));
   can_MB_3GeV[2] = new TCanvas("can_MB_3GeV[2]", "can_MB_3GeV[2]");
   hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
   hCFMBdata[2] ->Draw();
   leg_MB_3GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_MB_3GeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
   leg_MB_3GeV[2]->AddEntry((TObject*)0, "", "");
   leg_MB_3GeV[2]->AddEntry((TObject*)0, st_label, "");
   leg_MB_3GeV[2]->Draw();
   can_MB_3GeV[2]->SaveAs(OutputFolder+TString::Format("MBCF_3GeV_LantiL_%s.pdf",addon.Data()));

   can_HM_3GeV[0] = new TCanvas("can_HM_3GeV[0]", "can_HM_3GeV[0]");
   hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
   hCFHMdata[0] ->Draw();
   leg_HM_3GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_HM_3GeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
   leg_HM_3GeV[0]->AddEntry((TObject*)0, "", "");
   leg_HM_3GeV[0]->AddEntry((TObject*)0, st_label, "");
   leg_HM_3GeV[0]->Draw();
   can_HM_3GeV[0]->SaveAs(OutputFolder+TString::Format("HMCF_3GeV_pantip_%s.pdf",addon.Data()));
   can_HM_3GeV[1] = new TCanvas("can_HM_3GeV[1]", "can_HM_3GeV[1]");
   hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFHMdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
   hCFHMdata[1] ->Draw();
   leg_HM_3GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_HM_3GeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
   leg_HM_3GeV[1]->AddEntry((TObject*)0, "", "");
   leg_HM_3GeV[1]->AddEntry((TObject*)0, st_label, "");
   leg_HM_3GeV[1]->Draw();
   can_HM_3GeV[1]->SaveAs(OutputFolder+TString::Format("HMCF_3GeV_pantiL_%s.pdf",addon.Data()));
   can_HM_3GeV[2] = new TCanvas("can_HM_3GeV[2]", "can_HM_3GeV[2]");
   hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
   hCFHMdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
   hCFHMdata[2] ->Draw();
   leg_HM_3GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
   leg_HM_3GeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
   leg_HM_3GeV[2]->AddEntry((TObject*)0, "", "");
   leg_HM_3GeV[2]->AddEntry((TObject*)0, st_label, "");
   leg_HM_3GeV[2]->Draw();
   can_HM_3GeV[2]->SaveAs(OutputFolder+TString::Format("HMCF_3GeV_LantiL_%s.pdf",addon.Data()));
//----------1GeV-------------------------
can_MB_1GeV[0] = new TCanvas("can_MB_1GeV[0]", "can_MB_1GeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
leg_MB_1GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_1GeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MB_1GeV[0]->AddEntry((TObject*)0, "", "");
leg_MB_1GeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MB_1GeV[0]->Draw();
can_MB_1GeV[0]->SaveAs(OutputFolder+TString::Format("MBCF_1GeV_pantip_%s.pdf",addon.Data()));
can_MB_1GeV[1] = new TCanvas("can_MB_1GeV[1]", "can_MB_1GeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[1] ->Draw();
leg_MB_1GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_1GeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MB_1GeV[1]->AddEntry((TObject*)0, "", "");
leg_MB_1GeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MB_1GeV[1]->Draw();
can_MB_1GeV[1]->SaveAs(OutputFolder+TString::Format("MBCF_1GeV_pantiL_%s.pdf",addon.Data()));
can_MB_1GeV[2] = new TCanvas("can_MB_1GeV[2]", "can_MB_1GeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[2] ->Draw();
leg_MB_1GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_1GeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MB_1GeV[2]->AddEntry((TObject*)0, "", "");
leg_MB_1GeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MB_1GeV[2]->Draw();
can_MB_1GeV[2]->SaveAs(OutputFolder+TString::Format("MBCF_1GeV_LantiL_%s.pdf",addon.Data()));

can_HM_1GeV[0] = new TCanvas("can_HM_1GeV[0]", "can_HM_1GeV[0]");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFHMdata[0] ->Draw();
leg_HM_1GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_1GeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_HM_1GeV[0]->AddEntry((TObject*)0, "", "");
leg_HM_1GeV[0]->AddEntry((TObject*)0, st_label, "");
leg_HM_1GeV[0]->Draw();
can_HM_1GeV[0]->SaveAs(OutputFolder+TString::Format("HMCF_1GeV_pantip_%s.pdf",addon.Data()));
can_HM_1GeV[1] = new TCanvas("can_HM_1GeV[1]", "can_HM_1GeV[1]");
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[1] ->Draw();
leg_HM_1GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_1GeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_HM_1GeV[1]->AddEntry((TObject*)0, "", "");
leg_HM_1GeV[1]->AddEntry((TObject*)0, st_label, "");
leg_HM_1GeV[1]->Draw();
can_HM_1GeV[1]->SaveAs(OutputFolder+TString::Format("HMCF_1GeV_pantiL_%s.pdf",addon.Data()));
can_HM_1GeV[2] = new TCanvas("can_HM_1GeV[2]", "can_HM_1GeV[2]");
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[2] ->Draw();
leg_HM_1GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_1GeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_HM_1GeV[2]->AddEntry((TObject*)0, "", "");
leg_HM_1GeV[2]->AddEntry((TObject*)0, st_label, "");
leg_HM_1GeV[2]->Draw();
can_HM_1GeV[2]->SaveAs(OutputFolder+TString::Format("HMCF_1GeV_LantiL_%s.pdf",addon.Data()));
//----------500MeV-------------------------
can_MB_500MeV[0] = new TCanvas("can_MB_500MeV[0]", "can_MB_500MeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
leg_MB_500MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_500MeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MB_500MeV[0]->AddEntry((TObject*)0, "", "");
leg_MB_500MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MB_500MeV[0]->Draw();
can_MB_500MeV[0]->SaveAs(OutputFolder+TString::Format("MBCF_500MeV_pantip_%s.pdf",addon.Data()));
can_MB_500MeV[1] = new TCanvas("can_MB_500MeV[1]", "can_MB_500MeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[1] ->Draw();
leg_MB_500MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_500MeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MB_500MeV[1]->AddEntry((TObject*)0, "", "");
leg_MB_500MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MB_500MeV[1]->Draw();
can_MB_500MeV[1]->SaveAs(OutputFolder+TString::Format("MBCF_500MeV_pantiL_%s.pdf",addon.Data()));
can_MB_500MeV[2] = new TCanvas("can_MB_500MeV[2]", "can_MB_500MeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[2] ->Draw();
leg_MB_500MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_500MeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MB_500MeV[2]->AddEntry((TObject*)0, "", "");
leg_MB_500MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MB_500MeV[2]->Draw();
can_MB_500MeV[2]->SaveAs(OutputFolder+TString::Format("MBCF_500MeV_LantiL_%s.pdf",addon.Data()));

can_HM_500MeV[0] = new TCanvas("can_HM_500MeV[0]", "can_HM_500MeV[0]");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFHMdata[0] ->Draw();
leg_HM_500MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_500MeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_HM_500MeV[0]->AddEntry((TObject*)0, "", "");
leg_HM_500MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_HM_500MeV[0]->Draw();
can_HM_500MeV[0]->SaveAs(OutputFolder+TString::Format("HMCF_500MeV_pantip_%s.pdf",addon.Data()));
can_HM_500MeV[1] = new TCanvas("can_HM_500MeV[1]", "can_HM_500MeV[1]");
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[1] ->Draw();
leg_HM_500MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_500MeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_HM_500MeV[1]->AddEntry((TObject*)0, "", "");
leg_HM_500MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_HM_500MeV[1]->Draw();
can_HM_500MeV[1]->SaveAs(OutputFolder+TString::Format("HMCF_500MeV_pantiL_%s.pdf",addon.Data()));
can_HM_500MeV[2] = new TCanvas("can_HM_500MeV[2]", "can_HM_500MeV[2]");
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[2] ->Draw();
leg_HM_500MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_500MeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_HM_500MeV[2]->AddEntry((TObject*)0, "", "");
leg_HM_500MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_HM_500MeV[2]->Draw();
can_HM_500MeV[2]->SaveAs(OutputFolder+TString::Format("HMCF_500MeV_LantiL_%s.pdf",addon.Data()));
//----------300MeV-------------------------
can_MB_300MeV[0] = new TCanvas("can_MB_300MeV[0]", "can_MB_300MeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
leg_MB_300MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_300MeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MB_300MeV[0]->AddEntry((TObject*)0, "", "");
leg_MB_300MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MB_300MeV[0]->Draw();
can_MB_300MeV[0]->SaveAs(OutputFolder+TString::Format("MBCF_300MeV_pantip_%s.pdf",addon.Data()));
can_MB_300MeV[1] = new TCanvas("can_MB_300MeV[1]", "can_MB_300MeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[1] ->Draw();
leg_MB_300MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_300MeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MB_300MeV[1]->AddEntry((TObject*)0, "", "");
leg_MB_300MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MB_300MeV[1]->Draw();
can_MB_300MeV[1]->SaveAs(OutputFolder+TString::Format("MBCF_300MeV_pantiL_%s.pdf",addon.Data()));
can_MB_300MeV[2] = new TCanvas("can_MB_300MeV[2]", "can_MB_300MeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFMBdata[2] ->Draw();
leg_MB_300MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MB_300MeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MB_300MeV[2]->AddEntry((TObject*)0, "", "");
leg_MB_300MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MB_300MeV[2]->Draw();
can_MB_300MeV[2]->SaveAs(OutputFolder+TString::Format("MBCF_300MeV_LantiL_%s.pdf",addon.Data()));

can_HM_300MeV[0] = new TCanvas("can_HM_300MeV[0]", "can_HM_300MeV[0]");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFHMdata[0] ->Draw();
leg_HM_300MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_300MeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_HM_300MeV[0]->AddEntry((TObject*)0, "", "");
leg_HM_300MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_HM_300MeV[0]->Draw();
can_HM_300MeV[0]->SaveAs(OutputFolder+TString::Format("HMCF_300MeV_pantip_%s.pdf",addon.Data()));
can_HM_300MeV[1] = new TCanvas("can_HM_300MeV[1]", "can_HM_300MeV[1]");
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[1] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[1] ->Draw();
leg_HM_300MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_300MeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_HM_300MeV[1]->AddEntry((TObject*)0, "", "");
leg_HM_300MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_HM_300MeV[1]->Draw();
can_HM_300MeV[1]->SaveAs(OutputFolder+TString::Format("HMCF_300MeV_pantiL_%s.pdf",addon.Data()));
can_HM_300MeV[2] = new TCanvas("can_HM_300MeV[2]", "can_HM_300MeV[2]");
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[2] ->GetYaxis()->SetRangeUser(0.4,1.5);
hCFHMdata[2] ->Draw();
leg_HM_300MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HM_300MeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_HM_300MeV[2]->AddEntry((TObject*)0, "", "");
leg_HM_300MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_HM_300MeV[2]->Draw();
can_HM_300MeV[2]->SaveAs(OutputFolder+TString::Format("HMCF_300MeV_LantiL_%s.pdf",addon.Data()));

TCanvas *can_MBtot_3GeV;   TCanvas *can_MBtot_1GeV;   TCanvas *can_MBtot_500MeV;   TCanvas *can_MBtot_300MeV;
TLegend *leg_MBtot_3GeV;   TLegend *leg_MBtot_1GeV;   TLegend *leg_MBtot_500MeV;   TLegend *leg_MBtot_300MeV;
TCanvas *can_HMtot_3GeV;   TCanvas *can_HMtot_1GeV;   TCanvas *can_HMtot_500MeV;   TCanvas *can_HMtot_300MeV;
TLegend *leg_HMtot_3GeV;   TLegend *leg_HMtot_1GeV;   TLegend *leg_HMtot_500MeV;   TLegend *leg_HMtot_300MeV;

can_MBtot_3GeV = new TCanvas("can_MBtot_3GeV", "can_MBtot_3GeV");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFMBdata[0] ->Draw();
hCFMBdata[1] ->Draw("same");
hCFMBdata[2] ->Draw("same");
leg_MBtot_3GeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBtot_3GeV->AddEntry((TObject*)0, "MB", "");
leg_MBtot_3GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_3GeV->AddEntry(hCFMBdata[0],"p-#bar{p}");
leg_MBtot_3GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_3GeV->AddEntry(hCFMBdata[1],"p-#bar{#Lambda}");
leg_MBtot_3GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_3GeV->AddEntry(hCFMBdata[2],"#Lambda-#bar{#Lambda}");
leg_MBtot_3GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_3GeV->AddEntry((TObject*)0, st_label, "");
leg_MBtot_3GeV->Draw();
can_MBtot_3GeV->SaveAs(OutputFolder+TString::Format("MBCF_3GeV_tot_%s.pdf",addon.Data()));

can_MBtot_1GeV = new TCanvas("can_MBtot_1GeV", "can_MBtot_1GeV");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFMBdata[0] ->Draw();
hCFMBdata[1] ->Draw("same");
hCFMBdata[2] ->Draw("same");
leg_MBtot_1GeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBtot_1GeV->AddEntry((TObject*)0, "MB", "");
leg_MBtot_1GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_1GeV->AddEntry(hCFMBdata[0],"p-#bar{p}");
leg_MBtot_1GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_1GeV->AddEntry(hCFMBdata[1],"p-#bar{#Lambda}");
leg_MBtot_1GeV->AddEntry((TObject*)0, "", "");
leg_MBtot_1GeV->AddEntry(hCFMBdata[2],"#Lambda-#bar{#Lambda}");
leg_MBtot_1GeV->AddEntry((TObject*)0, st_label, "");
leg_MBtot_1GeV->Draw();
can_MBtot_1GeV->SaveAs(OutputFolder+TString::Format("MBCF_1GeV_tot_%s.pdf",addon.Data()));

can_MBtot_500MeV = new TCanvas("can_MBtot_500MeV", "can_MBtot_500MeV");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFMBdata[0] ->Draw();
hCFMBdata[1] ->Draw("same");
hCFMBdata[2] ->Draw("same");
leg_MBtot_500MeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBtot_500MeV->AddEntry((TObject*)0, "MB", "");
leg_MBtot_500MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_500MeV->AddEntry(hCFMBdata[0],"p-#bar{p}");
leg_MBtot_500MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_500MeV->AddEntry(hCFMBdata[1],"p-#bar{#Lambda}");
leg_MBtot_500MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_500MeV->AddEntry(hCFMBdata[2],"#Lambda-#bar{#Lambda}");
leg_MBtot_500MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_500MeV->AddEntry((TObject*)0, st_label, "");
leg_MBtot_500MeV->Draw();
can_MBtot_500MeV->SaveAs(OutputFolder+TString::Format("MBCF_500MeV_tot_%s.pdf",addon.Data()));

can_MBtot_300MeV = new TCanvas("can_MBtot_300MeV", "can_MBtot_300MeV");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFMBdata[0] ->Draw();
hCFMBdata[1] ->Draw("same");
hCFMBdata[2] ->Draw("same");
leg_MBtot_300MeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBtot_300MeV->AddEntry((TObject*)0, "MB", "");
leg_MBtot_300MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_300MeV->AddEntry(hCFMBdata[0],"p-#bar{p}");
leg_MBtot_300MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_300MeV->AddEntry(hCFMBdata[1],"p-#bar{#Lambda}");
leg_MBtot_300MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_300MeV->AddEntry(hCFMBdata[2],"#Lambda-#bar{#Lambda}");
leg_MBtot_300MeV->AddEntry((TObject*)0, "", "");
leg_MBtot_300MeV->AddEntry((TObject*)0, st_label, "");
leg_MBtot_300MeV->Draw();
can_MBtot_300MeV->SaveAs(OutputFolder+TString::Format("MBCF_300MeV_tot_%s.pdf",addon.Data()));


can_HMtot_3GeV = new TCanvas("can_HMtot_3GeV", "can_HMtot_3GeV");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFHMdata[0] ->Draw();
hCFHMdata[1] ->Draw("same");
hCFHMdata[2] ->Draw("same");
leg_HMtot_3GeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HMtot_3GeV->AddEntry((TObject*)0, "HM", "");
leg_HMtot_3GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_3GeV->AddEntry(hCFHMdata[0],"p-#bar{p}");
leg_HMtot_3GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_3GeV->AddEntry(hCFHMdata[1],"p-#bar{#Lambda}");
leg_HMtot_3GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_3GeV->AddEntry(hCFHMdata[2],"#Lambda-#bar{#Lambda}");
leg_HMtot_3GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_3GeV->AddEntry((TObject*)0, st_label, "");
leg_HMtot_3GeV->Draw();
can_HMtot_3GeV->SaveAs(OutputFolder+TString::Format("HMCF_3GeV_tot_%s.pdf",addon.Data()));

can_HMtot_1GeV = new TCanvas("can_HMtot_1GeV", "can_HMtot_1GeV");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFHMdata[0] ->Draw();
hCFHMdata[1] ->Draw("same");
hCFHMdata[2] ->Draw("same");
leg_HMtot_1GeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HMtot_1GeV->AddEntry((TObject*)0, "HM", "");
leg_HMtot_1GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_1GeV->AddEntry(hCFHMdata[0],"p-#bar{p}");
leg_HMtot_1GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_1GeV->AddEntry(hCFHMdata[1],"p-#bar{#Lambda}");
leg_HMtot_1GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_1GeV->AddEntry(hCFHMdata[2],"#Lambda-#bar{#Lambda}");
leg_HMtot_1GeV->AddEntry((TObject*)0, "", "");
leg_HMtot_1GeV->AddEntry((TObject*)0, st_label, "");
leg_HMtot_1GeV->Draw();
can_HMtot_1GeV->SaveAs(OutputFolder+TString::Format("HMCF_1GeV_tot_%s.pdf",addon.Data()));

can_HMtot_500MeV = new TCanvas("can_HMtot_500MeV", "can_HMtot_500MeV");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFHMdata[0] ->Draw();
hCFHMdata[1] ->Draw("same");
hCFHMdata[2] ->Draw("same");
leg_HMtot_500MeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HMtot_500MeV->AddEntry((TObject*)0, "HM", "");
leg_HMtot_500MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_500MeV->AddEntry(hCFHMdata[0],"p-#bar{p}");
leg_HMtot_500MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_500MeV->AddEntry(hCFHMdata[1],"p-#bar{#Lambda}");
leg_HMtot_500MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_500MeV->AddEntry(hCFHMdata[2],"#Lambda-#bar{#Lambda}");
leg_HMtot_500MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_500MeV->AddEntry((TObject*)0, st_label, "");
leg_HMtot_500MeV->Draw();
can_HMtot_500MeV->SaveAs(OutputFolder+TString::Format("HMCF_500MeV_tot_%s.pdf",addon.Data()));

can_HMtot_300MeV = new TCanvas("can_HMtot_300MeV", "can_HMtot_300MeV");
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[0] ->GetYaxis()->SetRangeUser(0.1,2.5);
hCFHMdata[0] ->Draw();
hCFHMdata[1] ->Draw("same");
hCFHMdata[2] ->Draw("same");
leg_HMtot_300MeV = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_HMtot_300MeV->AddEntry((TObject*)0, "HM", "");
leg_HMtot_300MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_300MeV->AddEntry(hCFHMdata[0],"p-#bar{p}");
leg_HMtot_300MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_300MeV->AddEntry(hCFHMdata[1],"p-#bar{#Lambda}");
leg_HMtot_300MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_300MeV->AddEntry(hCFHMdata[2],"#Lambda-#bar{#Lambda}");
leg_HMtot_300MeV->AddEntry((TObject*)0, "", "");
leg_HMtot_300MeV->AddEntry((TObject*)0, st_label, "");
leg_HMtot_300MeV->Draw();
can_HMtot_300MeV->SaveAs(OutputFolder+TString::Format("HMCF_300MeV_tot_%s.pdf",addon.Data()));

TCanvas *can_MBHM_3GeV[3];   TCanvas *can_MBHM_1GeV[3];   TCanvas *can_MBHM_500MeV[3];   TCanvas *can_MBHM_300MeV[3];
TLegend *leg_MBHM_3GeV[3];   TLegend *leg_MBHM_1GeV[3];   TLegend *leg_MBHM_500MeV[3];   TLegend *leg_MBHM_300MeV[3];

// p-antip
can_MBHM_3GeV[0] = new TCanvas("can_MBHM_3GeV[0]", "can_MBHM_3GeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.1,3.);
hCFMBdata[0] ->Draw();
hCFHMdata[0] ->Draw("same");
leg_MBHM_3GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_3GeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MBHM_3GeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_MBHM_3GeV[0]->AddEntry((TObject*)0, "", "");
leg_MBHM_3GeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_3GeV[0]->Draw();
can_MBHM_3GeV[0]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_3GeV_pantip_%s.pdf",addon.Data()));

can_MBHM_1GeV[0] = new TCanvas("can_MBHM_1GeV[0]", "can_MBHM_1GeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
hCFHMdata[0] ->Draw("same");
leg_MBHM_1GeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_1GeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MBHM_1GeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_MBHM_1GeV[0]->AddEntry((TObject*)0, "", "");
leg_MBHM_1GeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_1GeV[0]->Draw();
can_MBHM_1GeV[0]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_1GeV_pantip_%s.pdf",addon.Data()));

can_MBHM_500MeV[0] = new TCanvas("can_MBHM_500MeV[0]", "can_MBHM_500MeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
hCFHMdata[0] ->Draw("same");
leg_MBHM_500MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_500MeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MBHM_500MeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_MBHM_500MeV[0]->AddEntry((TObject*)0, "", "");
leg_MBHM_500MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_500MeV[0]->Draw();
can_MBHM_500MeV[0]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_500MeV_pantip_%s.pdf",addon.Data()));

can_MBHM_300MeV[0] = new TCanvas("can_MBHM_300MeV[0]", "can_MBHM_300MeV[0]");
hCFMBdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[0] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata[0] ->Draw();
hCFHMdata[0] ->Draw("same");
leg_MBHM_300MeV[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_300MeV[0]->AddEntry(hCFMBdata[0],"MB p-#bar{p}");
leg_MBHM_300MeV[0]->AddEntry(hCFHMdata[0],"HM p-#bar{p}");
leg_MBHM_300MeV[0]->AddEntry((TObject*)0, "", "");
leg_MBHM_300MeV[0]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_300MeV[0]->Draw();
can_MBHM_300MeV[0]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_300MeV_pantip_%s.pdf",addon.Data()));

//p-antiLambda
can_MBHM_3GeV[1] = new TCanvas("can_MBHM_3GeV[1]", "can_MBHM_3GeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.3,2.2);
hCFMBdata[1] ->Draw();
hCFHMdata[1] ->Draw("same");
leg_MBHM_3GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_3GeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MBHM_3GeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_3GeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_MBHM_3GeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_3GeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_3GeV[1]->Draw();
can_MBHM_3GeV[1]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_3GeV_pantiL_%s.pdf",addon.Data()));

can_MBHM_1GeV[1] = new TCanvas("can_MBHM_1GeV[1]", "can_MBHM_1GeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.3,2.2);
hCFMBdata[1] ->Draw();
hCFHMdata[1] ->Draw("same");
leg_MBHM_1GeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_1GeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MBHM_1GeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_1GeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_MBHM_1GeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_1GeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_1GeV[1]->Draw();
can_MBHM_1GeV[1]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_1GeV_pantiL_%s.pdf",addon.Data()));

can_MBHM_500MeV[1] = new TCanvas("can_MBHM_500MeV[1]", "can_MBHM_500MeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.3,1.2);
hCFMBdata[1] ->Draw();
hCFHMdata[1] ->Draw("same");
leg_MBHM_500MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_500MeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MBHM_500MeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_500MeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_MBHM_500MeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_500MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_500MeV[1]->Draw();
can_MBHM_500MeV[1]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_500MeV_pantiL_%s.pdf",addon.Data()));

can_MBHM_300MeV[1] = new TCanvas("can_MBHM_300MeV[1]", "can_MBHM_300MeV[1]");
hCFMBdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[1] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[1] ->GetYaxis()->SetRangeUser(0.3,1.8);
hCFMBdata[1] ->Draw();
hCFHMdata[1] ->Draw("same");
leg_MBHM_300MeV[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_300MeV[1]->AddEntry(hCFMBdata[1],"MB p-#bar{#Lambda}");
leg_MBHM_300MeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_300MeV[1]->AddEntry(hCFHMdata[1],"HM p-#bar{#Lambda}");
leg_MBHM_300MeV[1]->AddEntry((TObject*)0, "", "");
leg_MBHM_300MeV[1]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_300MeV[1]->Draw();
can_MBHM_300MeV[1]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_300MeV_pantiL_%s.pdf",addon.Data()));

//Lambda-antiLambda
can_MBHM_3GeV[2] = new TCanvas("can_MBHM_3GeV[2]", "can_MBHM_3GeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,3000);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.3,2.);
hCFMBdata[2] ->Draw();
hCFHMdata[2] ->Draw("same");
leg_MBHM_3GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_3GeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MBHM_3GeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_3GeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_MBHM_3GeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_3GeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_3GeV[2]->Draw();
can_MBHM_3GeV[2]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_3GeV_LantiL_%s.pdf",addon.Data()));

can_MBHM_1GeV[2] = new TCanvas("can_MBHM_1GeV[2]", "can_MBHM_1GeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,1000);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.3,2.);
hCFMBdata[2] ->Draw();
hCFHMdata[2] ->Draw("same");
leg_MBHM_1GeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_1GeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MBHM_1GeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_1GeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_MBHM_1GeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_1GeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_1GeV[2]->Draw();
can_MBHM_1GeV[2]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_1GeV_LantiL_%s.pdf",addon.Data()));

can_MBHM_500MeV[2] = new TCanvas("can_MBHM_500MeV[2]", "can_MBHM_500MeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,500);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.3,1.8);
hCFMBdata[2] ->Draw();
hCFHMdata[2] ->Draw("same");
leg_MBHM_500MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_500MeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MBHM_500MeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_500MeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_MBHM_500MeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_500MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_500MeV[2]->Draw();
can_MBHM_500MeV[2]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_500MeV_LantiL_%s.pdf",addon.Data()));

can_MBHM_300MeV[2] = new TCanvas("can_MBHM_300MeV[2]", "can_MBHM_300MeV[2]");
hCFMBdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFHMdata[2] ->GetXaxis()->SetRangeUser(0.,300);
hCFMBdata[2] ->GetYaxis()->SetRangeUser(0.3,2.);
hCFMBdata[2] ->Draw();
hCFHMdata[2] ->Draw("same");
leg_MBHM_300MeV[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
leg_MBHM_300MeV[2]->AddEntry(hCFMBdata[2],"MB #Lambda-#bar{#Lambda}");
leg_MBHM_300MeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_300MeV[2]->AddEntry(hCFHMdata[2],"HM #Lambda-#bar{#Lambda}");
leg_MBHM_300MeV[2]->AddEntry((TObject*)0, "", "");
leg_MBHM_300MeV[2]->AddEntry((TObject*)0, st_label, "");
leg_MBHM_300MeV[2]->Draw();
can_MBHM_300MeV[2]->SaveAs(OutputFolder+TString::Format("CompHMMBCF_300MeV_LantiL_%s.pdf",addon.Data()));

//Reading the root file
TFile* filecomp=TFile::Open(OutputFolder+"ComparisonHMMB.root");

TH1F* hCFMBdata_st1[3];TH1F* hCFHMdata_st1[3];
TH1F* hCFMBdata_st2[3];TH1F* hCFHMdata_st2[3];
TH1F* hCFMBdata_st3[3];TH1F* hCFHMdata_st3[3];
TH1F* hCFMBdata_st4[3];TH1F* hCFHMdata_st4[3];

for(int i=0;i<3;i++){
  hCFMBdata_st1[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFMBdata_%i_1",i)));
  hCFMBdata_st1[i]->Sumw2();
  hCFMBdata_st2[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFMBdata_%i_2",i)));
  hCFMBdata_st2[i]->Sumw2();
  hCFMBdata_st3[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFMBdata_%i_3",i)));
  hCFMBdata_st3[i]->Sumw2();
  hCFMBdata_st4[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFMBdata_%i_4",i)));
  hCFMBdata_st4[i]->Sumw2();

  hCFHMdata_st1[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFHMdata_%i_1",i)));
  hCFHMdata_st1[i]->Sumw2();
  hCFHMdata_st2[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFHMdata_%i_2",i)));
  hCFHMdata_st2[i]->Sumw2();
  hCFHMdata_st3[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFHMdata_%i_3",i)));
  hCFHMdata_st3[i]->Sumw2();
  hCFHMdata_st4[i] = (TH1F*)(filecomp->FindObjectAny(TString::Format("hCFHMdata_%i_4",i)));
  hCFHMdata_st4[i]->Sumw2();
}

TCanvas *can_MBtotST[3];
TLegend *leg_MBtotST[3];
TCanvas *can_HMtotST[3];
TLegend *leg_HMtotST[3];

double suckmycock = 0.02;

for(int i=0; i<3; i++){
DreamPlot::SetStyleHistoCF(hCFMBdata_st1[i], 8, kRed, 20);//st1 red
DreamPlot::SetStyleHistoCF(hCFMBdata_st2[i], 8, kMagenta, 20);//st2 magenta
DreamPlot::SetStyleHistoCF(hCFMBdata_st3[i], 8, kGreen, 20);//st3 green
DreamPlot::SetStyleHistoCF(hCFMBdata_st4[i], 8, kBlue, 20);//st4 blue
DreamPlot::SetStyleHistoCF(hCFHMdata_st1[i], 8, kRed, 20);//st1 red
DreamPlot::SetStyleHistoCF(hCFHMdata_st2[i], 8, kMagenta, 20);//st2 magenta
DreamPlot::SetStyleHistoCF(hCFHMdata_st3[i], 8, kGreen, 20);//st3 green
DreamPlot::SetStyleHistoCF(hCFHMdata_st4[i], 8, kBlue, 20);//st4 blue
}

can_MBtotST[0] = new TCanvas("can_MBtotST[0]", "can_MBtotST[0]");
hCFMBdata_st1[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st2[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st3[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st4[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st1[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFMBdata_st1[0] ->Draw();
hCFMBdata_st2[0] ->Draw("same");
hCFMBdata_st3[0] ->Draw("same");
hCFMBdata_st4[0] ->Draw("same");
leg_MBtotST[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_MBtotST[0]->AddEntry((TObject*)0, "MB", "");
leg_MBtotST[0]->AddEntry((TObject*)0, "", "");
leg_MBtotST[0]->AddEntry(hCFMBdata_st1[0],"p-#bar{p} 0 < s_{T} < 0.3");
leg_MBtotST[0]->AddEntry((TObject*)0, "", "");
leg_MBtotST[0]->AddEntry(hCFMBdata_st2[0],"p-#bar{p} 0.3 < s_{T} < 0.7");
leg_MBtotST[0]->AddEntry((TObject*)0, "", "");
leg_MBtotST[0]->AddEntry(hCFMBdata_st3[0],"p-#bar{p} 0.7 < s_{T} < 1");
leg_MBtotST[0]->AddEntry((TObject*)0, "", "");
leg_MBtotST[0]->AddEntry(hCFMBdata_st4[0],"p-#bar{p} 0 < s_{T} < 1");
leg_MBtotST[0]->Draw();
can_MBtotST[0]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_MBCF_pantip_%s.pdf",maxkstar));

can_MBtotST[1] = new TCanvas("can_MBtotST[1]", "can_MBtotST[1]");
hCFMBdata_st1[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st2[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st3[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st4[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st1[1] ->GetYaxis()->SetRangeUser(0.1,1.5);
hCFMBdata_st1[1] ->Draw();
hCFMBdata_st2[1] ->Draw("same");
hCFMBdata_st3[1] ->Draw("same");
hCFMBdata_st4[1] ->Draw("same");
leg_MBtotST[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_MBtotST[1]->AddEntry((TObject*)0, "MB", "");
leg_MBtotST[1]->AddEntry((TObject*)0, "", "");
leg_MBtotST[1]->AddEntry(hCFMBdata_st1[1],"p-#bar{#Lambda} 0 < s_{T} < 0.3");
leg_MBtotST[1]->AddEntry((TObject*)0, "", "");
leg_MBtotST[1]->AddEntry(hCFMBdata_st2[1],"p-#bar{#Lambda} 0.3 < s_{T} < 0.7");
leg_MBtotST[1]->AddEntry((TObject*)0, "", "");
leg_MBtotST[1]->AddEntry(hCFMBdata_st3[1],"p-#bar{#Lambda} 0.7 < s_{T} < 1");
leg_MBtotST[1]->AddEntry((TObject*)0, "", "");
leg_MBtotST[1]->AddEntry(hCFMBdata_st4[1],"p-#bar{#Lambda} 0 < s_{T} < 1");
leg_MBtotST[1]->Draw();
can_MBtotST[1]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_MBCF_pantiL_%s.pdf",maxkstar));

can_MBtotST[2] = new TCanvas("can_MBtotST[2]", "can_MBtotST[2]");
hCFMBdata_st1[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st2[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st3[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st4[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFMBdata_st1[2] ->GetYaxis()->SetRangeUser(0.1,1.5);
hCFMBdata_st1[2] ->Draw();
hCFMBdata_st2[2] ->Draw("same");
hCFMBdata_st3[2] ->Draw("same");
hCFMBdata_st4[2] ->Draw("same");
leg_MBtotST[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_MBtotST[2]->AddEntry((TObject*)0, "MB", "");
leg_MBtotST[2]->AddEntry((TObject*)0, "", "");
leg_MBtotST[2]->AddEntry(hCFMBdata_st1[2],"#Lambda-#bar{#Lambda} 0 < s_{T} < 0.3");
leg_MBtotST[2]->AddEntry((TObject*)0, "", "");
leg_MBtotST[2]->AddEntry(hCFMBdata_st2[2],"#Lambda-#bar{#Lambda} 0.3 < s_{T} < 0.7");
leg_MBtotST[2]->AddEntry((TObject*)0, "", "");
leg_MBtotST[2]->AddEntry(hCFMBdata_st3[2],"#Lambda-#bar{#Lambda} 0.7 < s_{T} < 1");
leg_MBtotST[2]->AddEntry((TObject*)0, "", "");
leg_MBtotST[2]->AddEntry(hCFMBdata_st4[2],"#Lambda-#bar{#Lambda} 0 < s_{T} < 1");
leg_MBtotST[2]->Draw();
can_MBtotST[2]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_MBCF_LantiL_%s.pdf",maxkstar));




can_HMtotST[0] = new TCanvas("can_HMtotST[0]", "can_HMtotST[0]");
hCFHMdata_st1[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st2[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st3[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st4[0] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st1[0] ->GetYaxis()->SetRangeUser(0.3,3.);
hCFHMdata_st1[0] ->Draw();
hCFHMdata_st2[0] ->Draw("same");
hCFHMdata_st3[0] ->Draw("same");
hCFHMdata_st4[0] ->Draw("same");
leg_HMtotST[0] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_HMtotST[0]->AddEntry((TObject*)0, "HM", "");
leg_HMtotST[0]->AddEntry((TObject*)0, "", "");
leg_HMtotST[0]->AddEntry(hCFHMdata_st1[0],"p-#bar{p} 0 < s_{T} < 0.3");
leg_HMtotST[0]->AddEntry((TObject*)0, "", "");
leg_HMtotST[0]->AddEntry(hCFHMdata_st2[0],"p-#bar{p} 0.3 < s_{T} < 0.7");
leg_HMtotST[0]->AddEntry((TObject*)0, "", "");
leg_HMtotST[0]->AddEntry(hCFHMdata_st3[0],"p-#bar{p} 0.7 < s_{T} < 1");
leg_HMtotST[0]->AddEntry((TObject*)0, "", "");
leg_HMtotST[0]->AddEntry(hCFHMdata_st4[0],"p-#bar{p} 0 < s_{T} < 1");
leg_HMtotST[0]->Draw();
can_HMtotST[0]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_HMCF_pantip_%s.pdf",maxkstar));

can_HMtotST[1] = new TCanvas("can_HMtotST[1]", "can_HMtotST[1]");
hCFHMdata_st1[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st2[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st3[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st4[1] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st1[1] ->GetYaxis()->SetRangeUser(0.1,1.5);
hCFHMdata_st1[1] ->Draw();
hCFHMdata_st2[1] ->Draw("same");
hCFHMdata_st3[1] ->Draw("same");
hCFHMdata_st4[1] ->Draw("same");
leg_HMtotST[1] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_HMtotST[1]->AddEntry((TObject*)0, "HM", "");
leg_HMtotST[1]->AddEntry((TObject*)0, "", "");
leg_HMtotST[1]->AddEntry(hCFHMdata_st1[1],"p-#bar{#Lambda} 0 < s_{T} < 0.3");
leg_HMtotST[1]->AddEntry((TObject*)0, "", "");
leg_HMtotST[1]->AddEntry(hCFHMdata_st2[1],"p-#bar{#Lambda} 0.3 < s_{T} < 0.7");
leg_HMtotST[1]->AddEntry((TObject*)0, "", "");
leg_HMtotST[1]->AddEntry(hCFHMdata_st3[1],"p-#bar{#Lambda} 0.7 < s_{T} < 1");
leg_HMtotST[1]->AddEntry((TObject*)0, "", "");
leg_HMtotST[1]->AddEntry(hCFHMdata_st4[1],"p-#bar{#Lambda} 0 < s_{T} < 1");
leg_HMtotST[1]->Draw();
can_HMtotST[1]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_HMCF_pantiL_%s.pdf",maxkstar));

can_HMtotST[2] = new TCanvas("can_HMtotST[2]", "can_HMtotST[2]");
hCFHMdata_st1[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st2[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st3[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st4[2] ->GetXaxis()->SetRangeUser(0.,charmask);
hCFHMdata_st1[2] ->GetYaxis()->SetRangeUser(0.1,1.5);
hCFHMdata_st1[2] ->Draw();
hCFHMdata_st2[2] ->Draw("same");
hCFHMdata_st3[2] ->Draw("same");
hCFHMdata_st4[2] ->Draw("same");
leg_HMtotST[2] = new TLegend(0.6, 0.7, 0.88, 0.8);
gStyle->SetLegendTextSize(suckmycock);
leg_HMtotST[2]->AddEntry((TObject*)0, "HM", "");
leg_HMtotST[2]->AddEntry((TObject*)0, "", "");
leg_HMtotST[2]->AddEntry(hCFHMdata_st1[2],"#Lambda-#bar{#Lambda} 0 < s_{T} < 0.3");
leg_HMtotST[2]->AddEntry((TObject*)0, "", "");
leg_HMtotST[2]->AddEntry(hCFHMdata_st2[2],"#Lambda-#bar{#Lambda} 0.3 < s_{T} < 0.7");
leg_HMtotST[2]->AddEntry((TObject*)0, "", "");
leg_HMtotST[2]->AddEntry(hCFHMdata_st3[2],"#Lambda-#bar{#Lambda} 0.7 < s_{T} < 1");
leg_HMtotST[2]->AddEntry((TObject*)0, "", "");
leg_HMtotST[2]->AddEntry(hCFHMdata_st4[2],"#Lambda-#bar{#Lambda} 0 < s_{T} < 1");
leg_HMtotST[2]->Draw();
can_HMtotST[2]->SaveAs(OutputFolder+TString::Format("Comp_sTclasses_HMCF_LantiL_%s.pdf",maxkstar));

return 0;
}
