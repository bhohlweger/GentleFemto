#include "DLM_Source.h"
#include "TLatex.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include "stdlib.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "SideBandFit.h"
#include "TMinuit.h"
#include "TDatabasePDG.h"
#include "TLegend.h"
#include <vector>
#include "TColor.h"
#include "TStyle.h"
#include "TApplication.h"
#include "stdlib.h"
#include <chrono>
#include <ctime>
#include <numeric>
#include "TLine.h"

static double radCutoff = 0;

static std::vector<int> fFillColors = { kGray + 1, kRed - 10, kBlue - 9, kGreen
    - 8, kMagenta - 9, kOrange - 9, kCyan - 3, kYellow - 7 };
static std::vector<int> fColors = { kBlack, kRed + 1,         //1
kBlue + 2,         //2
kGreen + 3,         //3
kMagenta - 8,         //4
kOrange - 3,         //5
kCyan + 2,         //6
kYellow + 2,         //7
kWhite,         //8
    kGreen - 5,         //9
};
static std::vector<int> fMarkers = { kFullCircle, kFullSquare, kOpenCircle,
    kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar,
    kOpenStar };

void SetStyleHisto(TH1 *histo) {

  histo->GetXaxis()->SetLabelFont(42);

  histo->GetXaxis()->SetLabelSize(0.11);
  histo->GetYaxis()->SetLabelSize(0.11);
  histo->GetXaxis()->SetTitleSize(0.11);
  histo->GetYaxis()->SetTitleSize(0.11);

  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(0.95);
  histo->GetYaxis()->SetTitleOffset(0.6);

  histo->SetMarkerSize(1.2);
  histo->SetLineWidth(3);
}

void SetStyleGraph(TGraph *histo) {
  histo->GetXaxis()->SetLabelFont(42);

  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetYaxis()->SetTitleOffset(1.0);

  histo->GetXaxis()->SetLabelSize(0.06);
  histo->GetYaxis()->SetLabelSize(0.06);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetYaxis()->SetTitleSize(0.06);
}

void PlotPotentials() {
//  TFile *CFs = TFile::Open("InletGoodCFs.root", "read");
  TFile *CFs = TFile::Open("InletCFs.root", "read");
  TH1F* CF_I0S0 = (TH1F*) CFs->Get("I0S0");
  CF_I0S0->SetLineStyle(7);
  TH1F* CF_I0S1 = (TH1F*) CFs->Get("I0S1");
  CF_I0S1->SetLineStyle(3);
  TH1F* CF_I1S0 = (TH1F*) CFs->Get("I1S0");
  CF_I1S0->SetLineStyle(1);
  TH1F* CF_I1S1 = (TH1F*) CFs->Get("I1S1");
  CF_I1S1->SetLineStyle(4);

  TFile *OutFile = TFile::Open("PotentialXi.root", "RECREATE");
  TH1F* I0S0[3];
  TList *ListI0S0 = new TList();
  ListI0S0->SetOwner();
  ListI0S0->SetName("I0S0");

  TH1F* I0S1[3];
  TList *ListI0S1 = new TList();
  ListI0S1->SetOwner();
  ListI0S1->SetName("I0S1");

  TH1F* I1S0[3];
  TList *ListI1S0 = new TList();
  ListI1S0->SetOwner();
  ListI1S0->SetName("I1S0");

  TH1F* I1S1[3];
  TList *ListI1S1 = new TList();
  ListI1S1->SetOwner();
  ListI1S1->SetName("I1S1");

  TLegend leg(0.55, 0.65, 0.95, 0.95);
  leg.SetLineColor(kWhite);
  leg.SetFillStyle(4000);
  leg.SetTextSize(0.07);

  TLegend leg2(0.55, 0.65, 0.95, 0.95);
  leg2.SetLineColor(kWhite);
  leg2.SetFillStyle(4000);
  leg2.SetTextSize(0.07);

  int nRadBins = 600;
  double dRad = 0.005;
  double radMax = dRad * nRadBins;
  double times[3] = {0, 23, 46}; 
  for (int iTime = 0; iTime < 3; ++iTime) {
    double QCDTime = times[iTime]; 

    TString I0S0Name = Form("I0S0_Time%u", iTime);
    I0S0[iTime] = new TH1F(I0S0Name.Data(), I0S0Name.Data(), nRadBins,
                           dRad / 2., radMax + dRad / 2.);

    TString I0S1Name = Form("I0S1_Time%u", iTime);
    I0S1[iTime] = new TH1F(I0S1Name.Data(), I0S1Name.Data(), nRadBins,
                           dRad / 2., radMax + dRad / 2.);

    TString I1S0Name = Form("I1S0_Time%u", iTime);
    I1S0[iTime] = new TH1F(I1S0Name.Data(), I1S0Name.Data(), nRadBins,
                           dRad / 2., radMax + dRad / 2.);

    TString I1S1Name = Form("I1S1_Time%u", iTime);
    I1S1[iTime] = new TH1F(I1S1Name.Data(), I1S1Name.Data(), nRadBins,
                           dRad / 2., radMax + dRad / 2.);

    for (int iRad = 1; iRad < nRadBins; ++iRad) {
      double pXimPotParsI0S0[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 0,
          -1, 1, 0, 0, 0};
      double pXimPotParsI0S1[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 0,
          -1, 1, 1, 0, 1};
      double pXimPotParsI1S0[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 1,
          1, 1, 0, 0, 0 };
      double pXimPotParsI1S1[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 1,
          1, 1, 1, 0, 1 };

      I0S0[iTime]->SetBinContent(I0S0[iTime]->FindBin(iRad * dRad),
                                 fDlmPot(pXimPotParsI0S0));

      I0S1[iTime]->SetBinContent(I0S1[iTime]->FindBin(iRad * dRad),
                                 fDlmPot(pXimPotParsI0S1));

      I1S0[iTime]->SetBinContent(I1S0[iTime]->FindBin(iRad * dRad),
                                 fDlmPot(pXimPotParsI1S0));

      I1S1[iTime]->SetBinContent(I1S1[iTime]->FindBin(iRad * dRad),
                                 fDlmPot(pXimPotParsI1S1));
    }
    OutFile->cd();
    ListI0S0->Add(I0S0[iTime]);
    ListI0S1->Add(I0S1[iTime]);
    ListI1S0->Add(I1S0[iTime]);
    ListI1S1->Add(I1S1[iTime]);
  }
  TGraphErrors *grI0S0 = new TGraphErrors();
  TGraph*grI0S0Min = new TGraph();
  TGraph*grI0S0Max = new TGraph();

  TGraphErrors *grI0S1 = new TGraphErrors();
  TGraph*grI0S1Min = new TGraph();
  TGraph*grI0S1Max = new TGraph();

  TGraphErrors *grI1S0 = new TGraphErrors();
  TGraph*grI1S0Min = new TGraph();
  TGraph*grI1S0Max = new TGraph();

  TGraphErrors *grI1S1 = new TGraphErrors();
  TGraph*grI1S1Min = new TGraph();
  TGraph*grI1S1Max = new TGraph();

  for (int iRad = 1; iRad < nRadBins; ++iRad) {
    std::vector<float> I0S0Val;
    std::vector<float> I0S1Val;
    std::vector<float> I1S0Val;
    std::vector<float> I1S1Val;
    for (int iTime = 0; iTime < 3; ++iTime) {
      I0S0Val.push_back(I0S0[iTime]->GetBinContent(iRad));
      I0S1Val.push_back(I0S1[iTime]->GetBinContent(iRad));
      I1S0Val.push_back(I1S0[iTime]->GetBinContent(iRad));
      I1S1Val.push_back(I1S1[iTime]->GetBinContent(iRad));
    }
    std::sort(I0S0Val.begin(), I0S0Val.end());

    double DiffI0S0 = TMath::Abs((I0S0Val[0] - I0S0Val[2]) / 2.);
    double potValI0S0 = I0S0Val[0] + DiffI0S0;

    grI0S0->SetPoint(iRad - 1, I0S0[0]->GetBinCenter(iRad), potValI0S0);
    grI0S0->SetPointError(iRad - 1, 0, DiffI0S0);
    grI0S0Min->SetPoint(iRad - 1, I0S0[0]->GetBinCenter(iRad), I0S0Val[0]);
    grI0S0Max->SetPoint(iRad - 1, I0S0[0]->GetBinCenter(iRad), I0S0Val[2]);

    std::sort(I0S1Val.begin(), I0S1Val.end());

    double DiffI0S1 = TMath::Abs((I0S1Val[0] - I0S1Val[2]) / 2.);
    double potValI0S1 = I0S1Val[0] + DiffI0S1;

    grI0S1->SetPoint(iRad - 1, I0S1[0]->GetBinCenter(iRad), potValI0S1);
    grI0S1->SetPointError(iRad - 1, 0, DiffI0S1);
    grI0S1Min->SetPoint(iRad - 1, I0S1[0]->GetBinCenter(iRad), I0S1Val[0]);
    grI0S1Max->SetPoint(iRad - 1, I0S1[0]->GetBinCenter(iRad), I0S1Val[2]);

    std::sort(I1S0Val.begin(), I1S0Val.end());

    double DiffI1S0 = TMath::Abs((I1S0Val[0] - I1S0Val[2]) / 2.);
    double potValI1S0 = I1S0Val[0] + DiffI1S0;

    grI1S0->SetPoint(iRad - 1, I1S0[0]->GetBinCenter(iRad), potValI1S0);
    grI1S0->SetPointError(iRad - 1, 0, DiffI1S0);
    grI1S0Min->SetPoint(iRad - 1, I1S0[0]->GetBinCenter(iRad), I1S0Val[0]);
    grI1S0Max->SetPoint(iRad - 1, I1S0[0]->GetBinCenter(iRad), I1S0Val[2]);

    std::sort(I1S1Val.begin(), I1S1Val.end());

    double DiffI1S1 = TMath::Abs((I1S1Val[0] - I1S1Val[2]) / 2.);
    double potValI1S1 = I1S1Val[0] + DiffI1S1;

    grI1S1->SetPoint(iRad - 1, I1S1[0]->GetBinCenter(iRad), potValI1S1);
    grI1S1->SetPointError(iRad - 1, 0, DiffI1S1);
    grI1S1Min->SetPoint(iRad - 1, I1S1[0]->GetBinCenter(iRad), I1S1Val[0]);
    grI1S1Max->SetPoint(iRad - 1, I1S1[0]->GetBinCenter(iRad), I1S1Val[2]);

  }
  grI0S0->SetName("grI0S0");
  grI0S0->GetXaxis()->SetTitle("#it{r} (fm)");
  grI0S0->GetYaxis()->SetRangeUser(-75, 59);
  grI0S0->GetXaxis()->SetRangeUser(0, 3);
  grI0S0->GetYaxis()->SetTitle("#it{V(r)} (MeV)");
//  grI0S0->GetYaxis()->CenterTitle(true);

  grI0S0->SetFillColor(fColors[1]);
  grI0S0->SetLineColor(fColors[1]);
  grI0S0->SetFillStyle(3013);

  grI0S0Min->SetLineColor(fColors[1]);
//  grI0S0Min->SetLineStyle(7);
  grI0S0Min->SetName("grI0S0Min");
  grI0S0Max->SetLineColor(fColors[1]);
//  grI0S0Max->SetLineStyle(7);
  grI0S0Max->SetName("grI0S0Max");

  ListI0S0->Add(grI0S0);
  ListI0S0->Add(grI0S0Min);
  ListI0S0->Add(grI0S0Max);
  ListI0S0->Write("I0S0", 1);

  grI0S1->SetName("grI0S1");
  grI0S1->SetLineColor(fColors[2]);
  grI0S1->SetFillColor(fColors[2]);
  grI0S1->SetFillStyle(3002);

  grI0S1Min->SetLineColor(fColors[2]);
  grI0S1Min->SetName("grI0S1Min");
//  grI0S1Min->SetLineStyle(3);
  grI0S1Max->SetLineColor(fColors[2]);
//  grI0S1Max->SetLineStyle(3);
  grI0S1Max->SetName("grI0S1Max");

  ListI0S1->Add(grI0S1);
  ListI0S1->Add(grI0S1Min);
  ListI0S1->Add(grI0S1Max);
  ListI0S1->Write("I0S1", 1);

  grI1S0->SetName("grI1S0");
  grI1S0->SetLineColor(fColors[3]);
  grI1S0->SetFillColor(fColors[3]);
  grI1S0->SetFillStyle(3004);

  grI1S0Min->SetLineColor(fColors[3]);
  grI1S0Min->SetName("grI1S0Min");
//  grI1S0Min->SetLineStyle(1);
  grI1S0Max->SetLineColor(fColors[3]);
  grI1S0Max->SetName("grI1S0Max");
//  grI1S0Max->SetLineStyle(1);

  ListI1S0->Add(grI1S0);
  ListI1S0->Add(grI1S0Min);
  ListI1S0->Add(grI1S0Max);
  ListI1S0->Write("I1S0", 1);

  grI1S1->SetName("grI1S1");
  grI1S1->SetLineColor(fColors[5]);
  grI1S1->SetFillColor(fColors[5]);
  grI1S1->SetFillStyle(3005);

  grI1S1Min->SetLineColor(fColors[5]);
  grI1S1Min->SetName("grI1S1Min");
//  grI1S1Min->SetLineStyle(4);
  grI1S1Max->SetLineColor(fColors[5]);
  grI1S1Max->SetName("grI1S1Max");
//  grI1S1Max->SetLineStyle(4);

  ListI1S1->Add(grI1S1);
  ListI1S1->Add(grI1S1Min);
  ListI1S1->Add(grI1S1Max);
  ListI1S1->Write("I1S1", 1);

  auto *c2 = new TCanvas("c2", "c2");
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.01);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.15);
  c2->cd();
  leg.AddEntry(grI0S0, "I = 0, S = 0", "f");
  leg.AddEntry(grI0S1, "I = 0, S = 1", "f");
  leg.AddEntry(grI1S0, "I = 1, S = 0", "f");
  leg.AddEntry(grI1S1, "I = 1, S = 1", "f");

  SetStyleGraph(grI0S0);
  SetStyleGraph(grI0S0Min);
  SetStyleGraph(grI0S0Max);

  SetStyleGraph(grI0S1);
  SetStyleGraph(grI0S1Min);
  SetStyleGraph(grI0S1Max);

  SetStyleGraph(grI1S0);
  SetStyleGraph(grI1S0Min);
  SetStyleGraph(grI1S0Max);

  SetStyleGraph(grI1S1);
  SetStyleGraph(grI1S1Min);
  SetStyleGraph(grI1S1Max);

  grI0S0->Draw("az3");
  grI0S0Min->Draw("lsame");
  grI0S0Max->Draw("lsame");

  grI0S1->Draw("3zsame");
  grI0S1Min->Draw("lsame");
  grI0S1Max->Draw("lsame");

  grI1S0->Draw("3zsame");
  grI1S0Min->Draw("lsame");
  grI1S0Max->Draw("lsame");

  grI1S1->Draw("3zsame");
  grI1S1Min->Draw("lsame");
  grI1S1Max->Draw("lsame");

  TLine lineOne = TLine(0, 0, 3., 0);
  lineOne.Draw("SAME");

  TLine lineTwo = TLine(0, 1., 175, 1.);
  lineOne.Draw("SAME");
  TPad *inset_pad = new TPad("inset", "inset", 0.5, 0.13, 0.97, 0.63);
  inset_pad->SetRightMargin(0.05);
  inset_pad->SetLeftMargin(0.15);
  inset_pad->SetBottomMargin(0.3);
  inset_pad->SetFillStyle(4000);
  inset_pad->Draw();
  inset_pad->cd();

  SetStyleHisto(CF_I0S0);
  SetStyleHisto(CF_I0S1);
  SetStyleHisto(CF_I1S0);
  SetStyleHisto(CF_I1S1);
  CF_I0S0->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_I0S0->GetYaxis()->CenterTitle(true);
  CF_I0S0->GetXaxis()->SetRangeUser(0, 175);
  CF_I0S0->GetYaxis()->SetRangeUser(0.8, 3.5);
  CF_I0S0->SetNdivisions(505);
  CF_I0S0->GetYaxis()->SetNdivisions(203);
  CF_I0S0->Draw("l");
  lineTwo.Draw("same");
  CF_I0S1->Draw("lsame");
  CF_I1S0->Draw("lsame");
  CF_I1S1->Draw("lsame");
  TLatex Numbering;
  Numbering.SetTextSize(gStyle->GetTextSize() * 2.0);
  Numbering.SetNDC(kTRUE);
  Numbering.DrawLatex(0.6, 0.8, "#it{r}_{0} = 1.4 fm");
  leg2.AddEntry(CF_I0S0, " ", "l");
  leg2.AddEntry(CF_I0S1, " ", "l");
  leg2.AddEntry(CF_I1S0, " ", "l");
  leg2.AddEntry(CF_I1S1, " ", "l");

  c2->cd();
  leg2.Draw("SAME");
  leg.Draw("same");
  OutFile->cd();
  c2->Write();
  c2->SaveAs("c2.pdf");
  OutFile->Close();
}

void JackknifePotentials() {
  TFile *OutFile = TFile::Open("PotentialXi.root", "RECREATE");
  OutFile->cd();

  TList *ListI0S0 = new TList();
  ListI0S0->SetOwner();
  ListI0S0->SetName("I0S0");

  TH1F* I0S1[3];
  TList *ListI0S1 = new TList();
  ListI0S1->SetOwner();
  ListI0S1->SetName("I0S1");

  TH1F* I1S0[3];
  TList *ListI1S0 = new TList();
  ListI1S0->SetOwner();
  ListI1S0->SetName("I1S0");

  TH1F* I1S1[3];
  TList *ListI1S1 = new TList();
  ListI1S1->SetOwner();
  ListI1S1->SetName("I1S1");

  int nRadBins = 1200;
  double dRad = 0.005;
  double radMax = dRad * nRadBins;

  std::vector<std::vector<float>> vecRad_I0S0_t11;
  std::vector<std::vector<float>> vecRad_I0S0_t12;
  std::vector<std::vector<float>> vecRad_I0S0_t13;

  std::vector<std::vector<float>> vecRad_I0S1_t11;
  std::vector<std::vector<float>> vecRad_I0S1_t12;
  std::vector<std::vector<float>> vecRad_I0S1_t13;
	
  std::vector<std::vector<float>> vecRad_I1S0_t11;
  std::vector<std::vector<float>> vecRad_I1S0_t12;
  std::vector<std::vector<float>> vecRad_I1S0_t13;

  std::vector<std::vector<float>> vecRad_I1S1_t11;
  std::vector<std::vector<float>> vecRad_I1S1_t12;
  std::vector<std::vector<float>> vecRad_I1S1_t13;   


  for (int iRad = 1; iRad < nRadBins; ++iRad) {
    std::vector<float> vec_I0S0_t11;
    std::vector<float> vec_I0S0_t12;
    std::vector<float> vec_I0S0_t13;
    
    std::vector<float> vec_I0S1_t11;
    std::vector<float> vec_I0S1_t12;
    std::vector<float> vec_I0S1_t13;
    
    std::vector<float> vec_I1S0_t11;
    std::vector<float> vec_I1S0_t12;
    std::vector<float> vec_I1S0_t13;
    
    std::vector<float> vec_I1S1_t11;
    std::vector<float> vec_I1S1_t12;
    std::vector<float> vec_I1S1_t13;
    
    for (int iTime = 0; iTime < 71; ++iTime) {
      if (iTime == 0 || iTime == 24 || iTime == 48) {
	//Default obtained by average 
	continue; 
      } 
      double QCDTime = iTime; 
      double pXimPotParsI0S0[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 0,
				     -1, 1, 0, 0, 0};
      double pXimPotParsI0S1[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 0,
				     -1, 1, 1, 0, 1};
      double pXimPotParsI1S0[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 1,
				     1, 1, 0, 0, 0 };
      double pXimPotParsI1S1[10] = { iRad * dRad, 0, pXim_HALQCDPaper2020, QCDTime, 1,
				     1, 1, 1, 0, 1 };
      if (iTime < 24) { 
	vec_I0S0_t11.push_back(fDlmPot(pXimPotParsI0S0));
	vec_I0S1_t11.push_back(fDlmPot(pXimPotParsI0S1));
	vec_I1S0_t11.push_back(fDlmPot(pXimPotParsI1S0));
	vec_I1S1_t11.push_back(fDlmPot(pXimPotParsI1S1)); 
      } else if (iTime < 48) {
	vec_I0S0_t12.push_back(fDlmPot(pXimPotParsI0S0));
	vec_I0S1_t12.push_back(fDlmPot(pXimPotParsI0S1));
	vec_I1S0_t12.push_back(fDlmPot(pXimPotParsI1S0));
	vec_I1S1_t12.push_back(fDlmPot(pXimPotParsI1S1)); 
      } else {
	vec_I0S0_t13.push_back(fDlmPot(pXimPotParsI0S0));
	vec_I0S1_t13.push_back(fDlmPot(pXimPotParsI0S1));
	vec_I1S0_t13.push_back(fDlmPot(pXimPotParsI1S0));
	vec_I1S1_t13.push_back(fDlmPot(pXimPotParsI1S1)); 
      }
    }

    std::sort(vec_I0S0_t11.begin(),vec_I0S0_t11.end(), std::greater<float>()); 
    std::sort(vec_I0S0_t12.begin(),vec_I0S0_t12.end(), std::greater<float>()); 
    std::sort(vec_I0S0_t13.begin(),vec_I0S0_t13.end(), std::greater<float>()); 
    
    std::sort(vec_I0S1_t11.begin(),vec_I0S1_t11.end(), std::greater<float>()); 
    std::sort(vec_I0S1_t12.begin(),vec_I0S1_t12.end(), std::greater<float>()); 
    std::sort(vec_I0S1_t13.begin(),vec_I0S1_t13.end(), std::greater<float>()); 

    std::sort(vec_I1S0_t11.begin(),vec_I1S0_t11.end(), std::greater<float>()); 
    std::sort(vec_I1S0_t12.begin(),vec_I1S0_t12.end(), std::greater<float>()); 
    std::sort(vec_I1S0_t13.begin(),vec_I1S0_t13.end(), std::greater<float>()); 

    std::sort(vec_I1S1_t11.begin(),vec_I1S1_t11.end(), std::greater<float>()); 
    std::sort(vec_I1S1_t12.begin(),vec_I1S1_t12.end(), std::greater<float>()); 
    std::sort(vec_I1S1_t13.begin(),vec_I1S1_t13.end(), std::greater<float>()); 

    vecRad_I0S0_t11.push_back(vec_I0S0_t11);
    vecRad_I0S0_t12.push_back(vec_I0S0_t12);
    vecRad_I0S0_t13.push_back(vec_I0S0_t13);
                   	                  
    vecRad_I0S1_t11.push_back(vec_I0S1_t11);
    vecRad_I0S1_t12.push_back(vec_I0S1_t12);
    vecRad_I0S1_t13.push_back(vec_I0S1_t13);
                   	                  
    vecRad_I1S0_t11.push_back(vec_I1S0_t11);
    vecRad_I1S0_t12.push_back(vec_I1S0_t12);
    vecRad_I1S0_t13.push_back(vec_I1S0_t13);
                   	                  
    vecRad_I1S1_t11.push_back(vec_I1S1_t11);
    vecRad_I1S1_t12.push_back(vec_I1S1_t12);
    vecRad_I1S1_t13.push_back(vec_I1S1_t13);
    
  }
  OutFile->cd(); 
  //I0S0
  TGraphErrors* grI0S0_t11 = new TGraphErrors();
  grI0S0_t11->SetName("I0S0_t11");
  int counter = 1; 
  for (auto it :  vecRad_I0S0_t11 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S0_t11->SetPoint(counter-1, dRad*counter, mean);
    grI0S0_t11->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S0_t11->Write("I0S0_t11");

  TGraphErrors* grI0S0_t12 = new TGraphErrors();
  grI0S0_t12->SetName("I0S0_t12");
  counter = 1; 
  for (auto it :  vecRad_I0S0_t12 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S0_t12->SetPoint(counter-1, dRad*counter, mean);
    grI0S0_t12->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S0_t12->Write("I0S0_t12");

  TGraphErrors* grI0S0_t13 = new TGraphErrors();
  grI0S0_t13->SetName("I0S0_t13");
  counter = 1; 
  for (auto it :  vecRad_I0S0_t13 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S0_t13->SetPoint(counter-1, dRad*counter, mean);
    grI0S0_t13->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S0_t13->Write("I0S0_t13");

  //I0S1
  TGraphErrors* grI0S1_t11 = new TGraphErrors();
  grI0S1_t11->SetName("I0S1_t11");
  counter = 1; 
  for (auto it :  vecRad_I0S1_t11 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S1_t11->SetPoint(counter-1, dRad*counter, mean);
    grI0S1_t11->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S1_t11->Write("I0S1_t11");

  TGraphErrors* grI0S1_t12 = new TGraphErrors();
  grI0S1_t12->SetName("I0S1_t12");
  counter = 1; 
  for (auto it :  vecRad_I0S1_t12 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S1_t12->SetPoint(counter-1, dRad*counter, mean);
    grI0S1_t12->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S1_t12->Write("I0S1_t12");

  TGraphErrors* grI0S1_t13 = new TGraphErrors();
  grI0S1_t13->SetName("I0S1_t13");
  counter = 1; 
  for (auto it :  vecRad_I0S1_t13 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI0S1_t13->SetPoint(counter-1, dRad*counter, mean);
    grI0S1_t13->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI0S1_t13->Write("I0S1_t13");

  //I1S0
  TGraphErrors* grI1S0_t11 = new TGraphErrors();
  grI1S0_t11->SetName("I1S0_t11");
  counter = 1; 
  for (auto it :  vecRad_I1S0_t11 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S0_t11->SetPoint(counter-1, dRad*counter, mean);
    grI1S0_t11->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S0_t11->Write("I1S0_t11");

  TGraphErrors* grI1S0_t12 = new TGraphErrors();
  grI1S0_t12->SetName("I1S0_t12");
  counter = 1; 
  for (auto it :  vecRad_I1S0_t12 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S0_t12->SetPoint(counter-1, dRad*counter, mean);
    grI1S0_t12->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S0_t12->Write("I1S0_t12");

  TGraphErrors* grI1S0_t13 = new TGraphErrors();
  grI1S0_t13->SetName("I1S0_t13");
  counter = 1; 
  for (auto it :  vecRad_I1S0_t13 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S0_t13->SetPoint(counter-1, dRad*counter, mean);
    grI1S0_t13->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S0_t13->Write("I1S0_t13");

  //I1S1
  TGraphErrors* grI1S1_t11 = new TGraphErrors();
  grI1S1_t11->SetName("I1S1_t11");
  counter = 1; 
  for (auto it :  vecRad_I1S1_t11 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S1_t11->SetPoint(counter-1, dRad*counter, mean);
    grI1S1_t11->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S1_t11->Write("I1S1_t11");

  TGraphErrors* grI1S1_t12 = new TGraphErrors();
  grI1S1_t12->SetName("I1S1_t12");
  counter = 1; 
  for (auto it :  vecRad_I1S1_t12 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S1_t12->SetPoint(counter-1, dRad*counter, mean);
    grI1S1_t12->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S1_t12->Write("I1S1_t12");

  TGraphErrors* grI1S1_t13 = new TGraphErrors();
  grI1S1_t13->SetName("I1S1_t13");
  counter = 1; 
  for (auto it :  vecRad_I1S1_t13 ) {
    double mean = (it.front() + it.back())/2.;
    double err = it.front() - mean; 
    grI1S1_t13->SetPoint(counter-1, dRad*counter, mean);
    grI1S1_t13->SetPointError(counter-1, 0, err);
    counter++; 
  } 
  grI1S1_t13->Write("I1S1_t13");


  OutFile->Close();
  return; 
} 

void JackknifeCFs() {
  TFile *OutFile = TFile::Open("CFsXi.root", "RECREATE");
  OutFile->cd();
  
  const int nBins = 305;
  double xmin = 0.5;
  double xmax = 305.5;

  std::vector<std::vector<float>> veckStar_tVar(nBins-1);
  
  std::vector<std::vector<float>> veckStar_t11(nBins-1);
  std::vector<std::vector<float>> veckStar_t11_RadVar(nBins-1);
  
  std::vector<std::vector<float>> veckStar_t12(nBins-1);
  std::vector<std::vector<float>> veckStar_t12_RadVar(nBins-1);

  std::vector<std::vector<float>> veckStar_t13(nBins-1);
  std::vector<std::vector<float>> veckStar_t13_RadVar(nBins-1);
  
  TidyCats* tidy = new TidyCats();
  std::vector<float> kStar;
  int total = 71;
  auto start = std::chrono::system_clock::now();
  for (int iTime = 0; iTime < total; ++iTime) {
    auto end = std::chrono::system_clock::now();
    if ( (iTime > 0 && iTime < 24) || (iTime > 48) ) {
      continue; 
    }
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "\r Processing progress: "
	      << TString::Format("%.1f %%", (double)iTime / total * 100.f).Data()
	      << " elapsed time: " << elapsed_seconds.count() / 60. << std::endl 
	      << std::flush;

    double QCDTime = iTime; 
    TidyCats::Sources TheSource = TidyCats::sGaussian;
    float ppRadii[3];
    ppRadii[0] = 0.97;
    ppRadii[1] = 1.02;
    ppRadii[2] = 1.07; 

    CATS HALLO;
    CATS HALLO_RadUp;
    CATS HALLO_RadDown; 
    //full calculation 
    tidy->GetCatsProtonXiMinus(&HALLO, nBins, xmin, xmax, TheSource,
			       TidyCats::pHALQCD, iTime);
    HALLO.SetAnaSource(0, ppRadii[1]);
    HALLO.KillTheCat();
    if ( iTime == 24 ) {
      tidy->GetCatsProtonXiMinus(&HALLO_RadDown, nBins, xmin, xmax, TheSource,
				 TidyCats::pHALQCD, iTime);
      HALLO_RadDown.SetAnaSource(0, ppRadii[0]);
      HALLO_RadDown.KillTheCat();
      
      tidy->GetCatsProtonXiMinus(&HALLO_RadUp, nBins, xmin, xmax, TheSource,
				 TidyCats::pHALQCD, iTime);
      HALLO_RadUp.SetAnaSource(0, ppRadii[2]);
      HALLO_RadUp.KillTheCat();
    } 
    for (auto it = 0; it < nBins-1; ++it) {
      if (iTime == 0) {
	kStar.push_back(HALLO.GetMomentum(it)); 
      } 
      if (iTime == 0 || iTime == 24 || iTime == 48) {
	//Default obtained by average 	
	veckStar_tVar[it].push_back(HALLO.GetCorrFun(it));
	/*
	  if (iTime == 0) {
	  veckStar_t11_RadVar[it].push_back(HALLO.GetCorrFun(it));
	  veckStar_t11_RadVar[it].push_back(HALLO_RadUp.GetCorrFun(it));
	  veckStar_t11_RadVar[it].push_back(HALLO_RadDown.GetCorrFun(it));
	  } else
	*/
	if (iTime == 24) {
	  veckStar_t12_RadVar[it].push_back(HALLO_RadDown.GetCorrFun(it));
	  veckStar_t12_RadVar[it].push_back(HALLO.GetCorrFun(it));
	  veckStar_t12_RadVar[it].push_back(HALLO_RadUp.GetCorrFun(it));	  
	}
	/*else {
	  veckStar_t13_RadVar[it].push_back(HALLO.GetCorrFun(it));
	  veckStar_t13_RadVar[it].push_back(HALLO_RadUp.GetCorrFun(it));
	  veckStar_t13_RadVar[it].push_back(HALLO_RadDown.GetCorrFun(it));
	  }
	*/
      } else { 
	if (iTime < 24) { 
	  veckStar_t11[it].push_back(HALLO.GetCorrFun(it));
	} else if (iTime < 48) {
	  veckStar_t12[it].push_back(HALLO.GetCorrFun(it));
	} else {
	  veckStar_t13[it].push_back(HALLO.GetCorrFun(it)); 
	}
      }
    }
  }
  
  OutFile->cd(); 

  // std::cout << "Doing t = 11\n"; 
  // TGraphErrors* gr_t11 = new TGraphErrors();
  // gr_t11->SetName("HAL_t11");
  
  // int counter = 0; 
  // for (auto it : veckStar_t11) {  
  //   double sum = std::accumulate(it.begin(), it.end(), 0.0);
  //   double sumsq = std::inner_product(it.begin(), it.end(),it.begin(), 0.0); 
  //   int nVar = (int)it.size(); 
  //   double err = std::sqrt(sumsq - (std::pow(sum,2)/double(nVar))); 
  //   double mean = sum/double(nVar);
  //   gr_t11->SetPoint(counter, kStar[counter], mean);
  //   gr_t11->SetPointError(counter, 0, err);
  //   counter++; 
  // } 
  // gr_t11->Write("HAL_t11"); 

  std::cout << "Doing t = 12\n"; 
  TGraphErrors* gr_t12 = new TGraphErrors();
  gr_t12->SetName("HAL_t12");
  
  int counter = 0; 
  for (auto it : veckStar_t12) {
    double sum = std::accumulate(it.begin(), it.end(), 0.0);
    double sumsq = std::inner_product(it.begin(), it.end(),it.begin(), 0.0); 
    int nVar = (int)it.size();
    double err = std::sqrt(sumsq - (std::pow(sum,2)/double(nVar))); 
    double mean = sum/double(nVar);
    gr_t12->SetPoint(counter, kStar[counter], mean);
    gr_t12->SetPointError(counter, 0, err);
    counter++; 
  } 
  gr_t12->Write("HAL_t12"); 

  // std::cout << "Doing t = 13\n"; 
  // TGraphErrors* gr_t13 = new TGraphErrors();
  // gr_t13->SetName("HAL_t13");
  // counter = 0; 
  // for (auto it : veckStar_t13) {
  //   double sum = std::accumulate(it.begin(), it.end(), 0.0);
  //   double sumsq = std::inner_product(it.begin(), it.end(),it.begin(), 0.0); 
  //   int nVar = (int)it.size(); 
  //   double err = std::sqrt(sumsq - (std::pow(sum,2)/double(nVar))); 
  //   double mean = sum/double(nVar);
  //   gr_t13->SetPoint(counter, kStar[counter], mean);
  //   gr_t13->SetPointError(counter, 0, err);
  //   counter++; 
  // } 
  // gr_t13->Write("HAL_t13"); 

  std::cout << "Doing t Vars\n"; 
  TGraphAsymmErrors* gr_tVar = new TGraphAsymmErrors();
  gr_tVar->SetName("HAL_tVar");
  counter = 0; 
  for (auto it : veckStar_tVar) {
    double errhigh = it[0] - it[1];
    double errlow = it[1] - it[2]; 
    double mean = it[1]; 
    gr_tVar->SetPoint(counter, kStar[counter], mean);
    gr_tVar->SetPointError(counter, 0, 0, errlow, errhigh);
    counter++; 
  } 
  gr_tVar->Write("HAL_tVar"); 
  
  std::cout << "Doing HAL Combined \n"; 
  TGraphAsymmErrors* gr_HALSum = new TGraphAsymmErrors();
  gr_HALSum->SetName("HAL_Sum");
  TGraphErrors* gr_HALSymSum = new TGraphErrors();
  gr_HALSymSum->SetName("HAL_SymSum");
  double mean_tVar, kstar_tVar;
  double mean_t12, kstar_t12; 
  for (int iPoint = 0; iPoint < gr_tVar->GetN(); ++iPoint) {
    gr_tVar->GetPoint(iPoint, kstar_tVar, mean_tVar);
    gr_t12->GetPoint(iPoint, kstar_t12, mean_t12); 
    if (std::fabs(mean_tVar - mean_t12)/mean_tVar > 1e-2 ) {
      std::cout << "mean_tVar = " << mean_tVar << " mean_t12: " << mean_t12 << std::endl; 
    } 
    //double err = std::sqrt( std::pow(gr_tVar->GetErrorY(iPoint),2) + std::pow(gr_t12->GetErrorY(iPoint),2) );
    double errlow = std::sqrt( std::pow(gr_tVar->GetErrorYlow(iPoint),2) + std::pow(gr_t12->GetErrorY(iPoint),2) );
    double errhigh = std::sqrt( std::pow(gr_tVar->GetErrorYhigh(iPoint),2) + std::pow(gr_t12->GetErrorY(iPoint),2) );
    gr_HALSum->SetPoint(iPoint, kstar_tVar, mean_tVar);
    gr_HALSum->SetPointError(iPoint, 0, 0, errlow, errhigh);

    double errlowSym = mean_tVar - errlow;
    double errhighSym = mean_tVar + errhigh;
    double meanSym = (errlowSym + errhighSym)/2.;
    double errsym_1 = meanSym - errlowSym;
    double errsym_2 = errhighSym - meanSym;

    if ( std::fabs(errsym_1 - errsym_2) > 1e-3 ) {
      std::cout <<" Error is fishy ... errsym_1 =  " << errsym_1 << " errsym_2 = " << errsym_2 << std::endl; 
    }
    gr_HALSymSum->SetPoint(iPoint, kstar_tVar, meanSym);
    gr_HALSymSum->SetPointError(iPoint, 0, errsym_1); 
    
  } 
  gr_HALSum->Write("HAL_Sum");
  gr_HALSymSum->Write("HAL_SymSum"); 
  
  
  std::cout << "Doing Rad Vars\n"; 
  TGraphAsymmErrors* gr_RadVar = new TGraphAsymmErrors();
  gr_RadVar->SetName("HAL_RadVar");
  counter = 0; 
  for (auto it : veckStar_t12_RadVar) {
    double errlow = it[0] - it[1]; 
    double errhigh = it[1] - it[2]; 
    double mean = it[1];
    gr_RadVar->SetPoint(counter, kStar[counter], mean);
    gr_RadVar->SetPointError(counter, 0, 0, errlow, errhigh);
    counter++; 
  } 
  gr_RadVar->Write("HAL_RadVar"); 

  std::cout << "Doing HAL + Rad Combined \n"; 
  TGraphAsymmErrors* gr_HALRadSum = new TGraphAsymmErrors();
  gr_HALRadSum->SetName("HALRad_Sum");
  double mean_HAL, kstar_HAL; 
  double mean_Rad, kstar_Rad;
  
  for (int iPoint = 0; iPoint < gr_HALSum->GetN(); ++iPoint) {
    gr_HALSum->GetPoint(iPoint, kstar_HAL, mean_HAL); 
    gr_RadVar->GetPoint(iPoint, kstar_Rad, mean_Rad); 
    if (std::fabs(mean_HAL - mean_Rad)/mean_HAL > 1e-2 ) {
      std::cout << "mean_HAL = " << mean_HAL << " mean_Rad: " << mean_Rad << std::endl; 
    } 
    double errlow = std::sqrt( std::pow(gr_HALSum->GetErrorYlow(iPoint),2) +
			       std::pow(gr_RadVar->GetErrorYlow(iPoint),2));
    double errhigh = std::sqrt( std::pow(gr_HALSum->GetErrorYhigh(iPoint),2) +
			       std::pow(gr_RadVar->GetErrorYhigh(iPoint),2));
    gr_HALRadSum->SetPoint(iPoint, kstar_HAL, mean_HAL);
    gr_HALRadSum->SetPointError(iPoint, 0, 0, errlow, errhigh);
  } 
  gr_HALRadSum->Write("HALRad_Sum"); 

  std::cout << "Doing Symmetric HAL + Rad Combined \n";
  TGraphErrors* gr_HALRadSymSum = new TGraphErrors();
  gr_HALRadSymSum->SetName("gr_HALRadSymSum");
  double mean, kstar; 
  for (int iPoint = 0; iPoint < gr_HALRadSum->GetN(); ++iPoint) {
    gr_HALRadSum->GetPoint(iPoint, kstar, mean); 
    double errlow = mean - gr_HALRadSum->GetErrorYlow(iPoint);
    double errhigh = mean + gr_HALRadSum->GetErrorYhigh(iPoint);
    double meanSym = (errlow + errhigh)/2.;
    double errsym_1 = meanSym - errlow;
    double errsym_2 = errhigh - meanSym;

    if ( std::fabs(errsym_1 - errsym_2) > 1e-3 ) {
      std::cout <<" Error is fishy ... errsym_1 =  " << errsym_1 << " errsym_2 = " << errsym_2 << std::endl; 
    }
    gr_HALRadSymSum->SetPoint(iPoint, kstar, meanSym);
    gr_HALRadSymSum->SetPointError(iPoint, 0, errsym_1); 
  }
  gr_HALRadSymSum->Write("gr_HALRadSymSum"); 
  
  OutFile->Close();
  return; 
} 


void PlotCF(double GaussSourceSize = 1.4) {

  double QCDTime = 24;
  //4th argument is the t parameter and can be:
  // 9, 10, 11, 12
  /* 
  double pXimPotParsI0S0[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };
  double pXimPotParsI0S1[8] = { pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1};
  double pXimPotParsI1S0[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0};
  double pXimPotParsI1S1[8] = { pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1};
  */
  double pXimPotParsI0S0[8] = { pXim_HALQCDPaper2020, QCDTime, 0, -1, 1, 0, 0, 0 };
  double pXimPotParsI0S1[8] = { pXim_HALQCDPaper2020, QCDTime, 0, -1, 1, 1, 0, 1};
  double pXimPotParsI1S0[8] = { pXim_HALQCDPaper2020, QCDTime, 1, 1, 1, 0, 0, 0};
  double pXimPotParsI1S1[8] = { pXim_HALQCDPaper2020, QCDTime, 1, 1, 1, 1, 0, 1};
  CATSparameters cPotParsI0S0(CATSparameters::tPotential, 8, true);
  cPotParsI0S0.SetParameters(pXimPotParsI0S0);

  CATSparameters cPotParsI0S1(CATSparameters::tPotential, 8, true);
  cPotParsI0S1.SetParameters(pXimPotParsI0S1);

  CATSparameters cPotParsI1S0(CATSparameters::tPotential, 8, true);
  cPotParsI1S0.SetParameters(pXimPotParsI1S0);

  CATSparameters cPotParsI1S1(CATSparameters::tPotential, 8, true);
  cPotParsI1S1.SetParameters(pXimPotParsI1S1);
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;
  int momBins = 55;
  float kMin = 0.;
  float kMax = 220;
  TFile* outfile = TFile::Open("InletCFs.root", "recreate");
  auto c1 = new TCanvas("c1", "c1");
  TLegend leg(0.7, 0.5, 0.9, 0.9);

  TString I0S0HistTitle = Form("I0S0");
  TH1F* CF_I0S0 = new TH1F(I0S0HistTitle.Data(), I0S0HistTitle.Data(), momBins,
                           kMin, kMax);
  TString I0S1HistTitle = Form("I0S1");
  TH1F* CF_I0S1 = new TH1F(I0S1HistTitle.Data(), I0S1HistTitle.Data(), momBins,
                           kMin, kMax);
  TString I1S0HistTitle = Form("I1S0");
  TH1F* CF_I1S0 = new TH1F(I1S0HistTitle.Data(), I1S0HistTitle.Data(), momBins,
                           kMin, kMax);
  TString I1S1HistTitle = Form("I1S1");
  TH1F* CF_I1S1 = new TH1F(I1S1HistTitle.Data(), I1S1HistTitle.Data(), momBins,
                           kMin, kMax);

  CATS AB_pXim_I0S0;
  CATSparameters cParsI0S0(CATSparameters::tSource, 1, true);
  cParsI0S0.SetParameter(0, GaussSourceSize);
  AB_pXim_I0S0.SetAnaSource(GaussSource, cParsI0S0);
  AB_pXim_I0S0.SetUseAnalyticSource(true);
  AB_pXim_I0S0.SetThetaDependentSource(false);
  AB_pXim_I0S0.SetExcludeFailedBins(false);
  AB_pXim_I0S0.SetMomBins(momBins, kMin, kMax);
  AB_pXim_I0S0.SetNumChannels(1);
  AB_pXim_I0S0.SetNumPW(0, 1);
  AB_pXim_I0S0.SetSpin(0, 0);    //I=0; S=0
  AB_pXim_I0S0.SetChannelWeight(0, 1);
  AB_pXim_I0S0.SetQ1Q2(0);
  AB_pXim_I0S0.SetPdgId(2212, 3122);
  AB_pXim_I0S0.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim_I0S0.SetMaxRad(64);
  AB_pXim_I0S0.SetMaxRho(32);
//  AB_pXim_I0S0.SetEpsilonConv(1e-8);
//  AB_pXim_I0S0.SetEpsilonProp(1e-8);

  AB_pXim_I0S0.SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S0);

  AB_pXim_I0S0.KillTheCat();

  CATS AB_pXim_I0S1;
  CATSparameters cParsI0S1(CATSparameters::tSource, 1, true);
  cParsI0S1.SetParameter(0, GaussSourceSize);
  AB_pXim_I0S1.SetAnaSource(GaussSource, cParsI0S1);
  AB_pXim_I0S1.SetUseAnalyticSource(true);
  AB_pXim_I0S1.SetThetaDependentSource(false);
  AB_pXim_I0S1.SetExcludeFailedBins(false);
  AB_pXim_I0S1.SetMomBins(momBins, kMin, kMax);
  AB_pXim_I0S1.SetNumChannels(1);
  AB_pXim_I0S1.SetNumPW(0, 1);
  AB_pXim_I0S1.SetSpin(0, 1);    //I=0; S=0
  AB_pXim_I0S1.SetChannelWeight(0, 1);
  AB_pXim_I0S1.SetQ1Q2(0);
  AB_pXim_I0S1.SetPdgId(2212, 3122);
  AB_pXim_I0S1.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim_I0S1.SetMaxRad(64);
  AB_pXim_I0S1.SetMaxRho(32);
//  AB_pXim_I0S1.SetEpsilonConv(1e-8);
//  AB_pXim_I0S1.SetEpsilonProp(1e-8);
  AB_pXim_I0S1.SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S1);

  AB_pXim_I0S1.KillTheCat();

  CATS AB_pXim_I1S0;
  CATSparameters cParsI1S0(CATSparameters::tSource, 1, true);
  cParsI1S0.SetParameter(0, GaussSourceSize);
  AB_pXim_I1S0.SetAnaSource(GaussSource, cParsI1S0);
  AB_pXim_I1S0.SetUseAnalyticSource(true);
  AB_pXim_I1S0.SetThetaDependentSource(false);
  AB_pXim_I1S0.SetExcludeFailedBins(false);
  AB_pXim_I1S0.SetMomBins(momBins, kMin, kMax);
  AB_pXim_I1S0.SetNumChannels(1);
  AB_pXim_I1S0.SetNumPW(0, 1);
  AB_pXim_I1S0.SetSpin(0, 0);    //I=0; S=0
  AB_pXim_I1S0.SetChannelWeight(0, 1);
  AB_pXim_I1S0.SetQ1Q2(0);
  AB_pXim_I1S0.SetPdgId(2212, 3122);
  AB_pXim_I1S0.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim_I1S0.SetMaxRad(64);
  AB_pXim_I1S0.SetMaxRho(32);
//  AB_pXim_I1S0.SetEpsilonConv(1e-8);
//  AB_pXim_I1S0.SetEpsilonProp(1e-8);
  AB_pXim_I1S0.SetShortRangePotential(0, 0, fDlmPot, cPotParsI1S0);

  AB_pXim_I1S0.KillTheCat();

  CATS AB_pXim_I1S1;
  CATSparameters cParsI1S1(CATSparameters::tSource, 1, true);
  cParsI1S1.SetParameter(0, GaussSourceSize);
  AB_pXim_I1S1.SetAnaSource(GaussSource, cParsI1S1);
  AB_pXim_I1S1.SetUseAnalyticSource(true);
  AB_pXim_I1S1.SetThetaDependentSource(false);
  AB_pXim_I1S1.SetExcludeFailedBins(false);
  AB_pXim_I1S1.SetMomBins(momBins, kMin, kMax);
  AB_pXim_I1S1.SetNumChannels(1);
  AB_pXim_I1S1.SetNumPW(0, 1);
  AB_pXim_I1S1.SetSpin(0, 1);    //I=0; S=0
  AB_pXim_I1S1.SetChannelWeight(0, 1);
  AB_pXim_I1S1.SetQ1Q2(0);
  AB_pXim_I1S1.SetPdgId(2212, 3122);
  AB_pXim_I1S1.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
  AB_pXim_I1S1.SetMaxRad(64);
  AB_pXim_I1S1.SetMaxRho(32);
//  AB_pXim_I1S1.SetEpsilonConv(1e-8);
//  AB_pXim_I1S1.SetEpsilonProp(1e-8);
  AB_pXim_I1S1.SetShortRangePotential(0, 0, fDlmPot, cPotParsI1S1);

  AB_pXim_I1S1.KillTheCat();

  for (int ikS = 0; ikS < momBins; ++ikS) {
    CF_I0S0->SetBinContent(ikS, AB_pXim_I0S0.GetCorrFun(ikS));
    CF_I0S1->SetBinContent(ikS, AB_pXim_I0S1.GetCorrFun(ikS));
    CF_I1S0->SetBinContent(ikS, AB_pXim_I1S0.GetCorrFun(ikS));
    CF_I1S1->SetBinContent(ikS, AB_pXim_I1S1.GetCorrFun(ikS));
  }
  c1->cd();
  CF_I0S0->GetXaxis()->SetRangeUser(0., 200.);
  CF_I0S0->SetStats(0);
  CF_I0S0->SetLineWidth(1.6);
  CF_I0S0->SetLineColor(fColors[1]);
  CF_I0S0->SetMarkerStyle(24);
  CF_I0S0->SetLineStyle(7);
  CF_I0S0->SetMarkerColor(fColors[1]);
  CF_I0S0->SetMarkerSize(1);
  CF_I0S0->Draw("L");

  outfile->cd();
  CF_I0S0->Write();

  c1->cd();
  CF_I0S1->GetXaxis()->SetRangeUser(0., 200.);
  CF_I0S1->SetStats(0);
  CF_I0S1->SetLineWidth(1.6);
  CF_I0S1->SetLineColor(fColors[2]);
  CF_I0S1->SetMarkerStyle(24);
  CF_I0S1->SetLineStyle(3);
  CF_I0S1->SetMarkerColor(fColors[2]);
  CF_I0S1->SetMarkerSize(1);
  CF_I0S1->Draw("LSame");

  outfile->cd();
  CF_I0S1->Write();

  c1->cd();
  CF_I1S0->GetXaxis()->SetRangeUser(0., 200.);
  CF_I1S0->SetStats(0);
  CF_I1S0->SetLineWidth(1.6);
  CF_I1S0->SetLineColor(fColors[3]);
  CF_I1S0->SetMarkerStyle(24);
  CF_I1S0->SetLineStyle(1);
  CF_I1S0->SetMarkerColor(fColors[3]);
  CF_I1S0->SetMarkerSize(1);
  CF_I1S0->Draw("LSame");

  outfile->cd();
  CF_I1S0->Write();

  c1->cd();
  CF_I1S1->GetXaxis()->SetRangeUser(0., 200.);
  CF_I1S1->SetStats(0);
  CF_I1S1->SetLineWidth(1.6);
  CF_I1S1->SetLineColor(fColors[5]);
  CF_I1S1->SetLineStyle(4);
  CF_I1S1->SetMarkerStyle(24);
  CF_I1S1->SetMarkerColor(fColors[5]);
  CF_I1S1->SetMarkerSize(1);
  CF_I1S1->Draw("LSame");

  outfile->cd();
  CF_I1S1->Write();

  c1->Write();
  outfile->Close();
}

void SetStyle(bool graypalette, bool title) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if (graypalette)
    gStyle->SetPalette(8, 0);
  else
    gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetLabelOffset(0.01, "x");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0.005);
}

int main(int argc, char *argv[]) {
  //TApplication app("TheApp", &argc, argv);
  SetStyle(false, false);

  radCutoff = 0.;

  int nBins = 305;
  double xmin = 0.5;
  double xmax = 305.5;
  TidyCats* tidy = new TidyCats();
  TidyCats::Sources TheSource = TidyCats::sGaussian;
  float ppRadii[3];
  ppRadii[0] = 0.97;
  ppRadii[1] = 1.02;
  ppRadii[2] = 1.07; 

  CATS Hal;
  
  tidy->GetCatsProtonXiMinus(&Hal, nBins, xmin, xmax, TheSource,
                             TidyCats::pHALQCD, 24);

  Hal.SetAnaSource(0, ppRadii[1]);
  Hal.KillTheCat();

  TGraphErrors* halandRad = new TGraphErrors();
  halandRad->SetName("HalAndRad");

  for (auto it = 0; it < nBins - 10; ++it) {
    double kStar = Hal.GetMomentum(it);
    double mean = Hal.GetCorrFun(it);
    halandRad->SetPoint(it, kStar, mean);
  }
  TFile* out = TFile::Open("newpXiSasaki.root", "recreate");
  //TFile* out = TFile::Open("oldpXiHatsuda.root", "recreate"); 
  out->cd();
  halandRad->Write("halandrad");
  out->Close(); 
  
  JackknifePotentials();
  JackknifeCFs();

  



  
  //PlotCF(0.95);
  //PlotPotentials();
  //app.Run(); 
  return 0;
}

