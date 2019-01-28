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
  histo->GetXaxis()->SetLabelSize(0.07);
  histo->GetXaxis()->SetTitleSize(0.07);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(0.95);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.07);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(0.95);
  histo->SetMarkerSize(1.2);
  histo->SetLineWidth(3);
}

void SetStyleGraph(TGraph *histo) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.0);
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

  TLegend leg(0.65, 0.65, 0.95, 0.95);
  leg.SetLineColor(kWhite);
  leg.SetFillStyle(4000);
  leg.SetTextSize(0.06);

  TLegend leg2(0.55, 0.65, 0.95, 0.95);
  leg2.SetLineColor(kWhite);
  leg2.SetFillStyle(4000);
  leg2.SetTextSize(0.06);

  int nRadBins = 600;
  double dRad = 0.005;
  double radMax = dRad * nRadBins;
  for (int iTime = 0; iTime < 3; ++iTime) {
    double QCDTime = 11 + iTime;

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
      double pXimPotParsI0S0[11] = { iRad * dRad, 0, pXim_HALQCD1, QCDTime, 0,
          -1, 1, 0, 0, 0, radCutoff };
      double pXimPotParsI0S1[11] = { iRad * dRad, 0, pXim_HALQCD1, QCDTime, 0,
          -1, 1, 1, 0, 1, radCutoff };
      double pXimPotParsI1S0[11] = { iRad * dRad, 0, pXim_HALQCD1, QCDTime, 1,
          1, 1, 0, 0, 0, radCutoff };
      double pXimPotParsI1S1[11] = { iRad * dRad, 0, pXim_HALQCD1, QCDTime, 1,
          1, 1, 1, 0, 1, radCutoff };
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
  grI0S0->GetXaxis()->SetTitle("r (fm)");
  grI0S0->GetYaxis()->SetRangeUser(-70, 50);
  grI0S0->GetXaxis()->SetRangeUser(0, 3);
  grI0S0->GetYaxis()->SetTitle("V(r) (MeV)");

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
  c2->SetLeftMargin(0.11);
  c2->SetRightMargin(0.01);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.12);
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


  TLine lineOne = TLine(0,0,3.,0);
  lineOne.Draw("SAME");

  TLine lineTwo = TLine(0,1.,200,1.);
  lineOne.Draw("SAME");
  TPad *inset_pad = new TPad("inset", "inset", 0.47, 0.13, 0.97, 0.59);
  inset_pad->SetRightMargin(0.05);
  inset_pad->SetFillStyle(4000);
  inset_pad->Draw();
  inset_pad->cd();

  SetStyleHisto(CF_I0S0);
  SetStyleHisto(CF_I0S1);
  SetStyleHisto(CF_I1S0);
  SetStyleHisto(CF_I1S1);
  CF_I0S0->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_I0S0->GetYaxis()->SetRangeUser(0.8,3.5);
  CF_I0S0->SetNdivisions(505);
  CF_I0S0->Draw("l");
  lineTwo.Draw("same");
  CF_I0S1->Draw("lsame");
  CF_I1S0->Draw("lsame");
  CF_I1S1->Draw("lsame");
  TLatex Numbering;
  Numbering.SetTextSize(gStyle->GetTextSize() * 2.0);
  Numbering.SetNDC(kTRUE);
  Numbering.DrawLatex( 0.5 , 0.7, "#it{r}_{Gauss} = 1.4 fm");
  leg2.AddEntry(CF_I0S0," ", "l");
  leg2.AddEntry(CF_I0S1," ", "l");
  leg2.AddEntry(CF_I1S0," ", "l");
  leg2.AddEntry(CF_I1S1," ", "l");

  c2->cd();
  leg2.Draw("SAME");
  leg.Draw("same");
  OutFile->cd();
  c2->Write();
  c2->SaveAs("c2.pdf");
  OutFile->Close();
}

void PlotCF(double GaussSourceSize = 1.4) {

  double QCDTime = 11;
  //4th argument is the t parameter and can be:
  // 9, 10, 11, 12
  double pXimPotParsI0S0[11] =
      { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0, radCutoff };
  double pXimPotParsI0S1[11] =
      { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1, radCutoff };
  double pXimPotParsI1S0[11] = { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0, radCutoff };
  double pXimPotParsI1S1[11] = { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1, radCutoff };

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
  double Pars_pXiI0S0[4] =
      { 0, 0, 0, GaussSourceSize};
  AB_pXim_I0S0.SetAnaSource(GaussSource, Pars_pXiI0S0);
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
  AB_pXim_I0S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S0);

  AB_pXim_I0S0.KillTheCat();

  CATS AB_pXim_I0S1;
  double Pars_pXiI0S1[6] =
      { 0, 0, 0, GaussSourceSize, GaussSourceSize * 2, 0.5 };
  AB_pXim_I0S1.SetAnaSource(GaussSource, Pars_pXiI0S1);
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
  AB_pXim_I0S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S1);

  AB_pXim_I0S1.KillTheCat();

  CATS AB_pXim_I1S0;
  double Pars_pXiI1S0[6] =
      { 0, 0, 0, GaussSourceSize, GaussSourceSize * 2, 0.5 };
  AB_pXim_I1S0.SetAnaSource(GaussSource, Pars_pXiI1S0);
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
  AB_pXim_I1S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S0);

  AB_pXim_I1S0.KillTheCat();

  CATS AB_pXim_I1S1;
  double Pars_pXiI1S1[6] =
      { 0, 0, 0, GaussSourceSize, GaussSourceSize * 2, 0.5 };
  AB_pXim_I1S1.SetAnaSource(GaussSource, Pars_pXiI1S1);
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
  AB_pXim_I1S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S1);

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
  SetStyle(false, false);
  radCutoff =  0.;
  PlotCF();
  PlotPotentials();
  return 0;
}

