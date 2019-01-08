#include "DLM_Source.h"
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

void PlotPotentials() {
  TidyCats* tidy = new TidyCats();
  CATSInput *CATSinput = new CATSInput();
  std::vector<float> cutOffValues = { 0.01 };
//  std::vector<float> cutOffValues = { 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
//      1.1, 1.3, 1.5, 1.7, 5.0, 15.2 };
//  std::vector<float> cutOffValues = { 0.1, 0.3, 0.5, 0.6, 0.7, 0.9, 1.1, 1.3, 1.8 };
  const int nCutOffs = cutOffValues.size();
  TH1F* CF_I0S0[nCutOffs];
  TH1F* CF_I0S1[nCutOffs];
  TH1F* CF_I1S0[nCutOffs];
  TH1F* CF_I1S1[nCutOffs];

  double QCDTime = 11;
  //4th argument is the t parameter and can be:
  // 9, 10, 11, 12
  double pXimPotParsI0S0[10] =
      { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0, 0 };
  double pXimPotParsI0S1[10] =
      { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0, 1 };
  double pXimPotParsI1S0[10] = { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0 };
  double pXimPotParsI1S1[10] = { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1 };

  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;
  int momBins = 44;
  float kMin = 0.;
  float kMax = 220;
  int iCutOff = 0;
  TFile* outfile = TFile::Open("output.root", "recreate");
  auto c1 = new TCanvas("c1", "c1");
  c1->Divide(2, 2);
  TLegend leg(0.7, 0.5, 0.9, 0.9);
  for (auto it : cutOffValues) {
    TString I0S0HistName = Form("hI0S0_%u", iCutOff);
    TString I0S0HistTitle = Form("I0S0");
    CF_I0S0[iCutOff] = new TH1F(I0S0HistName.Data(), I0S0HistTitle.Data(),
                                momBins, kMin, kMax );
    TString I0S1HistName = Form("hI0S1_%u", iCutOff);
    TString I0S1HistTitle = Form("I0S1");
    CF_I0S1[iCutOff] = new TH1F(I0S1HistName.Data(), I0S1HistTitle.Data(),
                                momBins, kMin, kMax );
    TString I1S0HistName = Form("hI1S0_%u", iCutOff);
    TString I1S0HistTitle = Form("I1S0");
    CF_I1S0[iCutOff] = new TH1F(I1S0HistName.Data(), I1S0HistTitle.Data(),
                                momBins, kMin, kMax );
    TString I1S1HistName = Form("hI1S1_%u", iCutOff);
    TString I1S1HistTitle = Form("I1S1");
    CF_I1S1[iCutOff] = new TH1F(I1S1HistName.Data(), I1S1HistTitle.Data(),
                                momBins, kMin, kMax );
    double GaussSourceSize = 1.4;

    CATS AB_pXim_I0S0;
    double Pars_pXiI0S0[6] = { 0, 0, 0, GaussSourceSize, it, 0.5 };
    AB_pXim_I0S0.SetAnaSource(GaussSource, Pars_pXiI0S0);
    AB_pXim_I0S0.SetUseAnalyticSource(true);
    AB_pXim_I0S0.SetThetaDependentSource(false);
    AB_pXim_I0S0.SetExcludeFailedBins(false);
    AB_pXim_I0S0.SetMomBins(momBins, kMin, kMax);
    AB_pXim_I0S0.SetNumChannels(4);
    AB_pXim_I0S0.SetNumPW(0, 1);
    AB_pXim_I0S0.SetSpin(0, 0);    //I=0; S=0
    AB_pXim_I0S0.SetChannelWeight(0, 1);
    AB_pXim_I0S0.SetQ1Q2(-1);
    AB_pXim_I0S0.SetPdgId(2212, 3122);
    AB_pXim_I0S0.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
    AB_pXim_I0S0.SetMaxRad(64);
    AB_pXim_I0S0.SetMaxRho(32);
    AB_pXim_I0S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S0);

    AB_pXim_I0S0.KillTheCat();

    CATS AB_pXim_I0S1;
    double Pars_pXiI0S1[6] = { 0, 0, 0, GaussSourceSize, it, 0.5 };
    AB_pXim_I0S1.SetAnaSource(GaussSource, Pars_pXiI0S1);
    AB_pXim_I0S1.SetUseAnalyticSource(true);
    AB_pXim_I0S1.SetThetaDependentSource(false);
    AB_pXim_I0S1.SetExcludeFailedBins(false);
    AB_pXim_I0S1.SetMomBins(momBins, kMin, kMax);
    AB_pXim_I0S1.SetNumChannels(4);
    AB_pXim_I0S1.SetNumPW(0, 1);
    AB_pXim_I0S1.SetSpin(0, 1);    //I=0; S=0
    AB_pXim_I0S1.SetChannelWeight(0, 1);
    AB_pXim_I0S1.SetQ1Q2(-1);
    AB_pXim_I0S1.SetPdgId(2212, 3122);
    AB_pXim_I0S1.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
    AB_pXim_I0S1.SetMaxRad(64);
    AB_pXim_I0S1.SetMaxRho(32);
    AB_pXim_I0S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S1);

    AB_pXim_I0S1.KillTheCat();

    CATS AB_pXim_I1S0;
    double Pars_pXiI1S0[6] = { 0, 0, 0, GaussSourceSize, it, 0.5 };
    AB_pXim_I1S0.SetAnaSource(GaussSource, Pars_pXiI1S0);
    AB_pXim_I1S0.SetUseAnalyticSource(true);
    AB_pXim_I1S0.SetThetaDependentSource(false);
    AB_pXim_I1S0.SetExcludeFailedBins(false);
    AB_pXim_I1S0.SetMomBins(momBins, kMin, kMax);
    AB_pXim_I1S0.SetNumChannels(4);
    AB_pXim_I1S0.SetNumPW(0, 1);
    AB_pXim_I1S0.SetSpin(0, 0);    //I=0; S=0
    AB_pXim_I1S0.SetChannelWeight(0, 1);
    AB_pXim_I1S0.SetQ1Q2(-1);
    AB_pXim_I1S0.SetPdgId(2212, 3122);
    AB_pXim_I1S0.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
    AB_pXim_I1S0.SetMaxRad(64);
    AB_pXim_I1S0.SetMaxRho(32);
    AB_pXim_I1S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S0);

    AB_pXim_I1S0.KillTheCat();

    CATS AB_pXim_I1S1;
    double Pars_pXiI1S1[6] = { 0, 0, 0, GaussSourceSize, it, 0.5 };
    AB_pXim_I1S1.SetAnaSource(GaussSource, Pars_pXiI1S1);
    AB_pXim_I1S1.SetUseAnalyticSource(true);
    AB_pXim_I1S1.SetThetaDependentSource(false);
    AB_pXim_I1S1.SetExcludeFailedBins(false);
    AB_pXim_I1S1.SetMomBins(momBins, kMin, kMax);
    AB_pXim_I1S1.SetNumChannels(4);
    AB_pXim_I1S1.SetNumPW(0, 1);
    AB_pXim_I1S1.SetSpin(0, 1);    //I=0; S=0
    AB_pXim_I1S1.SetChannelWeight(0, 1);
    AB_pXim_I1S1.SetQ1Q2(-1);
    AB_pXim_I1S1.SetPdgId(2212, 3122);
    AB_pXim_I1S1.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));
    AB_pXim_I1S1.SetMaxRad(64);
    AB_pXim_I1S1.SetMaxRho(32);
    AB_pXim_I1S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S1);

    AB_pXim_I1S1.KillTheCat();

    for (int ikS = 0; ikS < momBins; ++ikS) {
      CF_I0S0[iCutOff]->SetBinContent(ikS, AB_pXim_I0S0.GetCorrFun(ikS));
      CF_I0S1[iCutOff]->SetBinContent(ikS, AB_pXim_I0S1.GetCorrFun(ikS));
      CF_I1S0[iCutOff]->SetBinContent(ikS, AB_pXim_I1S0.GetCorrFun(ikS));
      CF_I1S1[iCutOff]->SetBinContent(ikS, AB_pXim_I1S1.GetCorrFun(ikS));
    }
    c1->cd(1);
    CF_I0S0[iCutOff]->GetYaxis()->SetRangeUser(0.7, 10);
    CF_I0S0[iCutOff]->SetStats(0);
    CF_I0S0[iCutOff]->SetLineWidth(1);
    CF_I0S0[iCutOff]->SetLineColor(iCutOff + 1);
    CF_I0S0[iCutOff]->SetMarkerStyle(24);
    CF_I0S0[iCutOff]->SetMarkerColor(iCutOff + 1);
    CF_I0S0[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I0S0[iCutOff]->Draw("L");
    } else {
      CF_I0S0[iCutOff]->Draw("LSAME");
    }
    outfile->cd();
    CF_I0S0[iCutOff]->Write();

    c1->cd(2);
    CF_I0S1[iCutOff]->GetYaxis()->SetRangeUser(0.7, 2.7);
    CF_I0S1[iCutOff]->SetStats(0);
    CF_I0S1[iCutOff]->SetLineWidth(1);
    CF_I0S1[iCutOff]->SetLineColor(iCutOff + 1);
    CF_I0S1[iCutOff]->SetMarkerStyle(24);
    CF_I0S1[iCutOff]->SetMarkerColor(iCutOff + 1);
    CF_I0S1[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I0S1[iCutOff]->Draw("L");
    } else {
      CF_I0S1[iCutOff]->Draw("LSAME");
    }

    outfile->cd();
    CF_I0S1[iCutOff]->Write();

    c1->cd(3);
    CF_I1S0[iCutOff]->GetYaxis()->SetRangeUser(0.7, 1.9);
    CF_I1S0[iCutOff]->SetStats(0);
    CF_I1S0[iCutOff]->SetLineWidth(1);
    CF_I1S0[iCutOff]->SetLineColor(iCutOff + 1);
    CF_I1S0[iCutOff]->SetMarkerStyle(24);
    CF_I1S0[iCutOff]->SetMarkerColor(iCutOff + 1);
    CF_I1S0[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I1S0[iCutOff]->Draw("L");
    } else {
      CF_I1S0[iCutOff]->Draw("LSAME");
    }
    outfile->cd();
    CF_I1S0[iCutOff]->Write();

    c1->cd(4);
    CF_I1S0[iCutOff]->GetYaxis()->SetRangeUser(0.7, 2.3);
    CF_I1S1[iCutOff]->SetStats(0);
    CF_I1S1[iCutOff]->SetLineWidth(1);
    CF_I1S1[iCutOff]->SetLineColor(iCutOff + 1);
    CF_I1S1[iCutOff]->SetMarkerStyle(24);
    CF_I1S1[iCutOff]->SetMarkerColor(iCutOff + 1);
    CF_I1S1[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I1S1[iCutOff]->Draw("L");
    } else {
      CF_I1S1[iCutOff]->Draw("LSAME");
    }
    outfile->cd();
    CF_I1S1[iCutOff]->Write();

    leg.AddEntry(CF_I0S1[iCutOff], TString::Format("Cut Off %.2f fm", it),
                 "lp");

    iCutOff++;
  }
  outfile->cd();
  c1->cd(2);
  leg.Draw("same");
  c1->Write();
  c1->SaveAs("Frezze.pdf");
  outfile->Close();
  return;
}

int main(int argc, char *argv[]) {

  PlotPotentials();
  return 0;
}
