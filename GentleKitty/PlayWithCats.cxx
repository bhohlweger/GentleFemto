/*
 * PlayWithCats.cxx
 *
 *  Created on: Jan 14, 2019
 *      Author: schmollweger
 */
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
#include "PlayWithCats.h"

PlayWithCats::PlayWithCats()
    : fCFHistData(nullptr),
      fCFHistErrorFit(nullptr),
      fCFHistDefaultData(nullptr),
      fCFHistDefaultFit(nullptr),
      fCFCoulombOnly(nullptr),
      fQCDTime(11),
      fMass_p(TDatabasePDG::Instance()->GetParticle(2212)->Mass() * 1000),
      fMass_Xim(TDatabasePDG::Instance()->GetParticle(3312)->Mass() * 1000),
      fMomBins(80),
      fkMin(0.),
      fkMax(400.),
      fRadius(1.421) {
  // TODO Auto-generated constructor stub
  fOutFile = TFile::Open("output.root", "recreate");
}

PlayWithCats::~PlayWithCats() {
  // TODO Auto-generated destructor stub
}

void PlayWithCats::PlotPotentials() {
  std::vector<float> cutOffValues = { 0., 0.1, 0.3, 0.5, 0.55, 0.6, 0.65, 0.7,
      0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0, 5.0, 100. };
  const int nCutOffs = cutOffValues.size();
  TH1F* CF_I0S0[nCutOffs];
  TH1F* CF_I0S1[nCutOffs];
  TH1F* CF_I1S0[nCutOffs];
  TH1F* CF_I1S1[nCutOffs];
  TLine one = TLine(0., 1., 400., 1.);
//4th argument is the t parameter and can be:
// 9, 10, 11, 12

  int RadBins = 1000;
  float rMin = 0.;
  float rMax = 10.;
  float dRad = (rMax - rMin) / (float) RadBins;
  int iCutOff = 0;
  auto c1 = new TCanvas("c1", "c1");
  auto c2 = new TCanvas("c2", "c2");
  auto c3 = new TCanvas("c3", "c3");
  auto c4 = new TCanvas("c4", "c4");
  TLegend leg(0.6, 0.5, 0.9, 0.9);
  TList *ListI0S0 = new TList();
  ListI0S0->SetOwner();
  ListI0S0->SetName("I0S0");

  TList *ListI0S1 = new TList();
  ListI0S1->SetOwner();
  ListI0S1->SetName("I0S1");

  TList *ListI1S0 = new TList();
  ListI1S0->SetOwner();
  ListI1S0->SetName("I1S0");

  TList *ListI1S1 = new TList();
  ListI1S1->SetOwner();
  ListI1S1->SetName("I1S1");

  for (auto it : cutOffValues) {
    int COLORMEBLIND = 51 + 3 * iCutOff;
    double PotParsI0S0[9] = { pXim_HALQCD1, fQCDTime, 0, -1, 1, 0, 0, 0, it };
    double PotParsI0S1[9] = { pXim_HALQCD1, fQCDTime, 0, -1, 1, 1, 0, 1, it };
    double PotParsI1S0[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 0, 0, 0, it };
    double PotParsI1S1[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 1, 0, 1, it };

    CATSparameters pXimPotParsI0S0(CATSparameters::tPotential, 9, true);
    pXimPotParsI0S0.SetParameters(PotParsI0S0);
    CATSparameters pXimPotParsI0S1(CATSparameters::tPotential, 9, true);
    pXimPotParsI0S1.SetParameters(PotParsI0S1);
    CATSparameters pXimPotParsI1S0(CATSparameters::tPotential, 9, true);
    pXimPotParsI1S0.SetParameters(PotParsI1S0);
    CATSparameters pXimPotParsI1S1(CATSparameters::tPotential, 9, true);
    pXimPotParsI1S1.SetParameters(PotParsI1S1);

    TString I0S0HistName = Form("hI0S0_%.2f", it);
    TString I0S0HistTitle = Form("I0S0");
    CF_I0S0[iCutOff] = new TH1F(I0S0HistName.Data(), I0S0HistTitle.Data(),
                                fMomBins, fkMin, fkMax);

    TString I0S1HistName = Form("hI0S1_%.2f", it);
    TString I0S1HistTitle = Form("I0S1");
    CF_I0S1[iCutOff] = new TH1F(I0S1HistName.Data(), I0S1HistTitle.Data(),
                                fMomBins, fkMin, fkMax);

    TString I1S0HistName = Form("hI1S0_%.2f", it);
    TString I1S0HistTitle = Form("I1S0");
    CF_I1S0[iCutOff] = new TH1F(I1S0HistName.Data(), I1S0HistTitle.Data(),
                                fMomBins, fkMin, fkMax);

    TString I1S1HistName = Form("hI1S1_%.2f", it);
    TString I1S1HistTitle = Form("I1S1");
    CF_I1S1[iCutOff] = new TH1F(I1S1HistName.Data(), I1S1HistTitle.Data(),
                                fMomBins, fkMin, fkMax);

    CATS AB_pXim_I0S0;
    CATSparameters cParsI0S0(CATSparameters::tSource, 1, true);
    cParsI0S0.SetParameter(0, fRadius);
    AB_pXim_I0S0.SetAnaSource(GaussSource, cParsI0S0);
    AB_pXim_I0S0.SetUseAnalyticSource(true);
    AB_pXim_I0S0.SetThetaDependentSource(false);
    AB_pXim_I0S0.SetExcludeFailedBins(false);
    AB_pXim_I0S0.SetMomBins(fMomBins, fkMin, fkMax);
    AB_pXim_I0S0.SetNumChannels(1);
    AB_pXim_I0S0.SetNumPW(0, 1);
    AB_pXim_I0S0.SetSpin(0, 0);    //I=0; S=0
    AB_pXim_I0S0.SetChannelWeight(0, 1);
    AB_pXim_I0S0.SetQ1Q2(0);
    AB_pXim_I0S0.SetPdgId(2212, 3122);
    AB_pXim_I0S0.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
    AB_pXim_I0S0.SetMaxRad(64);
    AB_pXim_I0S0.SetMaxRho(32);
    AB_pXim_I0S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S0);

    AB_pXim_I0S0.KillTheCat();

    CATS AB_pXim_I0S1;
    CATSparameters cParsI0S1(CATSparameters::tSource, 1, true);
    cParsI0S1.SetParameter(0, fRadius);
    AB_pXim_I0S1.SetAnaSource(GaussSource, cParsI0S1);
    AB_pXim_I0S1.SetUseAnalyticSource(true);
    AB_pXim_I0S1.SetThetaDependentSource(false);
    AB_pXim_I0S1.SetExcludeFailedBins(false);
    AB_pXim_I0S1.SetMomBins(fMomBins, fkMin, fkMax);
    AB_pXim_I0S1.SetNumChannels(1);
    AB_pXim_I0S1.SetNumPW(0, 1);
    AB_pXim_I0S1.SetSpin(0, 1);    //I=0; S=0
    AB_pXim_I0S1.SetChannelWeight(0, 1);
    AB_pXim_I0S1.SetQ1Q2(0);
    AB_pXim_I0S1.SetPdgId(2212, 3122);
    AB_pXim_I0S1.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
    AB_pXim_I0S1.SetMaxRad(64);
    AB_pXim_I0S1.SetMaxRho(32);
    AB_pXim_I0S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S1);

    AB_pXim_I0S1.KillTheCat();

    CATS AB_pXim_I1S0;
    CATSparameters cParsI1S0(CATSparameters::tSource, 1, true);
    cParsI1S0.SetParameter(0, fRadius);
    AB_pXim_I1S0.SetAnaSource(GaussSource, cParsI1S0);
    AB_pXim_I1S0.SetUseAnalyticSource(true);
    AB_pXim_I1S0.SetThetaDependentSource(false);
    AB_pXim_I1S0.SetExcludeFailedBins(false);
    AB_pXim_I1S0.SetMomBins(fMomBins, fkMin, fkMax);
    AB_pXim_I1S0.SetNumChannels(1);
    AB_pXim_I1S0.SetNumPW(0, 1);
    AB_pXim_I1S0.SetSpin(0, 0);    //I=0; S=0
    AB_pXim_I1S0.SetChannelWeight(0, 1);
    AB_pXim_I1S0.SetQ1Q2(0);
    AB_pXim_I1S0.SetPdgId(2212, 3122);
    AB_pXim_I1S0.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
    AB_pXim_I1S0.SetMaxRad(64);
    AB_pXim_I1S0.SetMaxRho(32);
    AB_pXim_I1S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S0);

    AB_pXim_I1S0.KillTheCat();

    CATS AB_pXim_I1S1;
    CATSparameters cParsI1S1(CATSparameters::tSource, 1, true);
    cParsI1S1.SetParameter(0, fRadius);
    AB_pXim_I1S1.SetAnaSource(GaussSource, cParsI1S1);
    AB_pXim_I1S1.SetUseAnalyticSource(true);
    AB_pXim_I1S1.SetThetaDependentSource(false);
    AB_pXim_I1S1.SetExcludeFailedBins(false);
    AB_pXim_I1S1.SetMomBins(fMomBins, fkMin, fkMax);
    AB_pXim_I1S1.SetNumChannels(1);
    AB_pXim_I1S1.SetNumPW(0, 1);
    AB_pXim_I1S1.SetSpin(0, 1);    //I=0; S=0
    AB_pXim_I1S1.SetChannelWeight(0, 1);
    AB_pXim_I1S1.SetQ1Q2(0);
    AB_pXim_I1S1.SetPdgId(2212, 3122);
    AB_pXim_I1S1.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
    AB_pXim_I1S1.SetMaxRad(64);
    AB_pXim_I1S1.SetMaxRho(32);
    AB_pXim_I1S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S1);

    AB_pXim_I1S1.KillTheCat();

    for (int ikS = 0; ikS < fMomBins; ++ikS) {
      CF_I0S0[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I0S0.GetCorrFun(ikS));
      CF_I0S1[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I0S1.GetCorrFun(ikS));
      CF_I1S0[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I1S0.GetCorrFun(ikS));
      CF_I1S1[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I1S1.GetCorrFun(ikS));
    }
    c1->cd();
//    CF_I0S0[iCutOff]->GetYaxis()->SetRangeUser(0.8, 3.5);
    CF_I0S0[iCutOff]->GetXaxis()->SetRangeUser(0., 400.);
    CF_I0S0[iCutOff]->SetStats(0);
    CF_I0S0[iCutOff]->SetLineWidth(1);
    CF_I0S0[iCutOff]->SetLineColor(COLORMEBLIND);
    CF_I0S0[iCutOff]->SetMarkerStyle(24);
    CF_I0S0[iCutOff]->SetMarkerColor(COLORMEBLIND);
    CF_I0S0[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I0S0[iCutOff]->Draw("L");
      one.Draw("SAME");
    } else {
      CF_I0S0[iCutOff]->Draw("LSAME");
    }
    fOutFile->cd();
    ListI0S0->Add(CF_I0S0[iCutOff]);

    c2->cd();
    CF_I0S1[iCutOff]->GetYaxis()->SetRangeUser(0.95, 1.8);
    CF_I0S1[iCutOff]->GetXaxis()->SetRangeUser(0., 400.);
    CF_I0S1[iCutOff]->SetStats(0);
    CF_I0S1[iCutOff]->SetLineWidth(1);
    CF_I0S1[iCutOff]->SetLineColor(COLORMEBLIND);
    CF_I0S1[iCutOff]->SetMarkerStyle(24);
    CF_I0S1[iCutOff]->SetMarkerColor(COLORMEBLIND);
    CF_I0S1[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I0S1[iCutOff]->Draw("L");
      one.Draw("SAME");
    } else {
      CF_I0S1[iCutOff]->Draw("LSAME");
    }

    fOutFile->cd();
    ListI0S1->Add(CF_I0S1[iCutOff]);

    c3->cd();
    CF_I1S0[iCutOff]->GetYaxis()->SetRangeUser(0.95, 1.4);
    CF_I1S0[iCutOff]->GetXaxis()->SetRangeUser(0., 400.);
    CF_I1S0[iCutOff]->SetStats(0);
    CF_I1S0[iCutOff]->SetLineWidth(1);
    CF_I1S0[iCutOff]->SetLineColor(COLORMEBLIND);
    CF_I1S0[iCutOff]->SetMarkerStyle(24);
    CF_I1S0[iCutOff]->SetMarkerColor(COLORMEBLIND);
    CF_I1S0[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I1S0[iCutOff]->Draw("L");
      one.Draw("SAME");
    } else {
      CF_I1S0[iCutOff]->Draw("LSAME");
    }
    fOutFile->cd();
    ListI1S0->Add(CF_I1S0[iCutOff]);

    c4->cd();
    CF_I1S1[iCutOff]->GetYaxis()->SetRangeUser(0.95, 1.5);
    CF_I1S1[iCutOff]->GetXaxis()->SetRangeUser(0., 400.);
    CF_I1S1[iCutOff]->SetStats(0);
    CF_I1S1[iCutOff]->SetLineWidth(1);
    CF_I1S1[iCutOff]->SetLineColor(COLORMEBLIND);
    CF_I1S1[iCutOff]->SetMarkerStyle(24);
    CF_I1S1[iCutOff]->SetMarkerColor(COLORMEBLIND);
    CF_I1S1[iCutOff]->SetMarkerSize(1);
    if (iCutOff == 0) {
      CF_I1S1[iCutOff]->Draw("L");
      one.Draw("SAME");
    } else {
      CF_I1S1[iCutOff]->Draw("LSAME");
    }
    fOutFile->cd();
    ListI1S1->Add(CF_I1S1[iCutOff]);

    leg.AddEntry(CF_I0S1[iCutOff], TString::Format("Cut off %.2f", it), "lp");

    iCutOff++;
  }
  fOutFile->cd();
  c1->cd();
  leg.Draw("same");
  c2->cd();
  leg.Draw("same");
  c3->cd();
  leg.Draw("same");
  c4->cd();
  leg.Draw("same");
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c1->SaveAs("I0S0.pdf");
  c2->SaveAs("I0S1.pdf");
  c3->SaveAs("I1S0.pdf");
  c4->SaveAs("I1S1.pdf");
  fOutFile->cd();
  ListI0S0->Write("I0S0", 1);
  ListI0S1->Write("I0S1", 1);
  ListI1S0->Write("I1S0", 1);
  ListI1S1->Write("I1S1", 1);
  return;
}

void PlayWithCats::PlotPotentialSum() {
  std::vector<float> cutOffValues = { 0., 0.1, 0.3, 0.5, 0.55, 0.6, 0.65, 0.7,
      0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 2.0 };
  const int nCutOffs = cutOffValues.size();
  TH1F* CF_Sum[nCutOffs];
  TLine one = TLine(fkMin, 1., 150, 1.);
//4th argument is the t parameter and can be:
// 9, 10, 11, 12

  int RadBins = 1000;
  float rMin = 0.;
  float rMax = 10.;
  float dRad = (rMax - rMin) / (float) RadBins;
  int iCutOff = 0;
  auto c1 = new TCanvas("cSumData", "cSumData");
  c1->cd();
  fCFHistDefaultData->Draw();
  one.Draw("SAME");

  auto c2 = new TCanvas("cSumFit", "cSumFit");
  c2->cd();
  fCFHistDefaultFit->Draw();
  one.Draw("SAME");

  TList *ListSum = new TList();
  ListSum->SetOwner();
  ListSum->SetName("Sum");

  TLegend legData(0.6, 0.5, 0.9, 0.9);
  TLegend legFit(0.6, 0.5, 0.9, 0.9);

  for (auto it : cutOffValues) {
    int COLORMEBLIND = 51 + 3 * iCutOff;
    double pXimPotParsI0S0[9] =
        { pXim_HALQCD1, fQCDTime, 0, -1, 1, 0, 0, 0, it };
    double pXimPotParsI0S1[9] =
        { pXim_HALQCD1, fQCDTime, 0, -1, 1, 1, 0, 1, it };
    double pXimPotParsI1S0[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 0, 0, 0, it };
    double pXimPotParsI1S1[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 1, 0, 1, it };

    CATSparameters cPotParsI0S0(CATSparameters::tPotential, 9, true);
    cPotParsI0S0.SetParameters(pXimPotParsI0S0);

    CATSparameters cPotParsI0S1(CATSparameters::tPotential, 9, true);
    cPotParsI0S1.SetParameters(pXimPotParsI0S1);

    CATSparameters cPotParsI1S0(CATSparameters::tPotential, 9, true);
    cPotParsI1S0.SetParameters(pXimPotParsI1S0);

    CATSparameters cPotParsI1S1(CATSparameters::tPotential, 9, true);
    cPotParsI1S1.SetParameters(pXimPotParsI1S1);

    TString SumHistName = Form("Sum_%.2f", it);
    TString SumHistTitle = Form("Sum");
    CF_Sum[iCutOff] = new TH1F(SumHistName.Data(), SumHistTitle.Data(),
                               fMomBins, fkMin, fkMax);

    CATS AB_pXimSum;
    CATSparameters cPars(CATSparameters::tSource, 1, true);
    cPars.SetParameter(0, fRadius);
    AB_pXimSum.SetAnaSource(GaussSource, cPars);
    AB_pXimSum.SetUseAnalyticSource(true);
    AB_pXimSum.SetThetaDependentSource(false);
    AB_pXimSum.SetExcludeFailedBins(false);
    AB_pXimSum.SetMomBins(fMomBins, fkMin, fkMax);

    AB_pXimSum.SetNumChannels(4);
    AB_pXimSum.SetNumPW(0, 1);
    AB_pXimSum.SetNumPW(1, 1);
    AB_pXimSum.SetNumPW(2, 1);
    AB_pXimSum.SetNumPW(3, 1);
    AB_pXimSum.SetSpin(0, 0);    //I=0; S=0
    AB_pXimSum.SetSpin(1, 1);    //I=0; S=1
    AB_pXimSum.SetSpin(2, 0);    //I=1; S=0
    AB_pXimSum.SetSpin(3, 1);    //I=1; S=1
    AB_pXimSum.SetChannelWeight(0, 1. / 8.);
    AB_pXimSum.SetChannelWeight(1, 3. / 8.);
    AB_pXimSum.SetChannelWeight(2, 1. / 8.);
    AB_pXimSum.SetChannelWeight(3, 3. / 8.);

    AB_pXimSum.SetQ1Q2(-1);
    AB_pXimSum.SetPdgId(2212, 3122);
    AB_pXimSum.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
    AB_pXimSum.SetMaxRad(64);
    AB_pXimSum.SetMaxRho(32);
    AB_pXimSum.SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S0);
    AB_pXimSum.SetShortRangePotential(1, 0, fDlmPot, cPotParsI0S1);
    AB_pXimSum.SetShortRangePotential(2, 0, fDlmPot, cPotParsI1S0);
    AB_pXimSum.SetShortRangePotential(3, 0, fDlmPot, cPotParsI1S1);

    AB_pXimSum.KillTheCat();

    double lampXi = 0.513;
    for (int ikS = 0; ikS < fMomBins; ++ikS) {
      CF_Sum[iCutOff]->SetBinContent(
          ikS + 1, lampXi * (AB_pXimSum.GetCorrFun(ikS) - 1.) + 1.);
    }
    CF_Sum[iCutOff]->GetYaxis()->SetRangeUser(
        0.8, CF_Sum[iCutOff]->GetMaximum() * 1.5);
    CF_Sum[iCutOff]->GetXaxis()->SetRangeUser(fkMin, fkMax);
    CF_Sum[iCutOff]->SetStats(0);
    CF_Sum[iCutOff]->SetLineWidth(1);
    CF_Sum[iCutOff]->SetLineColor(COLORMEBLIND);
    CF_Sum[iCutOff]->SetMarkerStyle(24);
    CF_Sum[iCutOff]->SetMarkerColor(COLORMEBLIND);
    CF_Sum[iCutOff]->SetMarkerSize(1);
    c1->cd();
    CF_Sum[iCutOff]->Draw("LSAME");
    c2->cd();
    CF_Sum[iCutOff]->Draw("LSAME");
    fOutFile->cd();
    ListSum->Add(CF_Sum[iCutOff]);

    double Chi2_pXimDataErr = 0;
    double Chi2_pXimFitErr = 0;
    int maxkStarBin = fCFHistDefaultData->FindBin(200.);
    for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {
      //double dataX;
      double theoryY;

      theoryY = CF_Sum[iCutOff]->GetBinContent(uBin + 1);

      double CkData = fCFHistDefaultData->GetBinContent(uBin + 1);
      double dataErr = fCFHistDefaultData->GetBinError(uBin + 1);
      Chi2_pXimDataErr += (CkData - theoryY) * (CkData - theoryY)
          / (dataErr * dataErr);

      double CkFit = fCFHistDefaultFit->GetBinContent(uBin + 1);
      double fitErr = fCFHistDefaultFit->GetBinError(uBin + 1);
      Chi2_pXimFitErr += (CkFit - theoryY) * (CkFit - theoryY)
          / (fitErr * fitErr);
    }
    double pvalXiData = TMath::Prob(Chi2_pXimDataErr, round(maxkStarBin - 1));
    double nSigmaXiData = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiData);

    double pvalXiFit = TMath::Prob(Chi2_pXimFitErr, round(maxkStarBin - 1));
    double nSigmaXiFit = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiFit);

    legData.AddEntry(
        CF_Sum[iCutOff],
        TString::Format("Cut off %.2f (n#sigma = %.2f)", it, nSigmaXiData),
        "lp");
    legFit.AddEntry(
        CF_Sum[iCutOff],
        TString::Format("Cut off %.2f (n#sigma = %.2f)", it, nSigmaXiFit),
        "lp");

    iCutOff++;
  }
  legData.AddEntry(fCFCoulombOnly, TString::Format("Coulomb Only"), "l");
  legFit.AddEntry(fCFCoulombOnly, TString::Format("Coulomb Only"), "l");

  fOutFile->cd();
  c1->cd();
  legData.Draw("same");
  fCFCoulombOnly->Draw("lSame");
  fCFHistDefaultData->Draw("Same");
  c1->Write();
  c1->SaveAs("CkSumvsDataErr.pdf");

  c2->cd();
  legFit.Draw("same");
  fCFCoulombOnly->Draw("lSame");
  fCFHistDefaultFit->Draw("Same");
  c2->Write();
  c2->SaveAs("CkSumvsFitErr.pdf");

  ListSum->Write("SumWithCoulomb", 1);
  return;
}

void PlayWithCats::ExtractUncertaintyData(const char* inPath) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->ReadCorrelationFile(inPath);
  CATSinput->ObtainCFs(5, 240, 340);

  TString HistpXiName = "hCk_ReweightedpXiMeV_0";
  fCFHistData = CATSinput->GetCF("pXi", HistpXiName.Data());
  fkMin = fCFHistData->GetXaxis()->GetXmin();
}

void PlayWithCats::ExtractUncertaintyFit(const char* inFit) {
  double dummy;
  double yUpper;
  double yLower;
  TFile* catsFile = TFile::Open(inFit, "read");
  TGraph* upper = (TGraph*) catsFile->Get("pXimGraphUpperLim");
  TGraph* lower = (TGraph*) catsFile->Get("pXimGraphLowerLim");
  fMomBins = upper->GetN();
  fkMax = upper->GetXaxis()->GetXmax();

  fCFHistErrorFit = new TH1F("pXiFitErr", "pXiFitErr", fMomBins, fkMin, fkMax);
  for (int iBin = 0; iBin < fMomBins; ++iBin) {
    fCFHistErrorFit->SetBinContent(iBin + 1, 1);
    upper->GetPoint(iBin, dummy, yUpper);
    lower->GetPoint(iBin, dummy, yLower);
    fCFHistErrorFit->SetBinError(iBin + 1, yUpper - yLower);
  }
}

void PlayWithCats::GenerateDefault() {
  double pXimPotParsI0S0[9] = { pXim_HALQCD1, fQCDTime, 0, -1, 1, 0, 0, 0, 0. };
  double pXimPotParsI0S1[9] = { pXim_HALQCD1, fQCDTime, 0, -1, 1, 1, 0, 1, 0. };
  double pXimPotParsI1S0[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 0, 0, 0, 0. };
  double pXimPotParsI1S1[9] = { pXim_HALQCD1, fQCDTime, 1, 1, 1, 1, 0, 1, 0. };

  TString DataHistName = Form("DefaultErrData");
  TString DataHistTitle = Form("DefaultErrData");
  fCFHistDefaultData = new TH1F(DataHistName.Data(), DataHistTitle.Data(),
                                fMomBins, fkMin, fkMax);

  TString FitHistName = Form("DefaultErrFit");
  TString FitHistTitle = Form("DefaultErrFit");
  fCFHistDefaultFit = new TH1F(FitHistName.Data(), FitHistTitle.Data(),
                               fMomBins, fkMin, fkMax);

  CATS AB_pXimSum;
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, fRadius);

  AB_pXimSum.SetAnaSource(GaussSource, cPars);
  AB_pXimSum.SetUseAnalyticSource(true);
  AB_pXimSum.SetThetaDependentSource(false);
  AB_pXimSum.SetExcludeFailedBins(false);
  AB_pXimSum.SetMomBins(fMomBins, fkMin, fkMax);

  AB_pXimSum.SetNumChannels(4);
  AB_pXimSum.SetNumPW(0, 1);
  AB_pXimSum.SetNumPW(1, 1);
  AB_pXimSum.SetNumPW(2, 1);
  AB_pXimSum.SetNumPW(3, 1);
  AB_pXimSum.SetSpin(0, 0);    //I=0; S=0
  AB_pXimSum.SetSpin(1, 1);    //I=0; S=1
  AB_pXimSum.SetSpin(2, 0);    //I=1; S=0
  AB_pXimSum.SetSpin(3, 1);    //I=1; S=1
  AB_pXimSum.SetChannelWeight(0, 1. / 8.);
  AB_pXimSum.SetChannelWeight(1, 3. / 8.);
  AB_pXimSum.SetChannelWeight(2, 1. / 8.);
  AB_pXimSum.SetChannelWeight(3, 3. / 8.);

  AB_pXimSum.SetQ1Q2(-1);
  AB_pXimSum.SetPdgId(2212, 3122);
  AB_pXimSum.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
  AB_pXimSum.SetMaxRad(64);
  AB_pXimSum.SetMaxRho(32);

  CATSparameters cPotParsI0S0(CATSparameters::tPotential, 9, true);
  cPotParsI0S0.SetParameters(pXimPotParsI0S0);

  CATSparameters cPotParsI0S1(CATSparameters::tPotential, 9, true);
  cPotParsI0S1.SetParameters(pXimPotParsI0S1);

  CATSparameters cPotParsI1S0(CATSparameters::tPotential, 9, true);
  cPotParsI1S0.SetParameters(pXimPotParsI1S0);

  CATSparameters cPotParsI1S1(CATSparameters::tPotential, 9, true);
  cPotParsI1S1.SetParameters(pXimPotParsI1S1);

  AB_pXimSum.SetShortRangePotential(0, 0, fDlmPot, cPotParsI0S0);
  AB_pXimSum.SetShortRangePotential(1, 0, fDlmPot, cPotParsI0S1);
  AB_pXimSum.SetShortRangePotential(2, 0, fDlmPot, cPotParsI1S0);
  AB_pXimSum.SetShortRangePotential(3, 0, fDlmPot, cPotParsI1S1);

  AB_pXimSum.KillTheCat();

  double lampXi = 0.513;
  for (int ikS = 0; ikS < fMomBins; ++ikS) {
    fCFHistDefaultData->SetBinContent(
        ikS + 1, lampXi * (AB_pXimSum.GetCorrFun(ikS) - 1) + 1);
    fCFHistDefaultData->SetBinError(ikS + 1, fCFHistData->GetBinError(ikS + 1));
    fCFHistDefaultFit->SetBinContent(
        ikS + 1, lampXi * (AB_pXimSum.GetCorrFun(ikS) - 1) + 1);
    fCFHistDefaultFit->SetBinError(ikS + 1,
                                   fCFHistErrorFit->GetBinError(ikS + 1));
  }

  fCFHistDefaultData->GetYaxis()->SetRangeUser(
      0.8, fCFHistDefaultData->GetMaximum() * 1.5);
  fCFHistDefaultData->GetXaxis()->SetRangeUser(fkMin, fkMax);
  fCFHistDefaultData->SetStats(0);
  fCFHistDefaultData->SetLineWidth(1);
  fCFHistDefaultData->SetLineColor(1);
  fCFHistDefaultData->SetMarkerStyle(24);
  fCFHistDefaultData->SetMarkerColor(1);
  fCFHistDefaultData->SetMarkerSize(0.04);

  fCFHistDefaultFit->GetYaxis()->SetRangeUser(
      0.8, fCFHistDefaultFit->GetMaximum() * 1.5);
  fCFHistDefaultFit->GetXaxis()->SetRangeUser(fkMin, fkMax);
  fCFHistDefaultFit->SetStats(0);
  fCFHistDefaultFit->SetLineWidth(1);
  fCFHistDefaultFit->SetLineColor(1);
  fCFHistDefaultFit->SetMarkerStyle(24);
  fCFHistDefaultFit->SetMarkerColor(1);
  fCFHistDefaultFit->SetMarkerSize(0.04);

  fOutFile->cd();
  fCFHistDefaultData->Write();
  fCFHistDefaultFit->Write();
}

void PlayWithCats::GenerateCoulombOnly() {
  CATS AB_pXimCoulomb;
  double Pars_pXiSum[4] = { 0, 0, 0, fRadius };
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, fRadius);
  AB_pXimCoulomb.SetAnaSource(GaussSource, cPars);
  AB_pXimCoulomb.SetUseAnalyticSource(true);
  AB_pXimCoulomb.SetThetaDependentSource(false);
  AB_pXimCoulomb.SetExcludeFailedBins(false);
  AB_pXimCoulomb.SetMomBins(fMomBins, fkMin, fkMax);

  AB_pXimCoulomb.SetNumChannels(1);
  AB_pXimCoulomb.SetNumPW(0, 1);
  AB_pXimCoulomb.SetSpin(0, 0);    //I=0; S=0
  AB_pXimCoulomb.SetChannelWeight(0, 1);

  AB_pXimCoulomb.SetQ1Q2(-1);
  AB_pXimCoulomb.SetPdgId(2212, 3122);
  AB_pXimCoulomb.SetRedMass((fMass_p * fMass_Xim) / (fMass_p + fMass_Xim));
  AB_pXimCoulomb.SetMaxRad(64);
  AB_pXimCoulomb.SetMaxRho(32);

  AB_pXimCoulomb.KillTheCat();

  TString SumHistName = Form("CoulombOnly");
  TString SumHistTitle = Form("CoulombOnly");
  fCFCoulombOnly = new TH1F(SumHistName.Data(), SumHistTitle.Data(), fMomBins,
                            fkMin, fkMax);
  double lampXi = 0.513;
  for (int ikS = 0; ikS < fMomBins; ++ikS) {
    fCFCoulombOnly->SetBinContent(
        ikS + 1, lampXi * (AB_pXimCoulomb.GetCorrFun(ikS) - 1.) + 1.);
  }

  fCFCoulombOnly->GetXaxis()->SetRangeUser(0., 150.);
  fCFCoulombOnly->SetStats(0);
  fCFCoulombOnly->SetLineWidth(1);
  fCFCoulombOnly->SetLineColor(kRed);
  fCFCoulombOnly->SetMarkerColor(kRed);
  fCFCoulombOnly->SetMarkerSize(0.04);
  fOutFile->cd();
  fCFCoulombOnly->Write();
  return;
}
