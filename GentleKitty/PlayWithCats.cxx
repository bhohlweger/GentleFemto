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
#include "DreamPlot.h"

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

void PlayWithCats::GenerateSourceDistpxi(TFile* out) {
  std::cout << "Coulomb\n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  tidy->SetTau(1.65, 4.69);
  tidy->SetMass(1361.52, 1462.93);

  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 0.93;
  ppRadii[2] = 0.979453;

  CATS CoulombUp;

  tidy->GetCatsProtonXiMinus(&CoulombUp, nBins, xmin, xmax, TheSource,
                             TidyCats::pCoulomb, 0);

  CoulombUp.SetAnaSource(0, ppRadii[1]);
  CoulombUp.KillTheCat();

  auto grSource = new TGraph();
  DreamPlot::SetStyleGraph(grSource, 20, kBlue + 3);
  grSource->SetLineWidth(2);
  grSource->SetTitle(";#it{r} (fm); S(#it{r}) (fm^{-1})");
  grSource->SetName("pXiSource"); 
  for (double i = 0; i < 150; ++i) {
    grSource->SetPoint(i, i * 0.1, CoulombUp.EvaluateTheSource(0, i * 0.1, 0));
  }

  auto gaussFit =
    new TF1(
	    "gauspXi",
	    [&](double *x, double *p) {return
				       4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
				       std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
	    },
	    0, 10, 1);
  gaussFit->SetParameter(0, 1.3);
  gaussFit->SetNpx(1000);
  gaussFit->SetLineColor(kBlue + 2);
  gaussFit->SetLineWidth(2);
  gaussFit->SetLineStyle(2);

  grSource->Fit(gaussFit, "", "RQ", 0, 10);
  out->cd();
  grSource->Write("pXiSource");
  gaussFit->Write("gauspXi"); 
  return;
}

void PlayWithCats::GenerateSourceDistpp(TFile* out, double ks) {
  std::cout << "Coulomb\n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  tidy->SetTau(1.65, 4.69);
  tidy->SetMass(1361.52, 1462.93);
  tidy->SetkStarCutOff(ks); 
  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 1.2;
  ppRadii[2] = 0.979453;

  CATS CoulombUp;
  tidy->GetCatsProtonProton(&CoulombUp, nBins, xmin, xmax, TheSource); 

  CoulombUp.SetAnaSource(0, ppRadii[1]);
  CoulombUp.KillTheCat();

  auto grSource = new TGraph();
  DreamPlot::SetStyleGraph(grSource, 20, kBlue + 3);
  grSource->SetLineWidth(2);
  grSource->SetTitle(";#it{r} (fm); S(#it{r}) (fm^{-1})");
  grSource->SetName("gausPP"); 
  for (double i = 0; i < 150; ++i) {
    grSource->SetPoint(i, i * 0.1, CoulombUp.EvaluateTheSource(0, i * 0.1, 0));
  }

  auto gaussFit =
       new TF1(
           "gauspp",
           [&](double *x, double *p) {return
             4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
             std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
           },
           0, 10, 1);
   gaussFit->SetParameter(0, 1.3);
   gaussFit->SetNpx(1000);
   gaussFit->SetLineColor(kBlue + 2);
   gaussFit->SetLineWidth(2);
   gaussFit->SetLineStyle(2);
   
  grSource->Fit(gaussFit, "", "RQ", 0, 10);
  out->cd();
  grSource->Write("ppSource");
  gaussFit->Write("gauspp"); 
  return;
}

void PlayWithCats::GenerateSourceDistpL(TFile* out) {
  std::cout << "Coulomb\n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  tidy->SetTau(1.65, 4.69);
  tidy->SetMass(1361.52, 1462.93);

  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 1.2;
  ppRadii[2] = 0.979453;

  CATS CoulombUp;

  tidy->GetCatsProtonLambda(&CoulombUp, nBins, xmin, xmax, TheSource, TidyCats::pNLOWF); 

  CoulombUp.SetAnaSource(0, ppRadii[1]);
  CoulombUp.KillTheCat();

  auto grSource = new TGraph();
  DreamPlot::SetStyleGraph(grSource, 20, kBlue + 3);
  grSource->SetLineWidth(2);
  grSource->SetTitle(";#it{r} (fm); S(#it{r}) (fm^{-1})");
  grSource->SetName("pLSource"); 
  for (double i = 0; i < 150; ++i) {
    grSource->SetPoint(i, i * 0.1, CoulombUp.EvaluateTheSource(0, i * 0.1, 0));
  }

  auto gaussFit =
       new TF1(
           "gauspL",
           [&](double *x, double *p) {return
             4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
             std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
           },
           0, 10, 1);
   gaussFit->SetParameter(0, 1.3);
   gaussFit->SetNpx(1000);
   gaussFit->SetLineColor(kBlue + 2);
   gaussFit->SetLineWidth(2);
   gaussFit->SetLineStyle(2);
  
  grSource->Fit(gaussFit, "", "RQ", 0, 10);
  out->cd();
  grSource->Write("pLSource");
  gaussFit->Write("gauspL"); 
  return;
}



void PlayWithCats::GenerateYukiCurves(TFile* out) {
  out->cd(); 
  std::vector<std::vector<float>> Yuki_0_97 =
    {
      {1.00E+00,4.55E+01,2.16E+00,3.88E+01,5.27E+00,3.95E-01},
      {2.00E+00,2.26E+01,1.08E+00,1.87E+01,2.63E+00,1.97E-01},
      {3.00E+00,1.50E+01,7.14E-01,1.21E+01,1.75E+00,1.31E-01},
      {4.00E+00,1.12E+01,5.33E-01,8.83E+00,1.31E+00,9.82E-02},
      {5.00E+00,9.00E+00,4.24E-01,6.88E+00,1.05E+00,7.84E-02},
      {6.00E+00,7.53E+00,3.52E-01,5.59E+00,8.72E-01,6.54E-02},
      {7.00E+00,6.50E+00,3.02E-01,4.69E+00,7.49E-01,5.62E-02},
      {8.00E+00,5.74E+00,2.64E-01,4.03E+00,6.58E-01,4.93E-02},
      {9.00E+00,5.16E+00,2.36E-01,3.53E+00,5.88E-01,4.41E-02},
      {1.00E+01,4.71E+00,2.13E-01,3.13E+00,5.32E-01,3.99E-02},
      {1.10E+01,4.34E+00,1.95E-01,2.82E+00,4.87E-01,3.66E-02},
      {1.20E+01,4.04E+00,1.80E-01,2.56E+00,4.50E-01,3.38E-02},
      {1.30E+01,3.79E+00,1.67E-01,2.34E+00,4.18E-01,3.14E-02},
      {1.40E+01,3.58E+00,1.57E-01,2.16E+00,3.91E-01,2.94E-02},
      {1.50E+01,3.40E+00,1.47E-01,2.00E+00,3.68E-01,2.77E-02},
      {1.60E+01,3.24E+00,1.39E-01,1.87E+00,3.47E-01,2.62E-02},
      {1.70E+01,3.10E+00,1.32E-01,1.75E+00,3.29E-01,2.48E-02},
      {1.80E+01,2.98E+00,1.26E-01,1.64E+00,3.13E-01,2.36E-02},
      {1.90E+01,2.87E+00,1.20E-01,1.55E+00,2.98E-01,2.25E-02},
      {2.00E+01,2.77E+00,1.15E-01,1.46E+00,2.85E-01,2.15E-02},
      {2.10E+01,2.68E+00,1.10E-01,1.39E+00,2.73E-01,2.06E-02},
      {2.20E+01,2.60E+00,1.06E-01,1.32E+00,2.62E-01,1.98E-02},
      {2.30E+01,2.53E+00,1.02E-01,1.25E+00,2.51E-01,1.91E-02},
      {2.40E+01,2.46E+00,9.88E-02,1.20E+00,2.42E-01,1.84E-02},
      {2.50E+01,2.39E+00,9.55E-02,1.14E+00,2.33E-01,1.77E-02},
      {2.60E+01,2.34E+00,9.24E-02,1.09E+00,2.25E-01,1.71E-02},
      {2.70E+01,2.28E+00,8.95E-02,1.05E+00,2.17E-01,1.65E-02},
      {2.80E+01,2.23E+00,8.68E-02,1.00E+00,2.10E-01,1.60E-02},
      {2.90E+01,2.18E+00,8.43E-02,9.65E-01,2.03E-01,1.55E-02},
      {3.00E+01,2.14E+00,8.19E-02,9.27E-01,1.97E-01,1.50E-02},
      {3.10E+01,2.10E+00,7.96E-02,8.92E-01,1.91E-01,1.46E-02},
      {3.20E+01,2.06E+00,7.74E-02,8.59E-01,1.85E-01,1.42E-02},
      {3.30E+01,2.02E+00,7.54E-02,8.27E-01,1.79E-01,1.38E-02},
      {3.40E+01,1.99E+00,7.34E-02,7.97E-01,1.74E-01,1.34E-02},
      {3.50E+01,1.95E+00,7.16E-02,7.69E-01,1.69E-01,1.30E-02},
      {3.60E+01,1.92E+00,6.98E-02,7.43E-01,1.64E-01,1.27E-02},
      {3.70E+01,1.89E+00,6.81E-02,7.17E-01,1.60E-01,1.23E-02},
      {3.80E+01,1.86E+00,6.64E-02,6.93E-01,1.55E-01,1.20E-02},
      {3.90E+01,1.83E+00,6.48E-02,6.70E-01,1.51E-01,1.17E-02},
      {4.00E+01,1.81E+00,6.33E-02,6.48E-01,1.47E-01,1.14E-02},
      {4.10E+01,1.78E+00,6.18E-02,6.26E-01,1.43E-01,1.11E-02},
      {4.20E+01,1.76E+00,6.04E-02,6.06E-01,1.40E-01,1.09E-02},
      {4.30E+01,1.73E+00,5.90E-02,5.87E-01,1.36E-01,1.06E-02},
      {4.40E+01,1.71E+00,5.76E-02,5.68E-01,1.33E-01,1.04E-02},
      {4.50E+01,1.69E+00,5.63E-02,5.50E-01,1.29E-01,1.01E-02},
      {4.60E+01,1.67E+00,5.50E-02,5.33E-01,1.26E-01,9.91E-03},
      {4.70E+01,1.65E+00,5.38E-02,5.16E-01,1.23E-01,9.69E-03},
      {4.80E+01,1.63E+00,5.26E-02,5.00E-01,1.20E-01,9.48E-03},
      {4.90E+01,1.61E+00,5.14E-02,4.84E-01,1.17E-01,9.27E-03},
      {5.00E+01,1.59E+00,5.02E-02,4.69E-01,1.15E-01,9.07E-03},
      {5.10E+01,1.58E+00,4.91E-02,4.55E-01,1.12E-01,8.88E-03},
      {5.20E+01,1.56E+00,4.80E-02,4.41E-01,1.09E-01,8.69E-03},
      {5.30E+01,1.54E+00,4.69E-02,4.28E-01,1.07E-01,8.52E-03},
      {5.40E+01,1.53E+00,4.58E-02,4.15E-01,1.04E-01,8.34E-03},
      {5.50E+01,1.51E+00,4.48E-02,4.02E-01,1.02E-01,8.17E-03},
      {5.60E+01,1.50E+00,4.38E-02,3.90E-01,9.98E-02,8.01E-03},
      {5.70E+01,1.48E+00,4.28E-02,3.78E-01,9.77E-02,7.85E-03},
      {5.80E+01,1.47E+00,4.18E-02,3.66E-01,9.55E-02,7.70E-03},
      {5.90E+01,1.46E+00,4.09E-02,3.55E-01,9.35E-02,7.55E-03},
      {6.00E+01,1.44E+00,3.99E-02,3.44E-01,9.15E-02,7.41E-03},
      {6.10E+01,1.43E+00,3.90E-02,3.34E-01,8.95E-02,7.27E-03},
      {6.20E+01,1.42E+00,3.81E-02,3.24E-01,8.76E-02,7.14E-03},
      {6.30E+01,1.41E+00,3.72E-02,3.14E-01,8.58E-02,7.01E-03},
      {6.40E+01,1.40E+00,3.64E-02,3.04E-01,8.40E-02,6.88E-03},
      {6.50E+01,1.38E+00,3.55E-02,2.95E-01,8.23E-02,6.75E-03},
      {6.60E+01,1.37E+00,3.47E-02,2.86E-01,8.06E-02,6.63E-03},
      {6.70E+01,1.36E+00,3.39E-02,2.77E-01,7.90E-02,6.52E-03},
      {6.80E+01,1.35E+00,3.31E-02,2.68E-01,7.74E-02,6.40E-03},
      {6.90E+01,1.34E+00,3.23E-02,2.60E-01,7.59E-02,6.29E-03},
      {7.00E+01,1.33E+00,3.15E-02,2.52E-01,7.44E-02,6.19E-03},
      {7.10E+01,1.32E+00,3.08E-02,2.44E-01,7.29E-02,6.08E-03},
      {7.20E+01,1.31E+00,3.00E-02,2.36E-01,7.15E-02,5.98E-03},
      {7.30E+01,1.30E+00,2.93E-02,2.28E-01,7.01E-02,5.88E-03},
      {7.40E+01,1.30E+00,2.86E-02,2.21E-01,6.88E-02,5.79E-03},
      {7.50E+01,1.29E+00,2.79E-02,2.14E-01,6.75E-02,5.69E-03},
      {7.60E+01,1.28E+00,2.72E-02,2.07E-01,6.62E-02,5.60E-03},
      {7.70E+01,1.27E+00,2.65E-02,2.00E-01,6.50E-02,5.51E-03},
      {7.80E+01,1.26E+00,2.59E-02,1.94E-01,6.38E-02,5.42E-03},
      {7.90E+01,1.26E+00,2.52E-02,1.87E-01,6.26E-02,5.34E-03},
      {8.00E+01,1.25E+00,2.46E-02,1.81E-01,6.14E-02,5.26E-03},
      {8.10E+01,1.24E+00,2.40E-02,1.75E-01,6.03E-02,5.18E-03},
      {8.20E+01,1.23E+00,2.34E-02,1.69E-01,5.92E-02,5.10E-03},
      {8.30E+01,1.23E+00,2.28E-02,1.63E-01,5.82E-02,5.02E-03},
      {8.40E+01,1.22E+00,2.22E-02,1.57E-01,5.71E-02,4.95E-03},
      {8.50E+01,1.21E+00,2.16E-02,1.52E-01,5.61E-02,4.88E-03},
      {8.60E+01,1.21E+00,2.11E-02,1.47E-01,5.51E-02,4.80E-03},
      {8.70E+01,1.20E+00,2.05E-02,1.41E-01,5.42E-02,4.74E-03},
      {8.80E+01,1.19E+00,2.00E-02,1.36E-01,5.32E-02,4.67E-03},
      {8.90E+01,1.19E+00,1.95E-02,1.31E-01,5.23E-02,4.60E-03},
      {9.00E+01,1.18E+00,1.90E-02,1.27E-01,5.14E-02,4.54E-03},
      {9.10E+01,1.18E+00,1.85E-02,1.22E-01,5.05E-02,4.47E-03},
      {9.20E+01,1.17E+00,1.80E-02,1.17E-01,4.97E-02,4.41E-03},
      {9.30E+01,1.17E+00,1.75E-02,1.13E-01,4.88E-02,4.35E-03},
      {9.40E+01,1.16E+00,1.70E-02,1.08E-01,4.80E-02,4.29E-03},
      {9.50E+01,1.16E+00,1.66E-02,1.04E-01,4.72E-02,4.24E-03},
      {9.60E+01,1.15E+00,1.61E-02,1.00E-01,4.65E-02,4.18E-03},
      {9.70E+01,1.15E+00,1.57E-02,9.60E-02,4.57E-02,4.13E-03},
      {9.80E+01,1.14E+00,1.53E-02,9.21E-02,4.50E-02,4.07E-03},
      {9.90E+01,1.14E+00,1.48E-02,8.83E-02,4.42E-02,4.02E-03},
      {1.00E+02,1.13E+00,1.44E-02,8.46E-02,4.35E-02,3.97E-03},
      {1.01E+02,1.13E+00,1.40E-02,8.10E-02,4.28E-02,3.92E-03},
      {1.02E+02,1.12E+00,1.36E-02,7.75E-02,4.22E-02,3.87E-03},
      {1.03E+02,1.12E+00,1.33E-02,7.40E-02,4.15E-02,3.82E-03},
      {1.04E+02,1.12E+00,1.29E-02,7.07E-02,4.08E-02,3.77E-03},
      {1.05E+02,1.11E+00,1.25E-02,6.75E-02,4.02E-02,3.73E-03},
      {1.06E+02,1.11E+00,1.22E-02,6.43E-02,3.96E-02,3.68E-03},
      {1.07E+02,1.10E+00,1.18E-02,6.13E-02,3.90E-02,3.64E-03},
      {1.08E+02,1.10E+00,1.15E-02,5.83E-02,3.84E-02,3.59E-03},
      {1.09E+02,1.10E+00,1.12E-02,5.54E-02,3.78E-02,3.55E-03},
      {1.10E+02,1.09E+00,1.09E-02,5.25E-02,3.72E-02,3.51E-03},
      {1.11E+02,1.09E+00,1.06E-02,4.98E-02,3.67E-02,3.47E-03},
      {1.12E+02,1.09E+00,1.03E-02,4.71E-02,3.61E-02,3.43E-03},
      {1.13E+02,1.08E+00,9.96E-03,4.45E-02,3.56E-02,3.39E-03},
      {1.14E+02,1.08E+00,9.68E-03,4.20E-02,3.51E-02,3.35E-03},
      {1.15E+02,1.08E+00,9.40E-03,3.95E-02,3.46E-02,3.31E-03},
      {1.16E+02,1.07E+00,9.13E-03,3.72E-02,3.41E-02,3.28E-03},
      {1.17E+02,1.07E+00,8.86E-03,3.49E-02,3.36E-02,3.24E-03},
      {1.18E+02,1.07E+00,8.60E-03,3.26E-02,3.31E-02,3.21E-03},
      {1.19E+02,1.07E+00,8.35E-03,3.04E-02,3.26E-02,3.17E-03},
      {1.20E+02,1.06E+00,8.11E-03,2.83E-02,3.22E-02,3.14E-03},
      {1.21E+02,1.06E+00,7.87E-03,2.63E-02,3.17E-02,3.10E-03},
      {1.22E+02,1.06E+00,7.64E-03,2.43E-02,3.13E-02,3.07E-03},
      {1.23E+02,1.06E+00,7.42E-03,2.23E-02,3.08E-02,3.04E-03},
      {1.24E+02,1.05E+00,7.20E-03,2.04E-02,3.04E-02,3.01E-03},
      {1.25E+02,1.05E+00,6.99E-03,1.86E-02,3.00E-02,2.98E-03},
      {1.26E+02,1.05E+00,6.79E-03,1.69E-02,2.96E-02,2.95E-03},
      {1.27E+02,1.05E+00,6.59E-03,1.51E-02,2.92E-02,2.92E-03},
      {1.28E+02,1.05E+00,6.39E-03,1.35E-02,2.88E-02,2.89E-03},
      {1.29E+02,1.04E+00,6.20E-03,1.19E-02,2.84E-02,2.86E-03},
      {1.30E+02,1.04E+00,6.02E-03,1.03E-02,2.80E-02,2.83E-03},
      {1.31E+02,1.04E+00,5.85E-03,8.82E-03,2.76E-02,2.80E-03},
      {1.32E+02,1.04E+00,5.67E-03,7.36E-03,2.73E-02,2.77E-03},
      {1.33E+02,1.04E+00,5.51E-03,5.95E-03,2.69E-02,2.75E-03},
      {1.34E+02,1.03E+00,5.35E-03,4.58E-03,2.66E-02,2.72E-03},
      {1.35E+02,1.03E+00,5.19E-03,3.26E-03,2.62E-02,2.69E-03},
      {1.36E+02,1.03E+00,5.04E-03,1.98E-03,2.59E-02,2.67E-03},
      {1.37E+02,1.03E+00,4.89E-03,7.47E-04,2.56E-02,2.64E-03},
      {1.38E+02,1.03E+00,4.75E-03,-4.47E-04,2.52E-02,2.62E-03},
      {1.39E+02,1.03E+00,4.61E-03,-1.60E-03,2.49E-02,2.60E-03},
      {1.40E+02,1.02E+00,4.48E-03,-2.71E-03,2.46E-02,2.57E-03},
      {1.41E+02,1.02E+00,4.35E-03,-3.79E-03,2.43E-02,2.55E-03},
      {1.42E+02,1.02E+00,4.23E-03,-4.83E-03,2.40E-02,2.52E-03},
      {1.43E+02,1.02E+00,4.11E-03,-5.83E-03,2.37E-02,2.50E-03},
      {1.44E+02,1.02E+00,3.99E-03,-6.79E-03,2.34E-02,2.48E-03},
      {1.45E+02,1.02E+00,3.88E-03,-7.72E-03,2.31E-02,2.46E-03},
      {1.46E+02,1.02E+00,3.77E-03,-8.62E-03,2.28E-02,2.44E-03},
      {1.47E+02,1.02E+00,3.67E-03,-9.48E-03,2.25E-02,2.42E-03},
      {1.48E+02,1.01E+00,3.57E-03,-1.03E-02,2.23E-02,2.39E-03},
      {1.49E+02,1.01E+00,3.47E-03,-1.11E-02,2.20E-02,2.37E-03},
      {1.50E+02,1.01E+00,3.38E-03,-1.19E-02,2.17E-02,2.35E-03},
      {1.51E+02,1.01E+00,3.29E-03,-1.26E-02,2.15E-02,2.33E-03},
      {1.52E+02,1.01E+00,3.21E-03,-1.33E-02,2.12E-02,2.31E-03},
      {1.53E+02,1.01E+00,3.13E-03,-1.40E-02,2.10E-02,2.30E-03},
      {1.54E+02,1.01E+00,3.05E-03,-1.47E-02,2.07E-02,2.28E-03},
      {1.55E+02,1.01E+00,2.97E-03,-1.53E-02,2.05E-02,2.26E-03},
      {1.56E+02,1.01E+00,2.90E-03,-1.59E-02,2.03E-02,2.24E-03},
      {1.57E+02,1.01E+00,2.83E-03,-1.65E-02,2.00E-02,2.22E-03},
      {1.58E+02,1.01E+00,2.77E-03,-1.70E-02,1.98E-02,2.20E-03},
      {1.59E+02,1.00E+00,2.70E-03,-1.76E-02,1.96E-02,2.19E-03},
      {1.60E+02,1.00E+00,2.64E-03,-1.81E-02,1.94E-02,2.17E-03},
      {1.61E+02,1.00E+00,2.59E-03,-1.86E-02,1.91E-02,2.15E-03},
      {1.62E+02,1.00E+00,2.53E-03,-1.90E-02,1.89E-02,2.13E-03},
      {1.63E+02,1.00E+00,2.48E-03,-1.95E-02,1.87E-02,2.12E-03},
      {1.64E+02,1.00E+00,2.43E-03,-1.99E-02,1.85E-02,2.10E-03},
      {1.65E+02,1.00E+00,2.39E-03,-2.03E-02,1.83E-02,2.09E-03},
      {1.66E+02,9.99E-01,2.34E-03,-2.07E-02,1.81E-02,2.07E-03},
      {1.67E+02,9.99E-01,2.30E-03,-2.11E-02,1.79E-02,2.05E-03},
      {1.68E+02,9.98E-01,2.26E-03,-2.14E-02,1.77E-02,2.04E-03},
      {1.69E+02,9.98E-01,2.22E-03,-2.17E-02,1.75E-02,2.02E-03},
      {1.70E+02,9.97E-01,2.19E-03,-2.20E-02,1.73E-02,2.01E-03},
      {1.71E+02,9.97E-01,2.15E-03,-2.23E-02,1.72E-02,1.99E-03},
      {1.72E+02,9.96E-01,2.12E-03,-2.26E-02,1.70E-02,1.98E-03},
      {1.73E+02,9.96E-01,2.09E-03,-2.29E-02,1.68E-02,1.97E-03},
      {1.74E+02,9.95E-01,2.06E-03,-2.31E-02,1.66E-02,1.95E-03},
      {1.75E+02,9.95E-01,2.03E-03,-2.34E-02,1.64E-02,1.94E-03},
      {1.76E+02,9.95E-01,2.00E-03,-2.36E-02,1.63E-02,1.92E-03},
      {1.77E+02,9.94E-01,1.98E-03,-2.38E-02,1.61E-02,1.91E-03},
      {1.78E+02,9.94E-01,1.95E-03,-2.40E-02,1.59E-02,1.90E-03},
      {1.79E+02,9.93E-01,1.93E-03,-2.42E-02,1.58E-02,1.88E-03},
      {1.80E+02,9.93E-01,1.91E-03,-2.43E-02,1.56E-02,1.87E-03},
      {1.81E+02,9.93E-01,1.89E-03,-2.45E-02,1.55E-02,1.86E-03},
      {1.82E+02,9.93E-01,1.86E-03,-2.46E-02,1.53E-02,1.85E-03},
      {1.83E+02,9.92E-01,1.84E-03,-2.48E-02,1.52E-02,1.83E-03},
      {1.84E+02,9.92E-01,1.82E-03,-2.49E-02,1.50E-02,1.82E-03},
      {1.85E+02,9.92E-01,1.81E-03,-2.50E-02,1.49E-02,1.81E-03},
      {1.86E+02,9.91E-01,1.79E-03,-2.51E-02,1.47E-02,1.80E-03},
      {1.87E+02,9.91E-01,1.77E-03,-2.52E-02,1.46E-02,1.79E-03},
      {1.88E+02,9.91E-01,1.75E-03,-2.53E-02,1.44E-02,1.77E-03},
      {1.89E+02,9.91E-01,1.74E-03,-2.54E-02,1.43E-02,1.76E-03},
      {1.90E+02,9.90E-01,1.72E-03,-2.54E-02,1.41E-02,1.75E-03},
      {1.91E+02,9.90E-01,1.70E-03,-2.55E-02,1.40E-02,1.74E-03},
      {1.92E+02,9.90E-01,1.69E-03,-2.55E-02,1.39E-02,1.73E-03},
      {1.93E+02,9.90E-01,1.67E-03,-2.56E-02,1.37E-02,1.72E-03},
      {1.94E+02,9.90E-01,1.66E-03,-2.56E-02,1.36E-02,1.71E-03},
      {1.95E+02,9.90E-01,1.64E-03,-2.56E-02,1.35E-02,1.70E-03},
      {1.96E+02,9.89E-01,1.63E-03,-2.57E-02,1.34E-02,1.69E-03},
      {1.97E+02,9.89E-01,1.61E-03,-2.57E-02,1.32E-02,1.68E-03},
      {1.98E+02,9.89E-01,1.60E-03,-2.57E-02,1.31E-02,1.67E-03},
      {1.99E+02,9.89E-01,1.58E-03,-2.57E-02,1.30E-02,1.66E-03},
      {2.00E+02,9.89E-01,1.57E-03,-2.57E-02,1.29E-02,1.65E-03},
      {2.01E+02,9.89E-01,1.55E-03,-2.57E-02,1.28E-02,1.64E-03},
      {2.02E+02,9.89E-01,1.54E-03,-2.57E-02,1.26E-02,1.63E-03},
      {2.03E+02,9.89E-01,1.53E-03,-2.56E-02,1.25E-02,1.62E-03},
      {2.04E+02,9.88E-01,1.51E-03,-2.56E-02,1.24E-02,1.61E-03},
      {2.05E+02,9.88E-01,1.50E-03,-2.56E-02,1.23E-02,1.60E-03},
      {2.06E+02,9.88E-01,1.49E-03,-2.55E-02,1.22E-02,1.59E-03},
      {2.07E+02,9.88E-01,1.47E-03,-2.55E-02,1.21E-02,1.58E-03},
      {2.08E+02,9.88E-01,1.46E-03,-2.55E-02,1.20E-02,1.57E-03},
      {2.09E+02,9.88E-01,1.45E-03,-2.54E-02,1.19E-02,1.56E-03},
      {2.10E+02,9.88E-01,1.44E-03,-2.54E-02,1.18E-02,1.56E-03},
      {2.11E+02,9.88E-01,1.42E-03,-2.53E-02,1.17E-02,1.55E-03},
      {2.12E+02,9.88E-01,1.41E-03,-2.52E-02,1.16E-02,1.54E-03},
      {2.13E+02,9.88E-01,1.40E-03,-2.52E-02,1.15E-02,1.53E-03},
      {2.14E+02,9.88E-01,1.38E-03,-2.51E-02,1.14E-02,1.52E-03},
      {2.15E+02,9.88E-01,1.37E-03,-2.50E-02,1.13E-02,1.51E-03},
      {2.16E+02,9.88E-01,1.36E-03,-2.50E-02,1.12E-02,1.51E-03},
      {2.17E+02,9.88E-01,1.35E-03,-2.49E-02,1.11E-02,1.50E-03},
      {2.18E+02,9.88E-01,1.34E-03,-2.48E-02,1.10E-02,1.49E-03},
      {2.19E+02,9.88E-01,1.32E-03,-2.47E-02,1.09E-02,1.48E-03},
      {2.20E+02,9.88E-01,1.31E-03,-2.46E-02,1.08E-02,1.47E-03},
      {2.21E+02,9.88E-01,1.30E-03,-2.46E-02,1.07E-02,1.47E-03},
      {2.22E+02,9.88E-01,1.29E-03,-2.45E-02,1.06E-02,1.46E-03},
      {2.23E+02,9.88E-01,1.27E-03,-2.44E-02,1.05E-02,1.45E-03},
      {2.24E+02,9.88E-01,1.26E-03,-2.43E-02,1.05E-02,1.44E-03},
      {2.25E+02,9.88E-01,1.25E-03,-2.42E-02,1.04E-02,1.44E-03},
      {2.26E+02,9.88E-01,1.24E-03,-2.41E-02,1.03E-02,1.43E-03},
      {2.27E+02,9.88E-01,1.23E-03,-2.40E-02,1.02E-02,1.42E-03},
      {2.28E+02,9.88E-01,1.22E-03,-2.39E-02,1.01E-02,1.42E-03},
      {2.29E+02,9.88E-01,1.21E-03,-2.38E-02,1.01E-02,1.41E-03},
      {2.30E+02,9.88E-01,1.19E-03,-2.37E-02,9.97E-03,1.40E-03},
      {2.31E+02,9.88E-01,1.18E-03,-2.36E-02,9.89E-03,1.40E-03},
      {2.32E+02,9.88E-01,1.17E-03,-2.35E-02,9.82E-03,1.39E-03},
      {2.33E+02,9.88E-01,1.16E-03,-2.34E-02,9.74E-03,1.38E-03},
      {2.34E+02,9.88E-01,1.15E-03,-2.33E-02,9.67E-03,1.38E-03},
      {2.35E+02,9.88E-01,1.14E-03,-2.32E-02,9.59E-03,1.37E-03},
      {2.36E+02,9.88E-01,1.13E-03,-2.31E-02,9.52E-03,1.36E-03},
      {2.37E+02,9.88E-01,1.12E-03,-2.30E-02,9.45E-03,1.36E-03},
      {2.38E+02,9.88E-01,1.11E-03,-2.28E-02,9.37E-03,1.35E-03},
      {2.39E+02,9.88E-01,1.09E-03,-2.27E-02,9.30E-03,1.34E-03},
      {2.40E+02,9.88E-01,1.08E-03,-2.26E-02,9.23E-03,1.34E-03},
      {2.41E+02,9.88E-01,1.07E-03,-2.25E-02,9.16E-03,1.33E-03},
      {2.42E+02,9.88E-01,1.06E-03,-2.24E-02,9.10E-03,1.33E-03},
      {2.43E+02,9.88E-01,1.05E-03,-2.23E-02,9.03E-03,1.32E-03},
      {2.44E+02,9.88E-01,1.04E-03,-2.22E-02,8.96E-03,1.31E-03},
      {2.45E+02,9.88E-01,1.03E-03,-2.21E-02,8.90E-03,1.31E-03},
      {2.46E+02,9.88E-01,1.02E-03,-2.19E-02,8.83E-03,1.30E-03},
      {2.47E+02,9.88E-01,1.01E-03,-2.18E-02,8.77E-03,1.30E-03},
      {2.48E+02,9.88E-01,1.00E-03,-2.17E-02,8.70E-03,1.29E-03},
      {2.49E+02,9.88E-01,9.91E-04,-2.16E-02,8.64E-03,1.28E-03},
      {2.50E+02,9.88E-01,9.82E-04,-2.15E-02,8.58E-03,1.28E-03},
      {2.51E+02,9.88E-01,9.72E-04,-2.14E-02,8.52E-03,1.27E-03},
      {2.52E+02,9.88E-01,9.62E-04,-2.13E-02,8.46E-03,1.27E-03},
      {2.53E+02,9.89E-01,9.52E-04,-2.11E-02,8.40E-03,1.26E-03},
      {2.54E+02,9.89E-01,9.43E-04,-2.10E-02,8.34E-03,1.26E-03},
      {2.55E+02,9.89E-01,9.33E-04,-2.09E-02,8.28E-03,1.25E-03},
      {2.56E+02,9.89E-01,9.24E-04,-2.08E-02,8.22E-03,1.25E-03},
      {2.57E+02,9.89E-01,9.15E-04,-2.07E-02,8.16E-03,1.24E-03},
      {2.58E+02,9.89E-01,9.05E-04,-2.06E-02,8.11E-03,1.24E-03},
      {2.59E+02,9.89E-01,8.96E-04,-2.04E-02,8.05E-03,1.23E-03},
      {2.60E+02,9.89E-01,8.87E-04,-2.03E-02,8.00E-03,1.23E-03},
      {2.61E+02,9.89E-01,8.78E-04,-2.02E-02,7.94E-03,1.22E-03},
      {2.62E+02,9.89E-01,8.69E-04,-2.01E-02,7.89E-03,1.22E-03},
      {2.63E+02,9.89E-01,8.60E-04,-2.00E-02,7.83E-03,1.21E-03},
      {2.64E+02,9.89E-01,8.52E-04,-1.99E-02,7.78E-03,1.21E-03},
      {2.65E+02,9.89E-01,8.43E-04,-1.97E-02,7.73E-03,1.20E-03},
      {2.66E+02,9.89E-01,8.34E-04,-1.96E-02,7.68E-03,1.20E-03},
      {2.67E+02,9.89E-01,8.26E-04,-1.95E-02,7.63E-03,1.19E-03},
      {2.68E+02,9.89E-01,8.17E-04,-1.94E-02,7.58E-03,1.19E-03},
      {2.69E+02,9.89E-01,8.09E-04,-1.93E-02,7.53E-03,1.18E-03},
      {2.70E+02,9.89E-01,8.01E-04,-1.92E-02,7.48E-03,1.18E-03},
      {2.71E+02,9.90E-01,7.92E-04,-1.91E-02,7.43E-03,1.17E-03},
      {2.72E+02,9.90E-01,7.84E-04,-1.90E-02,7.38E-03,1.17E-03},
      {2.73E+02,9.90E-01,7.76E-04,-1.88E-02,7.33E-03,1.17E-03},
      {2.74E+02,9.90E-01,7.68E-04,-1.87E-02,7.28E-03,1.16E-03},
      {2.75E+02,9.90E-01,7.60E-04,-1.86E-02,7.24E-03,1.16E-03},
      {2.76E+02,9.90E-01,7.52E-04,-1.85E-02,7.19E-03,1.15E-03},
      {2.77E+02,9.90E-01,7.44E-04,-1.84E-02,7.14E-03,1.15E-03},
      {2.78E+02,9.90E-01,7.37E-04,-1.83E-02,7.10E-03,1.14E-03},
      {2.79E+02,9.90E-01,7.29E-04,-1.82E-02,7.05E-03,1.14E-03},
      {2.80E+02,9.90E-01,7.22E-04,-1.81E-02,7.01E-03,1.13E-03},
      {2.81E+02,9.90E-01,7.14E-04,-1.80E-02,6.97E-03,1.13E-03},
      {2.82E+02,9.90E-01,7.07E-04,-1.79E-02,6.92E-03,1.13E-03},
      {2.83E+02,9.90E-01,6.99E-04,-1.77E-02,6.88E-03,1.12E-03},
      {2.84E+02,9.90E-01,6.92E-04,-1.76E-02,6.84E-03,1.12E-03},
      {2.85E+02,9.90E-01,6.85E-04,-1.75E-02,6.79E-03,1.11E-03},
      {2.86E+02,9.90E-01,6.78E-04,-1.74E-02,6.75E-03,1.11E-03},
      {2.87E+02,9.91E-01,6.71E-04,-1.73E-02,6.71E-03,1.11E-03},
      {2.88E+02,9.91E-01,6.64E-04,-1.72E-02,6.67E-03,1.10E-03},
      {2.89E+02,9.91E-01,6.57E-04,-1.71E-02,6.63E-03,1.10E-03},
      {2.90E+02,9.91E-01,6.50E-04,-1.70E-02,6.59E-03,1.09E-03},
      {2.91E+02,9.91E-01,6.43E-04,-1.69E-02,6.55E-03,1.09E-03},
      {2.92E+02,9.91E-01,6.36E-04,-1.68E-02,6.51E-03,1.09E-03},
      {2.93E+02,9.91E-01,6.30E-04,-1.67E-02,6.47E-03,1.08E-03},
      {2.94E+02,9.91E-01,6.23E-04,-1.66E-02,6.43E-03,1.08E-03},
      {2.95E+02,9.91E-01,6.16E-04,-1.65E-02,6.39E-03,1.08E-03},
      {2.96E+02,9.91E-01,6.10E-04,-1.64E-02,6.36E-03,1.07E-03},
      {2.97E+02,9.91E-01,6.04E-04,-1.63E-02,6.32E-03,1.07E-03},
      {2.98E+02,9.91E-01,5.97E-04,-1.62E-02,6.28E-03,1.06E-03},
      {2.99E+02,9.91E-01,5.91E-04,-1.61E-02,6.25E-03,1.06E-03},
      {3.00E+02,9.91E-01,5.85E-04,-1.60E-02,6.21E-03,1.06E-03}
    };

  std::vector<std::vector<float>> Yuki_1_02 =
    {
      {1.00E+00,4.41E+01,2.07E+00,3.78E+01,4.90E+00,3.60E-01},
      {2.00E+00,2.19E+01,1.03E+00,1.83E+01,2.45E+00,1.80E-01},
      {3.00E+00,1.46E+01,6.84E-01,1.18E+01,1.63E+00,1.19E-01},
      {4.00E+00,1.09E+01,5.10E-01,8.61E+00,1.22E+00,8.94E-02},
      {5.00E+00,8.74E+00,4.06E-01,6.70E+00,9.72E-01,7.14E-02},
      {6.00E+00,7.31E+00,3.37E-01,5.44E+00,8.10E-01,5.95E-02},
      {7.00E+00,6.31E+00,2.88E-01,4.56E+00,6.96E-01,5.12E-02},
      {8.00E+00,5.57E+00,2.52E-01,3.92E+00,6.11E-01,4.49E-02},
      {9.00E+00,5.01E+00,2.25E-01,3.43E+00,5.46E-01,4.02E-02},
      {1.00E+01,4.57E+00,2.03E-01,3.04E+00,4.94E-01,3.64E-02},
      {1.10E+01,4.22E+00,1.86E-01,2.73E+00,4.52E-01,3.33E-02},
      {1.20E+01,3.93E+00,1.72E-01,2.48E+00,4.18E-01,3.08E-02},
      {1.30E+01,3.68E+00,1.60E-01,2.27E+00,3.88E-01,2.86E-02},
      {1.40E+01,3.48E+00,1.49E-01,2.09E+00,3.63E-01,2.68E-02},
      {1.50E+01,3.30E+00,1.41E-01,1.94E+00,3.42E-01,2.52E-02},
      {1.60E+01,3.15E+00,1.33E-01,1.80E+00,3.23E-01,2.38E-02},
      {1.70E+01,3.01E+00,1.26E-01,1.69E+00,3.06E-01,2.26E-02},
      {1.80E+01,2.90E+00,1.20E-01,1.58E+00,2.91E-01,2.15E-02},
      {1.90E+01,2.79E+00,1.15E-01,1.49E+00,2.77E-01,2.05E-02},
      {2.00E+01,2.69E+00,1.10E-01,1.41E+00,2.65E-01,1.96E-02},
      {2.10E+01,2.61E+00,1.05E-01,1.33E+00,2.53E-01,1.88E-02},
      {2.20E+01,2.53E+00,1.01E-01,1.27E+00,2.43E-01,1.81E-02},
      {2.30E+01,2.46E+00,9.78E-02,1.21E+00,2.34E-01,1.74E-02},
      {2.40E+01,2.39E+00,9.44E-02,1.15E+00,2.25E-01,1.67E-02},
      {2.50E+01,2.33E+00,9.12E-02,1.10E+00,2.17E-01,1.61E-02},
      {2.60E+01,2.27E+00,8.83E-02,1.05E+00,2.09E-01,1.56E-02},
      {2.70E+01,2.22E+00,8.56E-02,1.00E+00,2.02E-01,1.51E-02},
      {2.80E+01,2.17E+00,8.30E-02,9.62E-01,1.95E-01,1.46E-02},
      {2.90E+01,2.13E+00,8.06E-02,9.23E-01,1.89E-01,1.41E-02},
      {3.00E+01,2.08E+00,7.83E-02,8.87E-01,1.83E-01,1.37E-02},
      {3.10E+01,2.04E+00,7.62E-02,8.52E-01,1.77E-01,1.33E-02},
      {3.20E+01,2.00E+00,7.41E-02,8.20E-01,1.72E-01,1.29E-02},
      {3.30E+01,1.97E+00,7.22E-02,7.89E-01,1.66E-01,1.25E-02},
      {3.40E+01,1.93E+00,7.03E-02,7.61E-01,1.62E-01,1.22E-02},
      {3.50E+01,1.90E+00,6.85E-02,7.33E-01,1.57E-01,1.19E-02},
      {3.60E+01,1.87E+00,6.68E-02,7.07E-01,1.52E-01,1.15E-02},
      {3.70E+01,1.84E+00,6.51E-02,6.82E-01,1.48E-01,1.12E-02},
      {3.80E+01,1.81E+00,6.35E-02,6.59E-01,1.44E-01,1.09E-02},
      {3.90E+01,1.79E+00,6.20E-02,6.36E-01,1.40E-01,1.07E-02},
      {4.00E+01,1.76E+00,6.05E-02,6.15E-01,1.37E-01,1.04E-02},
      {4.10E+01,1.74E+00,5.91E-02,5.94E-01,1.33E-01,1.02E-02},
      {4.20E+01,1.71E+00,5.77E-02,5.74E-01,1.30E-01,9.91E-03},
      {4.30E+01,1.69E+00,5.64E-02,5.55E-01,1.26E-01,9.68E-03},
      {4.40E+01,1.67E+00,5.50E-02,5.37E-01,1.23E-01,9.45E-03},
      {4.50E+01,1.65E+00,5.38E-02,5.20E-01,1.20E-01,9.23E-03},
      {4.60E+01,1.63E+00,5.25E-02,5.03E-01,1.17E-01,9.03E-03},
      {4.70E+01,1.61E+00,5.13E-02,4.87E-01,1.14E-01,8.82E-03},
      {4.80E+01,1.59E+00,5.01E-02,4.71E-01,1.12E-01,8.63E-03},
      {4.90E+01,1.57E+00,4.90E-02,4.56E-01,1.09E-01,8.44E-03},
      {5.00E+01,1.56E+00,4.78E-02,4.42E-01,1.06E-01,8.26E-03},
      {5.10E+01,1.54E+00,4.67E-02,4.28E-01,1.04E-01,8.09E-03},
      {5.20E+01,1.52E+00,4.57E-02,4.14E-01,1.01E-01,7.92E-03},
      {5.30E+01,1.51E+00,4.46E-02,4.01E-01,9.91E-02,7.76E-03},
      {5.40E+01,1.49E+00,4.36E-02,3.89E-01,9.69E-02,7.60E-03},
      {5.50E+01,1.48E+00,4.26E-02,3.77E-01,9.47E-02,7.44E-03},
      {5.60E+01,1.46E+00,4.16E-02,3.65E-01,9.26E-02,7.30E-03},
      {5.70E+01,1.45E+00,4.06E-02,3.53E-01,9.06E-02,7.15E-03},
      {5.80E+01,1.44E+00,3.97E-02,3.42E-01,8.86E-02,7.01E-03},
      {5.90E+01,1.43E+00,3.87E-02,3.31E-01,8.67E-02,6.88E-03},
      {6.00E+01,1.41E+00,3.78E-02,3.21E-01,8.48E-02,6.75E-03},
      {6.10E+01,1.40E+00,3.69E-02,3.11E-01,8.30E-02,6.62E-03},
      {6.20E+01,1.39E+00,3.60E-02,3.01E-01,8.13E-02,6.50E-03},
      {6.30E+01,1.38E+00,3.52E-02,2.91E-01,7.96E-02,6.38E-03},
      {6.40E+01,1.37E+00,3.43E-02,2.82E-01,7.79E-02,6.26E-03},
      {6.50E+01,1.36E+00,3.35E-02,2.73E-01,7.63E-02,6.15E-03},
      {6.60E+01,1.35E+00,3.27E-02,2.64E-01,7.47E-02,6.04E-03},
      {6.70E+01,1.34E+00,3.19E-02,2.56E-01,7.32E-02,5.94E-03},
      {6.80E+01,1.33E+00,3.11E-02,2.48E-01,7.18E-02,5.83E-03},
      {6.90E+01,1.32E+00,3.04E-02,2.40E-01,7.03E-02,5.73E-03},
      {7.00E+01,1.31E+00,2.96E-02,2.32E-01,6.89E-02,5.63E-03},
      {7.10E+01,1.30E+00,2.89E-02,2.24E-01,6.76E-02,5.54E-03},
      {7.20E+01,1.29E+00,2.81E-02,2.17E-01,6.63E-02,5.45E-03},
      {7.30E+01,1.28E+00,2.74E-02,2.10E-01,6.50E-02,5.36E-03},
      {7.40E+01,1.27E+00,2.68E-02,2.03E-01,6.37E-02,5.27E-03},
      {7.50E+01,1.26E+00,2.61E-02,1.96E-01,6.25E-02,5.18E-03},
      {7.60E+01,1.26E+00,2.54E-02,1.89E-01,6.13E-02,5.10E-03},
      {7.70E+01,1.25E+00,2.48E-02,1.83E-01,6.02E-02,5.02E-03},
      {7.80E+01,1.24E+00,2.41E-02,1.77E-01,5.91E-02,4.94E-03},
      {7.90E+01,1.23E+00,2.35E-02,1.71E-01,5.80E-02,4.86E-03},
      {8.00E+01,1.23E+00,2.29E-02,1.65E-01,5.69E-02,4.79E-03},
      {8.10E+01,1.22E+00,2.23E-02,1.59E-01,5.59E-02,4.72E-03},
      {8.20E+01,1.21E+00,2.17E-02,1.53E-01,5.48E-02,4.64E-03},
      {8.30E+01,1.21E+00,2.11E-02,1.48E-01,5.39E-02,4.58E-03},
      {8.40E+01,1.20E+00,2.06E-02,1.42E-01,5.29E-02,4.51E-03},
      {8.50E+01,1.19E+00,2.00E-02,1.37E-01,5.19E-02,4.44E-03},
      {8.60E+01,1.19E+00,1.95E-02,1.32E-01,5.10E-02,4.38E-03},
      {8.70E+01,1.18E+00,1.90E-02,1.27E-01,5.01E-02,4.31E-03},
      {8.80E+01,1.18E+00,1.84E-02,1.22E-01,4.93E-02,4.25E-03},
      {8.90E+01,1.17E+00,1.79E-02,1.18E-01,4.84E-02,4.19E-03},
      {9.00E+01,1.17E+00,1.74E-02,1.13E-01,4.76E-02,4.13E-03},
      {9.10E+01,1.16E+00,1.70E-02,1.09E-01,4.68E-02,4.08E-03},
      {9.20E+01,1.15E+00,1.65E-02,1.05E-01,4.60E-02,4.02E-03},
      {9.30E+01,1.15E+00,1.60E-02,1.00E-01,4.52E-02,3.97E-03},
      {9.40E+01,1.14E+00,1.56E-02,9.64E-02,4.44E-02,3.91E-03},
      {9.50E+01,1.14E+00,1.52E-02,9.24E-02,4.37E-02,3.86E-03},
      {9.60E+01,1.14E+00,1.47E-02,8.86E-02,4.30E-02,3.81E-03},
      {9.70E+01,1.13E+00,1.43E-02,8.48E-02,4.23E-02,3.76E-03},
      {9.80E+01,1.13E+00,1.39E-02,8.12E-02,4.16E-02,3.71E-03},
      {9.90E+01,1.12E+00,1.35E-02,7.77E-02,4.09E-02,3.66E-03},
      {1.00E+02,1.12E+00,1.31E-02,7.42E-02,4.03E-02,3.61E-03},
      {1.01E+02,1.11E+00,1.28E-02,7.09E-02,3.96E-02,3.57E-03},
      {1.02E+02,1.11E+00,1.24E-02,6.77E-02,3.90E-02,3.52E-03},
      {1.03E+02,1.11E+00,1.20E-02,6.45E-02,3.84E-02,3.48E-03},
      {1.04E+02,1.10E+00,1.17E-02,6.14E-02,3.78E-02,3.44E-03},
      {1.05E+02,1.10E+00,1.13E-02,5.84E-02,3.72E-02,3.40E-03},
      {1.06E+02,1.10E+00,1.10E-02,5.56E-02,3.66E-02,3.35E-03},
      {1.07E+02,1.09E+00,1.07E-02,5.27E-02,3.60E-02,3.31E-03},
      {1.08E+02,1.09E+00,1.04E-02,5.00E-02,3.55E-02,3.27E-03},
      {1.09E+02,1.09E+00,1.01E-02,4.74E-02,3.49E-02,3.24E-03},
      {1.10E+02,1.08E+00,9.76E-03,4.48E-02,3.44E-02,3.20E-03},
      {1.11E+02,1.08E+00,9.47E-03,4.23E-02,3.39E-02,3.16E-03},
      {1.12E+02,1.08E+00,9.19E-03,3.98E-02,3.34E-02,3.12E-03},
      {1.13E+02,1.07E+00,8.92E-03,3.75E-02,3.29E-02,3.09E-03},
      {1.14E+02,1.07E+00,8.65E-03,3.52E-02,3.24E-02,3.05E-03},
      {1.15E+02,1.07E+00,8.39E-03,3.30E-02,3.19E-02,3.02E-03},
      {1.16E+02,1.07E+00,8.14E-03,3.08E-02,3.15E-02,2.99E-03},
      {1.17E+02,1.06E+00,7.89E-03,2.87E-02,3.10E-02,2.95E-03},
      {1.18E+02,1.06E+00,7.65E-03,2.67E-02,3.06E-02,2.92E-03},
      {1.19E+02,1.06E+00,7.42E-03,2.48E-02,3.01E-02,2.89E-03},
      {1.20E+02,1.06E+00,7.19E-03,2.29E-02,2.97E-02,2.86E-03},
      {1.21E+02,1.05E+00,6.98E-03,2.10E-02,2.93E-02,2.83E-03},
      {1.22E+02,1.05E+00,6.76E-03,1.92E-02,2.89E-02,2.80E-03},
      {1.23E+02,1.05E+00,6.56E-03,1.75E-02,2.85E-02,2.77E-03},
      {1.24E+02,1.05E+00,6.36E-03,1.58E-02,2.81E-02,2.74E-03},
      {1.25E+02,1.04E+00,6.16E-03,1.42E-02,2.77E-02,2.71E-03},
      {1.26E+02,1.04E+00,5.97E-03,1.27E-02,2.73E-02,2.68E-03},
      {1.27E+02,1.04E+00,5.79E-03,1.11E-02,2.69E-02,2.66E-03},
      {1.28E+02,1.04E+00,5.61E-03,9.68E-03,2.66E-02,2.63E-03},
      {1.29E+02,1.04E+00,5.44E-03,8.27E-03,2.62E-02,2.60E-03},
      {1.30E+02,1.04E+00,5.28E-03,6.90E-03,2.59E-02,2.58E-03},
      {1.31E+02,1.03E+00,5.12E-03,5.58E-03,2.55E-02,2.55E-03},
      {1.32E+02,1.03E+00,4.96E-03,4.30E-03,2.52E-02,2.53E-03},
      {1.33E+02,1.03E+00,4.81E-03,3.06E-03,2.48E-02,2.50E-03},
      {1.34E+02,1.03E+00,4.66E-03,1.87E-03,2.45E-02,2.48E-03},
      {1.35E+02,1.03E+00,4.52E-03,7.22E-04,2.42E-02,2.46E-03},
      {1.36E+02,1.03E+00,4.39E-03,-3.88E-04,2.39E-02,2.43E-03},
      {1.37E+02,1.02E+00,4.25E-03,-1.46E-03,2.36E-02,2.41E-03},
      {1.38E+02,1.02E+00,4.13E-03,-2.49E-03,2.33E-02,2.39E-03},
      {1.39E+02,1.02E+00,4.01E-03,-3.48E-03,2.30E-02,2.36E-03},
      {1.40E+02,1.02E+00,3.89E-03,-4.44E-03,2.27E-02,2.34E-03},
      {1.41E+02,1.02E+00,3.77E-03,-5.37E-03,2.24E-02,2.32E-03},
      {1.42E+02,1.02E+00,3.66E-03,-6.25E-03,2.21E-02,2.30E-03},
      {1.43E+02,1.02E+00,3.56E-03,-7.11E-03,2.18E-02,2.28E-03},
      {1.44E+02,1.02E+00,3.46E-03,-7.93E-03,2.16E-02,2.26E-03},
      {1.45E+02,1.01E+00,3.36E-03,-8.72E-03,2.13E-02,2.24E-03},
      {1.46E+02,1.01E+00,3.26E-03,-9.48E-03,2.10E-02,2.22E-03},
      {1.47E+02,1.01E+00,3.17E-03,-1.02E-02,2.08E-02,2.20E-03},
      {1.48E+02,1.01E+00,3.09E-03,-1.09E-02,2.05E-02,2.18E-03},
      {1.49E+02,1.01E+00,3.00E-03,-1.16E-02,2.03E-02,2.16E-03},
      {1.50E+02,1.01E+00,2.93E-03,-1.22E-02,2.01E-02,2.14E-03},
      {1.51E+02,1.01E+00,2.85E-03,-1.29E-02,1.98E-02,2.13E-03},
      {1.52E+02,1.01E+00,2.78E-03,-1.34E-02,1.96E-02,2.11E-03},
      {1.53E+02,1.01E+00,2.71E-03,-1.40E-02,1.93E-02,2.09E-03},
      {1.54E+02,1.01E+00,2.64E-03,-1.46E-02,1.91E-02,2.07E-03},
      {1.55E+02,1.01E+00,2.58E-03,-1.51E-02,1.89E-02,2.06E-03},
      {1.56E+02,1.01E+00,2.52E-03,-1.56E-02,1.87E-02,2.04E-03},
      {1.57E+02,1.00E+00,2.46E-03,-1.60E-02,1.85E-02,2.02E-03},
      {1.58E+02,1.00E+00,2.41E-03,-1.65E-02,1.83E-02,2.01E-03},
      {1.59E+02,1.00E+00,2.36E-03,-1.69E-02,1.80E-02,1.99E-03},
      {1.60E+02,1.00E+00,2.31E-03,-1.73E-02,1.78E-02,1.98E-03},
      {1.61E+02,1.00E+00,2.26E-03,-1.77E-02,1.76E-02,1.96E-03},
      {1.62E+02,1.00E+00,2.22E-03,-1.81E-02,1.74E-02,1.95E-03},
      {1.63E+02,1.00E+00,2.18E-03,-1.85E-02,1.72E-02,1.93E-03},
      {1.64E+02,1.00E+00,2.14E-03,-1.88E-02,1.71E-02,1.92E-03},
      {1.65E+02,1.00E+00,2.10E-03,-1.91E-02,1.69E-02,1.90E-03},
      {1.66E+02,9.99E-01,2.06E-03,-1.94E-02,1.67E-02,1.89E-03},
      {1.67E+02,9.99E-01,2.03E-03,-1.97E-02,1.65E-02,1.87E-03},
      {1.68E+02,9.98E-01,2.00E-03,-2.00E-02,1.63E-02,1.86E-03},
      {1.69E+02,9.98E-01,1.97E-03,-2.02E-02,1.61E-02,1.84E-03},
      {1.70E+02,9.97E-01,1.94E-03,-2.05E-02,1.60E-02,1.83E-03},
      {1.71E+02,9.97E-01,1.91E-03,-2.07E-02,1.58E-02,1.82E-03},
      {1.72E+02,9.97E-01,1.89E-03,-2.09E-02,1.56E-02,1.80E-03},
      {1.73E+02,9.96E-01,1.86E-03,-2.11E-02,1.55E-02,1.79E-03},
      {1.74E+02,9.96E-01,1.84E-03,-2.13E-02,1.53E-02,1.78E-03},
      {1.75E+02,9.95E-01,1.81E-03,-2.15E-02,1.51E-02,1.77E-03},
      {1.76E+02,9.95E-01,1.79E-03,-2.16E-02,1.50E-02,1.75E-03},
      {1.77E+02,9.95E-01,1.77E-03,-2.18E-02,1.48E-02,1.74E-03},
      {1.78E+02,9.94E-01,1.75E-03,-2.19E-02,1.47E-02,1.73E-03},
      {1.79E+02,9.94E-01,1.73E-03,-2.21E-02,1.45E-02,1.72E-03},
      {1.80E+02,9.94E-01,1.71E-03,-2.22E-02,1.44E-02,1.71E-03},
      {1.81E+02,9.94E-01,1.70E-03,-2.23E-02,1.42E-02,1.69E-03},
      {1.82E+02,9.93E-01,1.68E-03,-2.24E-02,1.41E-02,1.68E-03},
      {1.83E+02,9.93E-01,1.66E-03,-2.25E-02,1.39E-02,1.67E-03},
      {1.84E+02,9.93E-01,1.64E-03,-2.26E-02,1.38E-02,1.66E-03},
      {1.85E+02,9.93E-01,1.63E-03,-2.26E-02,1.37E-02,1.65E-03},
      {1.86E+02,9.92E-01,1.61E-03,-2.27E-02,1.35E-02,1.64E-03},
      {1.87E+02,9.92E-01,1.60E-03,-2.27E-02,1.34E-02,1.63E-03},
      {1.88E+02,9.92E-01,1.58E-03,-2.28E-02,1.33E-02,1.62E-03},
      {1.89E+02,9.92E-01,1.57E-03,-2.28E-02,1.31E-02,1.61E-03},
      {1.90E+02,9.92E-01,1.55E-03,-2.29E-02,1.30E-02,1.60E-03},
      {1.91E+02,9.92E-01,1.54E-03,-2.29E-02,1.29E-02,1.59E-03},
      {1.92E+02,9.91E-01,1.52E-03,-2.29E-02,1.28E-02,1.58E-03},
      {1.93E+02,9.91E-01,1.51E-03,-2.29E-02,1.26E-02,1.57E-03},
      {1.94E+02,9.91E-01,1.49E-03,-2.29E-02,1.25E-02,1.56E-03},
      {1.95E+02,9.91E-01,1.48E-03,-2.29E-02,1.24E-02,1.55E-03},
      {1.96E+02,9.91E-01,1.47E-03,-2.29E-02,1.23E-02,1.54E-03},
      {1.97E+02,9.91E-01,1.45E-03,-2.29E-02,1.22E-02,1.53E-03},
      {1.98E+02,9.91E-01,1.44E-03,-2.29E-02,1.21E-02,1.52E-03},
      {1.99E+02,9.91E-01,1.43E-03,-2.29E-02,1.20E-02,1.51E-03},
      {2.00E+02,9.90E-01,1.41E-03,-2.29E-02,1.18E-02,1.50E-03},
      {2.01E+02,9.90E-01,1.40E-03,-2.28E-02,1.17E-02,1.49E-03},
      {2.02E+02,9.90E-01,1.39E-03,-2.28E-02,1.16E-02,1.48E-03},
      {2.03E+02,9.90E-01,1.38E-03,-2.28E-02,1.15E-02,1.48E-03},
      {2.04E+02,9.90E-01,1.36E-03,-2.27E-02,1.14E-02,1.47E-03},
      {2.05E+02,9.90E-01,1.35E-03,-2.27E-02,1.13E-02,1.46E-03},
      {2.06E+02,9.90E-01,1.34E-03,-2.26E-02,1.12E-02,1.45E-03},
      {2.07E+02,9.90E-01,1.33E-03,-2.26E-02,1.11E-02,1.44E-03},
      {2.08E+02,9.90E-01,1.31E-03,-2.25E-02,1.10E-02,1.43E-03},
      {2.09E+02,9.90E-01,1.30E-03,-2.25E-02,1.09E-02,1.43E-03},
      {2.10E+02,9.90E-01,1.29E-03,-2.24E-02,1.08E-02,1.42E-03},
      {2.11E+02,9.90E-01,1.28E-03,-2.23E-02,1.07E-02,1.41E-03},
      {2.12E+02,9.90E-01,1.27E-03,-2.23E-02,1.06E-02,1.40E-03},
      {2.13E+02,9.90E-01,1.25E-03,-2.22E-02,1.05E-02,1.39E-03},
      {2.14E+02,9.90E-01,1.24E-03,-2.21E-02,1.05E-02,1.39E-03},
      {2.15E+02,9.90E-01,1.23E-03,-2.21E-02,1.04E-02,1.38E-03},
      {2.16E+02,9.90E-01,1.22E-03,-2.20E-02,1.03E-02,1.37E-03},
      {2.17E+02,9.90E-01,1.21E-03,-2.19E-02,1.02E-02,1.37E-03},
      {2.18E+02,9.90E-01,1.20E-03,-2.18E-02,1.01E-02,1.36E-03},
      {2.19E+02,9.90E-01,1.18E-03,-2.17E-02,1.00E-02,1.35E-03},
      {2.20E+02,9.90E-01,1.17E-03,-2.17E-02,9.94E-03,1.34E-03},
      {2.21E+02,9.90E-01,1.16E-03,-2.16E-02,9.86E-03,1.34E-03},
      {2.22E+02,9.90E-01,1.15E-03,-2.15E-02,9.77E-03,1.33E-03},
      {2.23E+02,9.90E-01,1.14E-03,-2.14E-02,9.69E-03,1.32E-03},
      {2.24E+02,9.90E-01,1.13E-03,-2.13E-02,9.61E-03,1.32E-03},
      {2.25E+02,9.90E-01,1.12E-03,-2.12E-02,9.54E-03,1.31E-03},
      {2.26E+02,9.90E-01,1.11E-03,-2.11E-02,9.46E-03,1.30E-03},
      {2.27E+02,9.90E-01,1.10E-03,-2.10E-02,9.38E-03,1.30E-03},
      {2.28E+02,9.90E-01,1.09E-03,-2.09E-02,9.31E-03,1.29E-03},
      {2.29E+02,9.90E-01,1.07E-03,-2.08E-02,9.23E-03,1.28E-03},
      {2.30E+02,9.90E-01,1.06E-03,-2.07E-02,9.16E-03,1.28E-03},
      {2.31E+02,9.90E-01,1.05E-03,-2.06E-02,9.09E-03,1.27E-03},
      {2.32E+02,9.90E-01,1.04E-03,-2.05E-02,9.02E-03,1.27E-03},
      {2.33E+02,9.90E-01,1.03E-03,-2.04E-02,8.95E-03,1.26E-03},
      {2.34E+02,9.90E-01,1.02E-03,-2.03E-02,8.88E-03,1.25E-03},
      {2.35E+02,9.90E-01,1.01E-03,-2.02E-02,8.81E-03,1.25E-03},
      {2.36E+02,9.90E-01,1.00E-03,-2.01E-02,8.74E-03,1.24E-03},
      {2.37E+02,9.90E-01,9.92E-04,-2.00E-02,8.68E-03,1.24E-03},
      {2.38E+02,9.90E-01,9.83E-04,-1.99E-02,8.61E-03,1.23E-03},
      {2.39E+02,9.90E-01,9.73E-04,-1.98E-02,8.54E-03,1.22E-03},
      {2.40E+02,9.90E-01,9.63E-04,-1.97E-02,8.48E-03,1.22E-03},
      {2.41E+02,9.90E-01,9.53E-04,-1.96E-02,8.42E-03,1.21E-03},
      {2.42E+02,9.90E-01,9.44E-04,-1.95E-02,8.35E-03,1.21E-03},
      {2.43E+02,9.90E-01,9.34E-04,-1.94E-02,8.29E-03,1.20E-03},
      {2.44E+02,9.90E-01,9.25E-04,-1.93E-02,8.23E-03,1.20E-03},
      {2.45E+02,9.90E-01,9.15E-04,-1.92E-02,8.17E-03,1.19E-03},
      {2.46E+02,9.90E-01,9.06E-04,-1.91E-02,8.11E-03,1.19E-03},
      {2.47E+02,9.90E-01,8.97E-04,-1.90E-02,8.05E-03,1.18E-03},
      {2.48E+02,9.90E-01,8.87E-04,-1.89E-02,7.99E-03,1.18E-03},
      {2.49E+02,9.90E-01,8.78E-04,-1.88E-02,7.93E-03,1.17E-03},
      {2.50E+02,9.90E-01,8.69E-04,-1.87E-02,7.88E-03,1.17E-03},
      {2.51E+02,9.90E-01,8.60E-04,-1.86E-02,7.82E-03,1.16E-03},
      {2.52E+02,9.90E-01,8.52E-04,-1.85E-02,7.76E-03,1.16E-03},
      {2.53E+02,9.90E-01,8.43E-04,-1.84E-02,7.71E-03,1.15E-03},
      {2.54E+02,9.91E-01,8.34E-04,-1.83E-02,7.65E-03,1.15E-03},
      {2.55E+02,9.91E-01,8.26E-04,-1.82E-02,7.60E-03,1.14E-03},
      {2.56E+02,9.91E-01,8.17E-04,-1.81E-02,7.55E-03,1.14E-03},
      {2.57E+02,9.91E-01,8.09E-04,-1.80E-02,7.49E-03,1.13E-03},
      {2.58E+02,9.91E-01,8.00E-04,-1.79E-02,7.44E-03,1.13E-03},
      {2.59E+02,9.91E-01,7.92E-04,-1.78E-02,7.39E-03,1.12E-03},
      {2.60E+02,9.91E-01,7.84E-04,-1.77E-02,7.34E-03,1.12E-03},
      {2.61E+02,9.91E-01,7.76E-04,-1.75E-02,7.29E-03,1.11E-03},
      {2.62E+02,9.91E-01,7.68E-04,-1.74E-02,7.24E-03,1.11E-03},
      {2.63E+02,9.91E-01,7.60E-04,-1.73E-02,7.19E-03,1.10E-03},
      {2.64E+02,9.91E-01,7.52E-04,-1.72E-02,7.14E-03,1.10E-03},
      {2.65E+02,9.91E-01,7.44E-04,-1.71E-02,7.09E-03,1.10E-03},
      {2.66E+02,9.91E-01,7.36E-04,-1.70E-02,7.04E-03,1.09E-03},
      {2.67E+02,9.91E-01,7.28E-04,-1.69E-02,7.00E-03,1.09E-03},
      {2.68E+02,9.91E-01,7.21E-04,-1.68E-02,6.95E-03,1.08E-03},
      {2.69E+02,9.91E-01,7.13E-04,-1.67E-02,6.91E-03,1.08E-03},
      {2.70E+02,9.91E-01,7.06E-04,-1.66E-02,6.86E-03,1.07E-03},
      {2.71E+02,9.91E-01,6.98E-04,-1.65E-02,6.81E-03,1.07E-03},
      {2.72E+02,9.91E-01,6.91E-04,-1.64E-02,6.77E-03,1.07E-03},
      {2.73E+02,9.91E-01,6.84E-04,-1.63E-02,6.73E-03,1.06E-03},
      {2.74E+02,9.92E-01,6.77E-04,-1.62E-02,6.68E-03,1.06E-03},
      {2.75E+02,9.92E-01,6.69E-04,-1.61E-02,6.64E-03,1.05E-03},
      {2.76E+02,9.92E-01,6.62E-04,-1.60E-02,6.60E-03,1.05E-03},
      {2.77E+02,9.92E-01,6.55E-04,-1.59E-02,6.55E-03,1.05E-03},
      {2.78E+02,9.92E-01,6.48E-04,-1.58E-02,6.51E-03,1.04E-03},
      {2.79E+02,9.92E-01,6.42E-04,-1.57E-02,6.47E-03,1.04E-03},
      {2.80E+02,9.92E-01,6.35E-04,-1.57E-02,6.43E-03,1.03E-03},
      {2.81E+02,9.92E-01,6.28E-04,-1.56E-02,6.39E-03,1.03E-03},
      {2.82E+02,9.92E-01,6.22E-04,-1.55E-02,6.35E-03,1.03E-03},
      {2.83E+02,9.92E-01,6.15E-04,-1.54E-02,6.31E-03,1.02E-03},
      {2.84E+02,9.92E-01,6.08E-04,-1.53E-02,6.27E-03,1.02E-03},
      {2.85E+02,9.92E-01,6.02E-04,-1.52E-02,6.23E-03,1.02E-03},
      {2.86E+02,9.92E-01,5.96E-04,-1.51E-02,6.19E-03,1.01E-03},
      {2.87E+02,9.92E-01,5.89E-04,-1.50E-02,6.15E-03,1.01E-03},
      {2.88E+02,9.92E-01,5.83E-04,-1.49E-02,6.12E-03,1.01E-03},
      {2.89E+02,9.92E-01,5.77E-04,-1.48E-02,6.08E-03,1.00E-03},
      {2.90E+02,9.92E-01,5.71E-04,-1.47E-02,6.04E-03,9.98E-04},
      {2.91E+02,9.92E-01,5.65E-04,-1.46E-02,6.01E-03,9.94E-04},
      {2.92E+02,9.92E-01,5.59E-04,-1.45E-02,5.97E-03,9.91E-04},
      {2.93E+02,9.92E-01,5.53E-04,-1.44E-02,5.93E-03,9.88E-04},
      {2.94E+02,9.93E-01,5.47E-04,-1.44E-02,5.90E-03,9.84E-04},
      {2.95E+02,9.93E-01,5.41E-04,-1.43E-02,5.86E-03,9.81E-04},
      {2.96E+02,9.93E-01,5.35E-04,-1.42E-02,5.83E-03,9.77E-04},
      {2.97E+02,9.93E-01,5.30E-04,-1.41E-02,5.79E-03,9.74E-04},
      {2.98E+02,9.93E-01,5.24E-04,-1.40E-02,5.76E-03,9.71E-04},
      {2.99E+02,9.93E-01,5.19E-04,-1.39E-02,5.73E-03,9.67E-04},
      {3.00E+02,9.93E-01,5.13E-04,-1.38E-02,5.69E-03,9.64E-04}
    };

  std::vector<std::vector<float>> Yuki_1_07 =
    {
      {1.00E+00,4.28E+01,1.99E+00,3.69E+01,4.56E+00,3.29E-01},
      {2.00E+00,2.13E+01,9.92E-01,1.78E+01,2.28E+00,1.64E-01},
      {3.00E+00,1.42E+01,6.58E-01,1.15E+01,1.51E+00,1.09E-01},
      {4.00E+00,1.06E+01,4.90E-01,8.40E+00,1.13E+00,8.18E-02},
      {5.00E+00,8.50E+00,3.89E-01,6.53E+00,9.05E-01,6.53E-02},
      {6.00E+00,7.11E+00,3.23E-01,5.30E+00,7.54E-01,5.45E-02},
      {7.00E+00,6.13E+00,2.76E-01,4.44E+00,6.48E-01,4.68E-02},
      {8.00E+00,5.42E+00,2.42E-01,3.81E+00,5.69E-01,4.11E-02},
      {9.00E+00,4.87E+00,2.15E-01,3.33E+00,5.08E-01,3.67E-02},
      {1.00E+01,4.45E+00,1.95E-01,2.95E+00,4.60E-01,3.33E-02},
      {1.10E+01,4.10E+00,1.78E-01,2.65E+00,4.21E-01,3.05E-02},
      {1.20E+01,3.82E+00,1.64E-01,2.40E+00,3.89E-01,2.81E-02},
      {1.30E+01,3.58E+00,1.53E-01,2.20E+00,3.62E-01,2.62E-02},
      {1.40E+01,3.38E+00,1.43E-01,2.02E+00,3.38E-01,2.45E-02},
      {1.50E+01,3.21E+00,1.35E-01,1.87E+00,3.18E-01,2.31E-02},
      {1.60E+01,3.06E+00,1.27E-01,1.74E+00,3.00E-01,2.18E-02},
      {1.70E+01,2.93E+00,1.21E-01,1.63E+00,2.85E-01,2.07E-02},
      {1.80E+01,2.82E+00,1.15E-01,1.53E+00,2.70E-01,1.97E-02},
      {1.90E+01,2.72E+00,1.10E-01,1.44E+00,2.58E-01,1.88E-02},
      {2.00E+01,2.62E+00,1.05E-01,1.36E+00,2.46E-01,1.79E-02},
      {2.10E+01,2.54E+00,1.01E-01,1.29E+00,2.36E-01,1.72E-02},
      {2.20E+01,2.46E+00,9.73E-02,1.22E+00,2.26E-01,1.65E-02},
      {2.30E+01,2.39E+00,9.38E-02,1.16E+00,2.17E-01,1.59E-02},
      {2.40E+01,2.33E+00,9.05E-02,1.10E+00,2.09E-01,1.53E-02},
      {2.50E+01,2.27E+00,8.75E-02,1.05E+00,2.02E-01,1.48E-02},
      {2.60E+01,2.22E+00,8.47E-02,1.01E+00,1.94E-01,1.42E-02},
      {2.70E+01,2.16E+00,8.21E-02,9.63E-01,1.88E-01,1.38E-02},
      {2.80E+01,2.12E+00,7.97E-02,9.22E-01,1.81E-01,1.33E-02},
      {2.90E+01,2.07E+00,7.74E-02,8.84E-01,1.76E-01,1.29E-02},
      {3.00E+01,2.03E+00,7.52E-02,8.49E-01,1.70E-01,1.25E-02},
      {3.10E+01,1.99E+00,7.31E-02,8.15E-01,1.65E-01,1.21E-02},
      {3.20E+01,1.96E+00,7.11E-02,7.84E-01,1.60E-01,1.18E-02},
      {3.30E+01,1.92E+00,6.92E-02,7.54E-01,1.55E-01,1.15E-02},
      {3.40E+01,1.89E+00,6.74E-02,7.26E-01,1.50E-01,1.11E-02},
      {3.50E+01,1.86E+00,6.57E-02,6.99E-01,1.46E-01,1.08E-02},
      {3.60E+01,1.83E+00,6.41E-02,6.74E-01,1.42E-01,1.06E-02},
      {3.70E+01,1.80E+00,6.25E-02,6.50E-01,1.38E-01,1.03E-02},
      {3.80E+01,1.77E+00,6.09E-02,6.27E-01,1.34E-01,1.00E-02},
      {3.90E+01,1.75E+00,5.94E-02,6.05E-01,1.31E-01,9.76E-03},
      {4.00E+01,1.72E+00,5.80E-02,5.84E-01,1.27E-01,9.52E-03},
      {4.10E+01,1.70E+00,5.66E-02,5.64E-01,1.24E-01,9.29E-03},
      {4.20E+01,1.67E+00,5.53E-02,5.45E-01,1.21E-01,9.06E-03},
      {4.30E+01,1.65E+00,5.39E-02,5.26E-01,1.17E-01,8.85E-03},
      {4.40E+01,1.63E+00,5.27E-02,5.09E-01,1.15E-01,8.64E-03},
      {4.50E+01,1.61E+00,5.14E-02,4.92E-01,1.12E-01,8.45E-03},
      {4.60E+01,1.59E+00,5.02E-02,4.75E-01,1.09E-01,8.25E-03},
      {4.70E+01,1.57E+00,4.90E-02,4.60E-01,1.06E-01,8.07E-03},
      {4.80E+01,1.56E+00,4.79E-02,4.45E-01,1.04E-01,7.89E-03},
      {4.90E+01,1.54E+00,4.67E-02,4.30E-01,1.01E-01,7.72E-03},
      {5.00E+01,1.52E+00,4.56E-02,4.16E-01,9.89E-02,7.56E-03},
      {5.10E+01,1.51E+00,4.46E-02,4.02E-01,9.66E-02,7.40E-03},
      {5.20E+01,1.49E+00,4.35E-02,3.89E-01,9.43E-02,7.24E-03},
      {5.30E+01,1.48E+00,4.25E-02,3.77E-01,9.22E-02,7.09E-03},
      {5.40E+01,1.46E+00,4.15E-02,3.65E-01,9.01E-02,6.95E-03},
      {5.50E+01,1.45E+00,4.05E-02,3.53E-01,8.81E-02,6.81E-03},
      {5.60E+01,1.43E+00,3.95E-02,3.41E-01,8.61E-02,6.67E-03},
      {5.70E+01,1.42E+00,3.86E-02,3.30E-01,8.42E-02,6.54E-03},
      {5.80E+01,1.41E+00,3.76E-02,3.20E-01,8.24E-02,6.42E-03},
      {5.90E+01,1.40E+00,3.67E-02,3.09E-01,8.06E-02,6.29E-03},
      {6.00E+01,1.38E+00,3.58E-02,2.99E-01,7.88E-02,6.17E-03},
      {6.10E+01,1.37E+00,3.49E-02,2.89E-01,7.72E-02,6.06E-03},
      {6.20E+01,1.36E+00,3.41E-02,2.80E-01,7.55E-02,5.94E-03},
      {6.30E+01,1.35E+00,3.32E-02,2.71E-01,7.39E-02,5.84E-03},
      {6.40E+01,1.34E+00,3.24E-02,2.62E-01,7.24E-02,5.73E-03},
      {6.50E+01,1.33E+00,3.16E-02,2.53E-01,7.09E-02,5.63E-03},
      {6.60E+01,1.32E+00,3.08E-02,2.45E-01,6.95E-02,5.53E-03},
      {6.70E+01,1.31E+00,3.00E-02,2.37E-01,6.81E-02,5.43E-03},
      {6.80E+01,1.30E+00,2.93E-02,2.29E-01,6.67E-02,5.34E-03},
      {6.90E+01,1.29E+00,2.85E-02,2.21E-01,6.53E-02,5.24E-03},
      {7.00E+01,1.28E+00,2.78E-02,2.14E-01,6.40E-02,5.15E-03},
      {7.10E+01,1.27E+00,2.71E-02,2.06E-01,6.28E-02,5.07E-03},
      {7.20E+01,1.27E+00,2.64E-02,1.99E-01,6.16E-02,4.98E-03},
      {7.30E+01,1.26E+00,2.57E-02,1.92E-01,6.04E-02,4.90E-03},
      {7.40E+01,1.25E+00,2.50E-02,1.86E-01,5.92E-02,4.82E-03},
      {7.50E+01,1.24E+00,2.44E-02,1.79E-01,5.81E-02,4.74E-03},
      {7.60E+01,1.23E+00,2.37E-02,1.73E-01,5.70E-02,4.67E-03},
      {7.70E+01,1.23E+00,2.31E-02,1.67E-01,5.59E-02,4.59E-03},
      {7.80E+01,1.22E+00,2.25E-02,1.61E-01,5.48E-02,4.52E-03},
      {7.90E+01,1.21E+00,2.19E-02,1.55E-01,5.38E-02,4.45E-03},
      {8.00E+01,1.21E+00,2.13E-02,1.50E-01,5.28E-02,4.38E-03},
      {8.10E+01,1.20E+00,2.07E-02,1.44E-01,5.19E-02,4.31E-03},
      {8.20E+01,1.19E+00,2.01E-02,1.39E-01,5.09E-02,4.25E-03},
      {8.30E+01,1.19E+00,1.96E-02,1.34E-01,5.00E-02,4.19E-03},
      {8.40E+01,1.18E+00,1.90E-02,1.29E-01,4.91E-02,4.12E-03},
      {8.50E+01,1.18E+00,1.85E-02,1.24E-01,4.82E-02,4.06E-03},
      {8.60E+01,1.17E+00,1.80E-02,1.19E-01,4.74E-02,4.00E-03},
      {8.70E+01,1.16E+00,1.75E-02,1.14E-01,4.65E-02,3.95E-03},
      {8.80E+01,1.16E+00,1.70E-02,1.10E-01,4.57E-02,3.89E-03},
      {8.90E+01,1.15E+00,1.65E-02,1.05E-01,4.49E-02,3.83E-03},
      {9.00E+01,1.15E+00,1.60E-02,1.01E-01,4.42E-02,3.78E-03},
      {9.10E+01,1.14E+00,1.56E-02,9.71E-02,4.34E-02,3.73E-03},
      {9.20E+01,1.14E+00,1.51E-02,9.31E-02,4.27E-02,3.68E-03},
      {9.30E+01,1.13E+00,1.47E-02,8.92E-02,4.19E-02,3.63E-03},
      {9.40E+01,1.13E+00,1.43E-02,8.55E-02,4.12E-02,3.58E-03},
      {9.50E+01,1.13E+00,1.39E-02,8.18E-02,4.05E-02,3.53E-03},
      {9.60E+01,1.12E+00,1.34E-02,7.82E-02,3.99E-02,3.48E-03},
      {9.70E+01,1.12E+00,1.31E-02,7.48E-02,3.92E-02,3.44E-03},
      {9.80E+01,1.11E+00,1.27E-02,7.14E-02,3.86E-02,3.39E-03},
      {9.90E+01,1.11E+00,1.23E-02,6.81E-02,3.79E-02,3.35E-03},
      {1.00E+02,1.11E+00,1.19E-02,6.50E-02,3.73E-02,3.31E-03},
      {1.01E+02,1.10E+00,1.16E-02,6.19E-02,3.67E-02,3.26E-03},
      {1.02E+02,1.10E+00,1.12E-02,5.89E-02,3.61E-02,3.22E-03},
      {1.03E+02,1.09E+00,1.09E-02,5.60E-02,3.56E-02,3.18E-03},
      {1.04E+02,1.09E+00,1.06E-02,5.32E-02,3.50E-02,3.14E-03},
      {1.05E+02,1.09E+00,1.02E-02,5.05E-02,3.45E-02,3.11E-03},
      {1.06E+02,1.08E+00,9.92E-03,4.78E-02,3.39E-02,3.07E-03},
      {1.07E+02,1.08E+00,9.62E-03,4.52E-02,3.34E-02,3.03E-03},
      {1.08E+02,1.08E+00,9.32E-03,4.27E-02,3.29E-02,3.00E-03},
      {1.09E+02,1.08E+00,9.04E-03,4.03E-02,3.24E-02,2.96E-03},
      {1.10E+02,1.07E+00,8.76E-03,3.80E-02,3.19E-02,2.93E-03},
      {1.11E+02,1.07E+00,8.49E-03,3.57E-02,3.14E-02,2.89E-03},
      {1.12E+02,1.07E+00,8.22E-03,3.35E-02,3.09E-02,2.86E-03},
      {1.13E+02,1.06E+00,7.97E-03,3.14E-02,3.05E-02,2.83E-03},
      {1.14E+02,1.06E+00,7.72E-03,2.93E-02,3.00E-02,2.79E-03},
      {1.15E+02,1.06E+00,7.48E-03,2.73E-02,2.96E-02,2.76E-03},
      {1.16E+02,1.06E+00,7.24E-03,2.54E-02,2.92E-02,2.73E-03},
      {1.17E+02,1.06E+00,7.02E-03,2.35E-02,2.87E-02,2.70E-03},
      {1.18E+02,1.05E+00,6.79E-03,2.17E-02,2.83E-02,2.67E-03},
      {1.19E+02,1.05E+00,6.58E-03,1.99E-02,2.79E-02,2.64E-03},
      {1.20E+02,1.05E+00,6.37E-03,1.83E-02,2.75E-02,2.62E-03},
      {1.21E+02,1.05E+00,6.17E-03,1.66E-02,2.71E-02,2.59E-03},
      {1.22E+02,1.04E+00,5.98E-03,1.50E-02,2.67E-02,2.56E-03},
      {1.23E+02,1.04E+00,5.79E-03,1.35E-02,2.64E-02,2.53E-03},
      {1.24E+02,1.04E+00,5.60E-03,1.20E-02,2.60E-02,2.51E-03},
      {1.25E+02,1.04E+00,5.43E-03,1.06E-02,2.56E-02,2.48E-03},
      {1.26E+02,1.04E+00,5.25E-03,9.20E-03,2.53E-02,2.46E-03},
      {1.27E+02,1.04E+00,5.09E-03,7.86E-03,2.49E-02,2.43E-03},
      {1.28E+02,1.03E+00,4.93E-03,6.58E-03,2.46E-02,2.41E-03},
      {1.29E+02,1.03E+00,4.77E-03,5.33E-03,2.43E-02,2.38E-03},
      {1.30E+02,1.03E+00,4.62E-03,4.14E-03,2.39E-02,2.36E-03},
      {1.31E+02,1.03E+00,4.47E-03,2.98E-03,2.36E-02,2.34E-03},
      {1.32E+02,1.03E+00,4.33E-03,1.87E-03,2.33E-02,2.31E-03},
      {1.33E+02,1.03E+00,4.20E-03,7.93E-04,2.30E-02,2.29E-03},
      {1.34E+02,1.02E+00,4.07E-03,-2.42E-04,2.27E-02,2.27E-03},
      {1.35E+02,1.02E+00,3.94E-03,-1.24E-03,2.24E-02,2.25E-03},
      {1.36E+02,1.02E+00,3.82E-03,-2.20E-03,2.21E-02,2.23E-03},
      {1.37E+02,1.02E+00,3.70E-03,-3.12E-03,2.18E-02,2.20E-03},
      {1.38E+02,1.02E+00,3.59E-03,-4.01E-03,2.15E-02,2.18E-03},
      {1.39E+02,1.02E+00,3.48E-03,-4.86E-03,2.13E-02,2.16E-03},
      {1.40E+02,1.02E+00,3.38E-03,-5.68E-03,2.10E-02,2.14E-03},
      {1.41E+02,1.02E+00,3.28E-03,-6.47E-03,2.07E-02,2.12E-03},
      {1.42E+02,1.02E+00,3.18E-03,-7.22E-03,2.05E-02,2.11E-03},
      {1.43E+02,1.01E+00,3.09E-03,-7.95E-03,2.02E-02,2.09E-03},
      {1.44E+02,1.01E+00,3.00E-03,-8.65E-03,2.00E-02,2.07E-03},
      {1.45E+02,1.01E+00,2.92E-03,-9.32E-03,1.97E-02,2.05E-03},
      {1.46E+02,1.01E+00,2.83E-03,-9.96E-03,1.95E-02,2.03E-03},
      {1.47E+02,1.01E+00,2.76E-03,-1.06E-02,1.92E-02,2.01E-03},
      {1.48E+02,1.01E+00,2.68E-03,-1.12E-02,1.90E-02,2.00E-03},
      {1.49E+02,1.01E+00,2.61E-03,-1.17E-02,1.88E-02,1.98E-03},
      {1.50E+02,1.01E+00,2.54E-03,-1.23E-02,1.85E-02,1.96E-03},
      {1.51E+02,1.01E+00,2.48E-03,-1.28E-02,1.83E-02,1.95E-03},
      {1.52E+02,1.01E+00,2.42E-03,-1.33E-02,1.81E-02,1.93E-03},
      {1.53E+02,1.01E+00,2.36E-03,-1.37E-02,1.79E-02,1.91E-03},
      {1.54E+02,1.01E+00,2.31E-03,-1.42E-02,1.77E-02,1.90E-03},
      {1.55E+02,1.00E+00,2.26E-03,-1.46E-02,1.75E-02,1.88E-03},
      {1.56E+02,1.00E+00,2.21E-03,-1.50E-02,1.73E-02,1.87E-03},
      {1.57E+02,1.00E+00,2.16E-03,-1.54E-02,1.71E-02,1.85E-03},
      {1.58E+02,1.00E+00,2.12E-03,-1.58E-02,1.69E-02,1.84E-03},
      {1.59E+02,1.00E+00,2.07E-03,-1.61E-02,1.67E-02,1.82E-03},
      {1.60E+02,1.00E+00,2.03E-03,-1.64E-02,1.65E-02,1.81E-03},
      {1.61E+02,1.00E+00,2.00E-03,-1.68E-02,1.63E-02,1.79E-03},
      {1.62E+02,1.00E+00,1.96E-03,-1.71E-02,1.61E-02,1.78E-03},
      {1.63E+02,1.00E+00,1.93E-03,-1.73E-02,1.59E-02,1.77E-03},
      {1.64E+02,1.00E+00,1.90E-03,-1.76E-02,1.58E-02,1.75E-03},
      {1.65E+02,9.99E-01,1.87E-03,-1.78E-02,1.56E-02,1.74E-03},
      {1.66E+02,9.99E-01,1.84E-03,-1.81E-02,1.54E-02,1.73E-03},
      {1.67E+02,9.99E-01,1.81E-03,-1.83E-02,1.52E-02,1.71E-03},
      {1.68E+02,9.98E-01,1.79E-03,-1.85E-02,1.51E-02,1.70E-03},
      {1.69E+02,9.98E-01,1.76E-03,-1.87E-02,1.49E-02,1.69E-03},
      {1.70E+02,9.98E-01,1.74E-03,-1.89E-02,1.48E-02,1.68E-03},
      {1.71E+02,9.97E-01,1.71E-03,-1.91E-02,1.46E-02,1.66E-03},
      {1.72E+02,9.97E-01,1.69E-03,-1.92E-02,1.44E-02,1.65E-03},
      {1.73E+02,9.97E-01,1.67E-03,-1.94E-02,1.43E-02,1.64E-03},
      {1.74E+02,9.96E-01,1.65E-03,-1.95E-02,1.41E-02,1.63E-03},
      {1.75E+02,9.96E-01,1.63E-03,-1.96E-02,1.40E-02,1.62E-03},
      {1.76E+02,9.96E-01,1.62E-03,-1.98E-02,1.38E-02,1.61E-03},
      {1.77E+02,9.95E-01,1.60E-03,-1.99E-02,1.37E-02,1.59E-03},
      {1.78E+02,9.95E-01,1.58E-03,-2.00E-02,1.36E-02,1.58E-03},
      {1.79E+02,9.95E-01,1.56E-03,-2.00E-02,1.34E-02,1.57E-03},
      {1.80E+02,9.95E-01,1.55E-03,-2.01E-02,1.33E-02,1.56E-03},
      {1.81E+02,9.95E-01,1.53E-03,-2.02E-02,1.31E-02,1.55E-03},
      {1.82E+02,9.94E-01,1.52E-03,-2.03E-02,1.30E-02,1.54E-03},
      {1.83E+02,9.94E-01,1.50E-03,-2.03E-02,1.29E-02,1.53E-03},
      {1.84E+02,9.94E-01,1.49E-03,-2.04E-02,1.27E-02,1.52E-03},
      {1.85E+02,9.94E-01,1.47E-03,-2.04E-02,1.26E-02,1.51E-03},
      {1.86E+02,9.94E-01,1.46E-03,-2.04E-02,1.25E-02,1.50E-03},
      {1.87E+02,9.93E-01,1.44E-03,-2.05E-02,1.24E-02,1.49E-03},
      {1.88E+02,9.93E-01,1.43E-03,-2.05E-02,1.23E-02,1.48E-03},
      {1.89E+02,9.93E-01,1.42E-03,-2.05E-02,1.21E-02,1.47E-03},
      {1.90E+02,9.93E-01,1.40E-03,-2.05E-02,1.20E-02,1.46E-03},
      {1.91E+02,9.93E-01,1.39E-03,-2.05E-02,1.19E-02,1.45E-03},
      {1.92E+02,9.93E-01,1.38E-03,-2.05E-02,1.18E-02,1.44E-03},
      {1.93E+02,9.93E-01,1.36E-03,-2.05E-02,1.17E-02,1.43E-03},
      {1.94E+02,9.92E-01,1.35E-03,-2.05E-02,1.16E-02,1.43E-03},
      {1.95E+02,9.92E-01,1.34E-03,-2.05E-02,1.14E-02,1.42E-03},
      {1.96E+02,9.92E-01,1.33E-03,-2.05E-02,1.13E-02,1.41E-03},
      {1.97E+02,9.92E-01,1.31E-03,-2.04E-02,1.12E-02,1.40E-03},
      {1.98E+02,9.92E-01,1.30E-03,-2.04E-02,1.11E-02,1.39E-03},
      {1.99E+02,9.92E-01,1.29E-03,-2.04E-02,1.10E-02,1.38E-03},
      {2.00E+02,9.92E-01,1.28E-03,-2.03E-02,1.09E-02,1.37E-03},
      {2.01E+02,9.92E-01,1.26E-03,-2.03E-02,1.08E-02,1.37E-03},
      {2.02E+02,9.92E-01,1.25E-03,-2.03E-02,1.07E-02,1.36E-03},
      {2.03E+02,9.92E-01,1.24E-03,-2.02E-02,1.06E-02,1.35E-03},
      {2.04E+02,9.92E-01,1.23E-03,-2.02E-02,1.05E-02,1.34E-03},
      {2.05E+02,9.92E-01,1.22E-03,-2.01E-02,1.04E-02,1.34E-03},
      {2.06E+02,9.92E-01,1.21E-03,-2.01E-02,1.03E-02,1.33E-03},
      {2.07E+02,9.92E-01,1.19E-03,-2.00E-02,1.03E-02,1.32E-03},
      {2.08E+02,9.92E-01,1.18E-03,-1.99E-02,1.02E-02,1.31E-03},
      {2.09E+02,9.92E-01,1.17E-03,-1.99E-02,1.01E-02,1.31E-03},
      {2.10E+02,9.91E-01,1.16E-03,-1.98E-02,9.98E-03,1.30E-03},
      {2.11E+02,9.91E-01,1.15E-03,-1.97E-02,9.90E-03,1.29E-03},
      {2.12E+02,9.91E-01,1.14E-03,-1.97E-02,9.81E-03,1.28E-03},
      {2.13E+02,9.91E-01,1.13E-03,-1.96E-02,9.72E-03,1.28E-03},
      {2.14E+02,9.91E-01,1.12E-03,-1.95E-02,9.64E-03,1.27E-03},
      {2.15E+02,9.91E-01,1.10E-03,-1.95E-02,9.56E-03,1.26E-03},
      {2.16E+02,9.91E-01,1.09E-03,-1.94E-02,9.48E-03,1.26E-03},
      {2.17E+02,9.91E-01,1.08E-03,-1.93E-02,9.40E-03,1.25E-03},
      {2.18E+02,9.91E-01,1.07E-03,-1.92E-02,9.32E-03,1.24E-03},
      {2.19E+02,9.91E-01,1.06E-03,-1.91E-02,9.24E-03,1.24E-03},
      {2.20E+02,9.91E-01,1.05E-03,-1.91E-02,9.16E-03,1.23E-03},
      {2.21E+02,9.91E-01,1.04E-03,-1.90E-02,9.09E-03,1.22E-03},
      {2.22E+02,9.91E-01,1.03E-03,-1.89E-02,9.01E-03,1.22E-03},
      {2.23E+02,9.91E-01,1.02E-03,-1.88E-02,8.94E-03,1.21E-03},
      {2.24E+02,9.91E-01,1.01E-03,-1.87E-02,8.86E-03,1.21E-03},
      {2.25E+02,9.91E-01,1.00E-03,-1.86E-02,8.79E-03,1.20E-03},
      {2.26E+02,9.91E-01,9.90E-04,-1.86E-02,8.72E-03,1.19E-03},
      {2.27E+02,9.91E-01,9.80E-04,-1.85E-02,8.65E-03,1.19E-03},
      {2.28E+02,9.91E-01,9.70E-04,-1.84E-02,8.58E-03,1.18E-03},
      {2.29E+02,9.91E-01,9.60E-04,-1.83E-02,8.51E-03,1.18E-03},
      {2.30E+02,9.91E-01,9.50E-04,-1.82E-02,8.44E-03,1.17E-03},
      {2.31E+02,9.91E-01,9.41E-04,-1.81E-02,8.38E-03,1.16E-03},
      {2.32E+02,9.91E-01,9.31E-04,-1.80E-02,8.31E-03,1.16E-03},
      {2.33E+02,9.91E-01,9.22E-04,-1.79E-02,8.25E-03,1.15E-03},
      {2.34E+02,9.92E-01,9.12E-04,-1.78E-02,8.18E-03,1.15E-03},
      {2.35E+02,9.92E-01,9.03E-04,-1.77E-02,8.12E-03,1.14E-03},
      {2.36E+02,9.92E-01,8.94E-04,-1.76E-02,8.06E-03,1.14E-03},
      {2.37E+02,9.92E-01,8.84E-04,-1.75E-02,7.99E-03,1.13E-03},
      {2.38E+02,9.92E-01,8.75E-04,-1.75E-02,7.93E-03,1.13E-03},
      {2.39E+02,9.92E-01,8.66E-04,-1.74E-02,7.87E-03,1.12E-03},
      {2.40E+02,9.92E-01,8.57E-04,-1.73E-02,7.81E-03,1.12E-03},
      {2.41E+02,9.92E-01,8.48E-04,-1.72E-02,7.75E-03,1.11E-03},
      {2.42E+02,9.92E-01,8.40E-04,-1.71E-02,7.70E-03,1.11E-03},
      {2.43E+02,9.92E-01,8.31E-04,-1.70E-02,7.64E-03,1.10E-03},
      {2.44E+02,9.92E-01,8.22E-04,-1.69E-02,7.58E-03,1.10E-03},
      {2.45E+02,9.92E-01,8.14E-04,-1.68E-02,7.53E-03,1.09E-03},
      {2.46E+02,9.92E-01,8.05E-04,-1.67E-02,7.47E-03,1.09E-03},
      {2.47E+02,9.92E-01,7.97E-04,-1.66E-02,7.42E-03,1.08E-03},
      {2.48E+02,9.92E-01,7.89E-04,-1.65E-02,7.36E-03,1.08E-03},
      {2.49E+02,9.92E-01,7.80E-04,-1.64E-02,7.31E-03,1.07E-03},
      {2.50E+02,9.92E-01,7.72E-04,-1.63E-02,7.25E-03,1.07E-03},
      {2.51E+02,9.92E-01,7.64E-04,-1.62E-02,7.20E-03,1.06E-03},
      {2.52E+02,9.92E-01,7.56E-04,-1.61E-02,7.15E-03,1.06E-03},
      {2.53E+02,9.92E-01,7.48E-04,-1.60E-02,7.10E-03,1.05E-03},
      {2.54E+02,9.92E-01,7.40E-04,-1.59E-02,7.05E-03,1.05E-03},
      {2.55E+02,9.92E-01,7.32E-04,-1.58E-02,7.00E-03,1.05E-03},
      {2.56E+02,9.92E-01,7.25E-04,-1.58E-02,6.95E-03,1.04E-03},
      {2.57E+02,9.92E-01,7.17E-04,-1.57E-02,6.90E-03,1.04E-03},
      {2.58E+02,9.92E-01,7.10E-04,-1.56E-02,6.85E-03,1.03E-03},
      {2.59E+02,9.92E-01,7.02E-04,-1.55E-02,6.81E-03,1.03E-03},
      {2.60E+02,9.92E-01,6.95E-04,-1.54E-02,6.76E-03,1.02E-03},
      {2.61E+02,9.92E-01,6.87E-04,-1.53E-02,6.71E-03,1.02E-03},
      {2.62E+02,9.92E-01,6.80E-04,-1.52E-02,6.67E-03,1.02E-03},
      {2.63E+02,9.93E-01,6.73E-04,-1.51E-02,6.62E-03,1.01E-03},
      {2.64E+02,9.93E-01,6.66E-04,-1.50E-02,6.58E-03,1.01E-03},
      {2.65E+02,9.93E-01,6.59E-04,-1.49E-02,6.53E-03,1.00E-03},
      {2.66E+02,9.93E-01,6.52E-04,-1.48E-02,6.49E-03,9.99E-04},
      {2.67E+02,9.93E-01,6.45E-04,-1.47E-02,6.44E-03,9.95E-04},
      {2.68E+02,9.93E-01,6.38E-04,-1.47E-02,6.40E-03,9.92E-04},
      {2.69E+02,9.93E-01,6.31E-04,-1.46E-02,6.36E-03,9.88E-04},
      {2.70E+02,9.93E-01,6.24E-04,-1.45E-02,6.32E-03,9.84E-04},
      {2.71E+02,9.93E-01,6.18E-04,-1.44E-02,6.27E-03,9.80E-04},
      {2.72E+02,9.93E-01,6.11E-04,-1.43E-02,6.23E-03,9.76E-04},
      {2.73E+02,9.93E-01,6.05E-04,-1.42E-02,6.19E-03,9.73E-04},
      {2.74E+02,9.93E-01,5.98E-04,-1.41E-02,6.15E-03,9.69E-04},
      {2.75E+02,9.93E-01,5.92E-04,-1.40E-02,6.11E-03,9.65E-04},
      {2.76E+02,9.93E-01,5.85E-04,-1.40E-02,6.07E-03,9.62E-04},
      {2.77E+02,9.93E-01,5.79E-04,-1.39E-02,6.03E-03,9.58E-04},
      {2.78E+02,9.93E-01,5.73E-04,-1.38E-02,5.99E-03,9.54E-04},
      {2.79E+02,9.93E-01,5.67E-04,-1.37E-02,5.96E-03,9.51E-04},
      {2.80E+02,9.93E-01,5.61E-04,-1.36E-02,5.92E-03,9.47E-04},
      {2.81E+02,9.93E-01,5.55E-04,-1.35E-02,5.88E-03,9.44E-04},
      {2.82E+02,9.93E-01,5.49E-04,-1.34E-02,5.84E-03,9.40E-04},
      {2.83E+02,9.93E-01,5.43E-04,-1.34E-02,5.81E-03,9.37E-04},
      {2.84E+02,9.93E-01,5.37E-04,-1.33E-02,5.77E-03,9.34E-04},
      {2.85E+02,9.93E-01,5.31E-04,-1.32E-02,5.73E-03,9.30E-04},
      {2.86E+02,9.94E-01,5.26E-04,-1.31E-02,5.70E-03,9.27E-04},
      {2.87E+02,9.94E-01,5.20E-04,-1.30E-02,5.66E-03,9.24E-04},
      {2.88E+02,9.94E-01,5.15E-04,-1.30E-02,5.63E-03,9.20E-04},
      {2.89E+02,9.94E-01,5.09E-04,-1.29E-02,5.59E-03,9.17E-04},
      {2.90E+02,9.94E-01,5.04E-04,-1.28E-02,5.56E-03,9.14E-04},
      {2.91E+02,9.94E-01,4.98E-04,-1.27E-02,5.53E-03,9.11E-04},
      {2.92E+02,9.94E-01,4.93E-04,-1.26E-02,5.49E-03,9.07E-04},
      {2.93E+02,9.94E-01,4.88E-04,-1.26E-02,5.46E-03,9.04E-04},
      {2.94E+02,9.94E-01,4.82E-04,-1.25E-02,5.43E-03,9.01E-04},
      {2.95E+02,9.94E-01,4.77E-04,-1.24E-02,5.39E-03,8.98E-04},
      {2.96E+02,9.94E-01,4.72E-04,-1.23E-02,5.36E-03,8.95E-04},
      {2.97E+02,9.94E-01,4.67E-04,-1.22E-02,5.33E-03,8.92E-04},
      {2.98E+02,9.94E-01,4.62E-04,-1.22E-02,5.30E-03,8.89E-04},
      {2.99E+02,9.94E-01,4.57E-04,-1.21E-02,5.27E-03,8.86E-04},
      {3.00E+02,9.94E-01,4.52E-04,-1.20E-02,5.24E-03,8.83E-04}
    };
  //k* [MeV/c](0) C_pXi- Added (1) Error (2) C_pXi- (3) C_nXi0->pXi0 (4) C_LL->pXi0 (5)
  TGraphErrors* GrYuki_0_97_Added = new TGraphErrors(Yuki_0_97.size());
  TGraphErrors* GrYuki_0_97_pXim  = new TGraphErrors(Yuki_0_97.size());
  TGraphErrors* GrYuki_0_97_nXi0  = new TGraphErrors(Yuki_0_97.size());
  TGraphErrors* GrYuki_0_97_LaLa  = new TGraphErrors(Yuki_0_97.size());
  int counter = 0;
  for (auto it : Yuki_0_97) {
    GrYuki_0_97_Added->SetPoint(counter,it[0],it[1]);
    GrYuki_0_97_Added->SetPointError(counter,0,it[2]);
    
    GrYuki_0_97_pXim->SetPoint(counter, it[0], it[3]+1); 
    GrYuki_0_97_nXi0->SetPoint(counter, it[0], it[4]+1); 
    GrYuki_0_97_LaLa->SetPoint(counter, it[0], it[5]+1); 
    
    counter++; 
  }
  TGraphErrors* GrYuki_1_02_Added = new TGraphErrors(Yuki_1_02.size());
  TGraphErrors* GrYuki_1_02_pXim  = new TGraphErrors(Yuki_1_02.size());
  TGraphErrors* GrYuki_1_02_nXi0  = new TGraphErrors(Yuki_1_02.size());
  TGraphErrors* GrYuki_1_02_LaLa  = new TGraphErrors(Yuki_1_02.size());
  counter = 0; 
  for (auto it : Yuki_1_02) {
    GrYuki_1_02_Added->SetPoint(counter,it[0],it[1]);
    GrYuki_1_02_Added->SetPointError(counter,0,it[2]);
    
    GrYuki_1_02_pXim->SetPoint(counter, it[0], it[3]+1); 
    GrYuki_1_02_nXi0->SetPoint(counter, it[0], it[4]+1); 
    GrYuki_1_02_LaLa->SetPoint(counter, it[0], it[5]+1); 
    
    counter++; 
  }
  TGraphErrors* GrYuki_1_07_Added = new TGraphErrors(Yuki_1_07.size());
  TGraphErrors* GrYuki_1_07_pXim  = new TGraphErrors(Yuki_1_07.size());
  TGraphErrors* GrYuki_1_07_nXi0  = new TGraphErrors(Yuki_1_07.size());
  TGraphErrors* GrYuki_1_07_LaLa  = new TGraphErrors(Yuki_1_07.size());
  counter = 0; 
  for (auto it : Yuki_1_07) {
    GrYuki_1_07_Added->SetPoint(counter,it[0],it[1]);
    GrYuki_1_07_Added->SetPointError(counter,0,it[2]);
    
    GrYuki_1_07_pXim->SetPoint(counter, it[0], it[3]+1); 
    GrYuki_1_07_nXi0->SetPoint(counter, it[0], it[4]+1); 
    GrYuki_1_07_LaLa->SetPoint(counter, it[0], it[5]+1); 
    
    counter++; 
  }

  out->cd();
  GrYuki_0_97_Added->Write("GrYuki_0_97_Added");
  GrYuki_0_97_pXim->Write("GrYuki_0_97_pXim");
  GrYuki_0_97_nXi0->Write("GrYuki_0_97_nXi0");
  GrYuki_0_97_LaLa->Write("GrYuki_0_97_LaLa");

  GrYuki_1_02_Added->Write("GrYuki_1_02_Added");
  GrYuki_1_02_pXim->Write("GrYuki_1_02_pXim");
  GrYuki_1_02_nXi0->Write("GrYuki_1_02_nXi0");
  GrYuki_1_02_LaLa->Write("GrYuki_1_02_LaLa");
  
  GrYuki_1_07_Added->Write("GrYuki_1_07_Added");
  GrYuki_1_07_pXim->Write("GrYuki_1_07_pXim");
  GrYuki_1_07_nXi0->Write("GrYuki_1_07_nXi0");
  GrYuki_1_07_LaLa->Write("GrYuki_1_07_LaLa");
 
  return;
} 


void PlayWithCats::GenerateYukiCurves_200515(TFile* out) {
 std:vector<std::vector<float>> Yuki_102 = {
    { 0.10000E+01, 0.76926E+02},
    { 0.20000E+01, 0.38169E+02},
    { 0.30000E+01, 0.25253E+02},
    { 0.40000E+01, 0.18742E+02},
    { 0.50000E+01, 0.14800E+02},
    { 0.60000E+01, 0.12162E+02},
    { 0.70000E+01, 0.10283E+02},
    { 0.80000E+01, 0.88858E+01},
    { 0.90000E+01, 0.78116E+01},
    { 0.10000E+02, 0.69647E+01},
    { 0.11000E+02, 0.62825E+01},
    { 0.12000E+02, 0.57231E+01},
    { 0.13000E+02, 0.52571E+01},
    { 0.14000E+02, 0.48640E+01},
    { 0.15000E+02, 0.45283E+01},
    { 0.16000E+02, 0.42387E+01},
    { 0.17000E+02, 0.39865E+01},
    { 0.18000E+02, 0.37652E+01},
    { 0.19000E+02, 0.35695E+01},
    { 0.20000E+02, 0.33954E+01},
    { 0.21000E+02, 0.32395E+01},
    { 0.22000E+02, 0.30992E+01},
    { 0.23000E+02, 0.29723E+01},
    { 0.24000E+02, 0.28570E+01},
    { 0.25000E+02, 0.27519E+01},
    { 0.26000E+02, 0.26556E+01},
    { 0.27000E+02, 0.25672E+01},
    { 0.28000E+02, 0.24857E+01},
    { 0.29000E+02, 0.24104E+01},
    { 0.30000E+02, 0.23407E+01},
    { 0.31000E+02, 0.22758E+01},
    { 0.32000E+02, 0.22155E+01},
    { 0.33000E+02, 0.21591E+01},
    { 0.34000E+02, 0.21065E+01},
    { 0.35000E+02, 0.20571E+01},
    { 0.36000E+02, 0.20108E+01},
    { 0.37000E+02, 0.19672E+01},
    { 0.38000E+02, 0.19262E+01},
    { 0.39000E+02, 0.18875E+01},
    { 0.40000E+02, 0.18510E+01},
    { 0.41000E+02, 0.18164E+01},
    { 0.42000E+02, 0.17837E+01},
    { 0.43000E+02, 0.17526E+01},
    { 0.44000E+02, 0.17232E+01},
    { 0.45000E+02, 0.16952E+01},
    { 0.46000E+02, 0.16685E+01},
    { 0.47000E+02, 0.16432E+01},
    { 0.48000E+02, 0.16190E+01},
    { 0.49000E+02, 0.15959E+01},
    { 0.50000E+02, 0.15739E+01},
    { 0.51000E+02, 0.15529E+01},
    { 0.52000E+02, 0.15327E+01},
    { 0.53000E+02, 0.15135E+01},
    { 0.54000E+02, 0.14950E+01},
    { 0.55000E+02, 0.14773E+01},
    { 0.56000E+02, 0.14604E+01},
    { 0.57000E+02, 0.14441E+01},
    { 0.58000E+02, 0.14285E+01},
    { 0.59000E+02, 0.14134E+01},
    { 0.60000E+02, 0.13990E+01},
    { 0.61000E+02, 0.13851E+01},
    { 0.62000E+02, 0.13717E+01},
    { 0.63000E+02, 0.13589E+01},
    { 0.64000E+02, 0.13465E+01},
    { 0.65000E+02, 0.13346E+01},
    { 0.66000E+02, 0.13230E+01},
    { 0.67000E+02, 0.13119E+01},
    { 0.68000E+02, 0.13012E+01},
    { 0.69000E+02, 0.12909E+01},
    { 0.70000E+02, 0.12809E+01},
    { 0.71000E+02, 0.12713E+01},
    { 0.72000E+02, 0.12620E+01},
    { 0.73000E+02, 0.12530E+01},
    { 0.74000E+02, 0.12443E+01},
    { 0.75000E+02, 0.12359E+01},
    { 0.76000E+02, 0.12277E+01},
    { 0.77000E+02, 0.12199E+01},
    { 0.78000E+02, 0.12123E+01},
    { 0.79000E+02, 0.12049E+01},
    { 0.80000E+02, 0.11978E+01},
    { 0.81000E+02, 0.11909E+01},
    { 0.82000E+02, 0.11842E+01},
    { 0.83000E+02, 0.11778E+01},
    { 0.84000E+02, 0.11715E+01},
    { 0.85000E+02, 0.11654E+01},
    { 0.86000E+02, 0.11595E+01},
    { 0.87000E+02, 0.11539E+01},
    { 0.88000E+02, 0.11483E+01},
    { 0.89000E+02, 0.11430E+01},
    { 0.90000E+02, 0.11378E+01},
    { 0.91000E+02, 0.11328E+01},
    { 0.92000E+02, 0.11279E+01},
    { 0.93000E+02, 0.11232E+01},
    { 0.94000E+02, 0.11186E+01},
    { 0.95000E+02, 0.11142E+01},
    { 0.96000E+02, 0.11099E+01},
    { 0.97000E+02, 0.11057E+01},
    { 0.98000E+02, 0.11017E+01},
    { 0.99000E+02, 0.10978E+01},
    { 0.10000E+03, 0.10940E+01},
    { 0.10100E+03, 0.10903E+01},
    { 0.10200E+03, 0.10868E+01},
    { 0.10300E+03, 0.10833E+01},
    { 0.10400E+03, 0.10799E+01},
    { 0.10500E+03, 0.10767E+01},
    { 0.10600E+03, 0.10735E+01},
    { 0.10700E+03, 0.10705E+01},
    { 0.10800E+03, 0.10675E+01},
    { 0.10900E+03, 0.10647E+01},
    { 0.11000E+03, 0.10619E+01},
    { 0.11100E+03, 0.10592E+01},
    { 0.11200E+03, 0.10566E+01},
    { 0.11300E+03, 0.10540E+01},
    { 0.11400E+03, 0.10516E+01},
    { 0.11500E+03, 0.10492E+01},
    { 0.11600E+03, 0.10469E+01},
    { 0.11700E+03, 0.10447E+01},
    { 0.11800E+03, 0.10425E+01},
    { 0.11900E+03, 0.10405E+01},
    { 0.12000E+03, 0.10384E+01},
    { 0.12100E+03, 0.10365E+01},
    { 0.12200E+03, 0.10346E+01},
    { 0.12300E+03, 0.10327E+01},
    { 0.12400E+03, 0.10310E+01},
    { 0.12500E+03, 0.10292E+01},
    { 0.12600E+03, 0.10276E+01},
    { 0.12700E+03, 0.10260E+01},
    { 0.12800E+03, 0.10244E+01},
    { 0.12900E+03, 0.10229E+01},
    { 0.13000E+03, 0.10215E+01},
    { 0.13100E+03, 0.10201E+01},
    { 0.13200E+03, 0.10187E+01},
    { 0.13300E+03, 0.10174E+01},
    { 0.13400E+03, 0.10161E+01},
    { 0.13500E+03, 0.10149E+01},
    { 0.13600E+03, 0.10137E+01},
    { 0.13700E+03, 0.10126E+01},
    { 0.13800E+03, 0.10115E+01},
    { 0.13900E+03, 0.10104E+01},
    { 0.14000E+03, 0.10094E+01},
    { 0.14100E+03, 0.10084E+01},
    { 0.14200E+03, 0.10074E+01},
    { 0.14300E+03, 0.10065E+01},
    { 0.14400E+03, 0.10056E+01},
    { 0.14500E+03, 0.10047E+01},
    { 0.14600E+03, 0.10039E+01},
    { 0.14700E+03, 0.10031E+01},
    { 0.14800E+03, 0.10024E+01},
    { 0.14900E+03, 0.10016E+01},
    { 0.15000E+03, 0.10009E+01},
    { 0.15100E+03, 0.10002E+01},
    { 0.15200E+03, 0.99956E+00},
    { 0.15300E+03, 0.99893E+00},
    { 0.15400E+03, 0.99832E+00},
    { 0.15500E+03, 0.99774E+00},
    { 0.15600E+03, 0.99717E+00},
    { 0.15700E+03, 0.99663E+00},
    { 0.15800E+03, 0.99612E+00},
    { 0.15900E+03, 0.99562E+00},
    { 0.16000E+03, 0.99514E+00},
    { 0.16100E+03, 0.99469E+00},
    { 0.16200E+03, 0.99425E+00},
    { 0.16300E+03, 0.99383E+00},
    { 0.16400E+03, 0.99343E+00},
    { 0.16500E+03, 0.99305E+00},
    { 0.16600E+03, 0.99268E+00},
    { 0.16700E+03, 0.99233E+00},
    { 0.16800E+03, 0.99199E+00},
    { 0.16900E+03, 0.99167E+00},
    { 0.17000E+03, 0.99137E+00},
    { 0.17100E+03, 0.99107E+00},
    { 0.17200E+03, 0.99080E+00},
    { 0.17300E+03, 0.99053E+00},
    { 0.17400E+03, 0.99028E+00},
    { 0.17500E+03, 0.99004E+00},
    { 0.17600E+03, 0.98982E+00},
    { 0.17700E+03, 0.98960E+00},
    { 0.17800E+03, 0.98940E+00},
    { 0.17900E+03, 0.98920E+00},
    { 0.18000E+03, 0.98902E+00},
    { 0.18100E+03, 0.98885E+00},
    { 0.18200E+03, 0.98868E+00},
    { 0.18300E+03, 0.98853E+00},
    { 0.18400E+03, 0.98838E+00},
    { 0.18500E+03, 0.98825E+00},
    { 0.18600E+03, 0.98812E+00},
    { 0.18700E+03, 0.98800E+00},
    { 0.18800E+03, 0.98789E+00},
    { 0.18900E+03, 0.98779E+00},
    { 0.19000E+03, 0.98769E+00},
    { 0.19100E+03, 0.98760E+00},
    { 0.19200E+03, 0.98752E+00},
    { 0.19300E+03, 0.98744E+00},
    { 0.19400E+03, 0.98737E+00},
    { 0.19500E+03, 0.98731E+00},
    { 0.19600E+03, 0.98725E+00},
    { 0.19700E+03, 0.98720E+00},
    { 0.19800E+03, 0.98715E+00},
    { 0.19900E+03, 0.98711E+00},
    { 0.20000E+03, 0.98707E+00},
    { 0.20100E+03, 0.98704E+00},
    { 0.20200E+03, 0.98701E+00},
    { 0.20300E+03, 0.98699E+00},
    { 0.20400E+03, 0.98697E+00},
    { 0.20500E+03, 0.98696E+00},
    { 0.20600E+03, 0.98695E+00},
    { 0.20700E+03, 0.98694E+00},
    { 0.20800E+03, 0.98694E+00},
    { 0.20900E+03, 0.98694E+00},
    { 0.21000E+03, 0.98694E+00},
    { 0.21100E+03, 0.98695E+00},
    { 0.21200E+03, 0.98696E+00},
    { 0.21300E+03, 0.98697E+00},
    { 0.21400E+03, 0.98699E+00},
    { 0.21500E+03, 0.98701E+00},
    { 0.21600E+03, 0.98703E+00},
    { 0.21700E+03, 0.98705E+00},
    { 0.21800E+03, 0.98708E+00},
    { 0.21900E+03, 0.98711E+00},
    { 0.22000E+03, 0.98714E+00},
    { 0.22100E+03, 0.98717E+00},
    { 0.22200E+03, 0.98720E+00},
    { 0.22300E+03, 0.98724E+00},
    { 0.22400E+03, 0.98728E+00},
    { 0.22500E+03, 0.98732E+00},
    { 0.22600E+03, 0.98736E+00},
    { 0.22700E+03, 0.98740E+00},
    { 0.22800E+03, 0.98744E+00},
    { 0.22900E+03, 0.98749E+00},
    { 0.23000E+03, 0.98754E+00},
    { 0.23100E+03, 0.98758E+00},
    { 0.23200E+03, 0.98763E+00},
    { 0.23300E+03, 0.98768E+00},
    { 0.23400E+03, 0.98774E+00},
    { 0.23500E+03, 0.98779E+00},
    { 0.23600E+03, 0.98784E+00},
    { 0.23700E+03, 0.98790E+00},
    { 0.23800E+03, 0.98795E+00},
    { 0.23900E+03, 0.98801E+00},
    { 0.24000E+03, 0.98807E+00},
    { 0.24100E+03, 0.98812E+00},
    { 0.24200E+03, 0.98818E+00},
    { 0.24300E+03, 0.98824E+00},
    { 0.24400E+03, 0.98830E+00},
    { 0.24500E+03, 0.98836E+00},
    { 0.24600E+03, 0.98842E+00},
    { 0.24700E+03, 0.98848E+00},
    { 0.24800E+03, 0.98854E+00},
    { 0.24900E+03, 0.98860E+00},
    { 0.25000E+03, 0.98866E+00},
    { 0.25100E+03, 0.98873E+00},
    { 0.25200E+03, 0.98879E+00},
    { 0.25300E+03, 0.98885E+00},
    { 0.25400E+03, 0.98892E+00},
    { 0.25500E+03, 0.98898E+00},
    { 0.25600E+03, 0.98904E+00},
    { 0.25700E+03, 0.98911E+00},
    { 0.25800E+03, 0.98917E+00},
    { 0.25900E+03, 0.98923E+00},
    { 0.26000E+03, 0.98930E+00},
    { 0.26100E+03, 0.98936E+00},
    { 0.26200E+03, 0.98943E+00},
    { 0.26300E+03, 0.98949E+00},
    { 0.26400E+03, 0.98956E+00},
    { 0.26500E+03, 0.98962E+00},
    { 0.26600E+03, 0.98968E+00},
    { 0.26700E+03, 0.98975E+00},
    { 0.26800E+03, 0.98981E+00},
    { 0.26900E+03, 0.98988E+00},
    { 0.27000E+03, 0.98994E+00},
    { 0.27100E+03, 0.99000E+00},
    { 0.27200E+03, 0.99007E+00},
    { 0.27300E+03, 0.99013E+00},
    { 0.27400E+03, 0.99020E+00},
    { 0.27500E+03, 0.99026E+00},
    { 0.27600E+03, 0.99032E+00},
    { 0.27700E+03, 0.99039E+00},
    { 0.27800E+03, 0.99045E+00},
    { 0.27900E+03, 0.99051E+00},
    { 0.28000E+03, 0.99058E+00},
    { 0.28100E+03, 0.99064E+00},
    { 0.28200E+03, 0.99070E+00},
    { 0.28300E+03, 0.99076E+00},
    { 0.28400E+03, 0.99083E+00},
    { 0.28500E+03, 0.99089E+00},
    { 0.28600E+03, 0.99095E+00},
    { 0.28700E+03, 0.99101E+00},
    { 0.28800E+03, 0.99107E+00},
    { 0.28900E+03, 0.99113E+00},
    { 0.29000E+03, 0.99119E+00},
    { 0.29100E+03, 0.99125E+00},
    { 0.29200E+03, 0.99131E+00},
    { 0.29300E+03, 0.99137E+00},
    { 0.29400E+03, 0.99143E+00},
    { 0.29500E+03, 0.99149E+00},
    { 0.29600E+03, 0.99155E+00},
    { 0.29700E+03, 0.99161E+00},
    { 0.29800E+03, 0.99167E+00},
    { 0.29900E+03, 0.99173E+00},
    { 0.30000E+03, 0.99179E+00}
  };
  TGraphErrors* GrYuki_1_02_Added = new TGraphErrors(Yuki_102.size());
  int counter = 0; 
  for (auto it : Yuki_102) {
    GrYuki_1_02_Added->SetPoint(counter,it[0],it[1]);
    //GrYuki_1_02_Added->SetPointError(counter,0,it[2]);    
    counter++; 
  }
  out->cd();
  GrYuki_1_02_Added->Write("GrYuki_1_02_Added"); 
  return; 
} 
