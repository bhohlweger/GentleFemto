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
//  std::vector<float> cutOffValues = { 0. };
  std::vector<float> cutOffValues = { 0., 0.1, 0.3, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0,
      1.1, 1.3, 1.5, 2.0};
//  std::vector<float> cutOffValues = { 0.3, 0.7, 1.0, 1.02, 1.06, 1.1, 1.3, 5.0 };
  const int nCutOffs = cutOffValues.size();
  TH1F* CF_I0S0[nCutOffs];
//  TH2F* CF_I0S0WaveFunction[nCutOffs];
  TH1F* CF_I0S1[nCutOffs];
//  TH2F* CF_I0S1WaveFunction[nCutOffs];
  TH1F* CF_I1S0[nCutOffs];
//  TH2F* CF_I1S0WaveFunction[nCutOffs];
  TH1F* CF_I1S1[nCutOffs];
//  TH2F* CF_I1S1WaveFunction[nCutOffs];
  TLine one = TLine(0., 1., 400., 1.);
  double QCDTime = 11;
//4th argument is the t parameter and can be:
// 9, 10, 11, 12

  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double Mass_Xim = TDatabasePDG::Instance()->GetParticle(3312)->Mass()
      * 1000;
  int RadBins = 1000;
  float rMin = 0.;
  float rMax = 10.;
  float dRad = (rMax - rMin)/(float)RadBins;
  int momBins = 80;
  float kMin = 0.;
  float kMax = 400;
  int iCutOff = 0;
  TFile* outfile = TFile::Open("output.root", "recreate");
  auto c1 = new TCanvas("c1", "c1");
  auto c2 = new TCanvas("c2", "c2");
  auto c3 = new TCanvas("c3", "c3");
  auto c4 = new TCanvas("c4", "c4");
  TLegend leg(0.7, 0.5, 0.9, 0.9);
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
    int COLORMEBLIND = 2 + 10 * iCutOff;
    if (COLORMEBLIND > 100) {
      COLORMEBLIND -= 99;
    }
//    int COLORMEBLIND = 1 + iCutOff;
    double pXimPotParsI0S0[11] = { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 0, 0,
        0, it };
    double pXimPotParsI0S1[11] = { 0, 0, pXim_HALQCD1, QCDTime, 0, -1, 1, 1, 0,
        1, it };
    double pXimPotParsI1S0[11] =
        { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 0, 0, 0, it };
    double pXimPotParsI1S1[11] =
        { 0, 0, pXim_HALQCD1, QCDTime, 1, 1, 1, 1, 0, 1, it };
//    static double pXimPotParsI0S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0,
//        0 };
//    static double pXimPotParsI0S1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
//        0 };
//    static double pXimPotParsI1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
//        1 };
//    static double pXimPotParsI1S1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
//        2 };

    TString I0S0HistName = Form("hI0S0_%.2f", it);
    TString I0S0HistTitle = Form("I0S0");
    CF_I0S0[iCutOff] = new TH1F(I0S0HistName.Data(), I0S0HistTitle.Data(),
                                momBins, kMin, kMax);
//    I0S0HistName+="WF";
//    CF_I0S0WaveFunction[iCutOff] = new TH2F(I0S0HistName.Data(), I0S0HistTitle.Data(),
//                                momBins, kMin, kMax,RadBins,rMin,rMax);

    TString I0S1HistName = Form("hI0S1_%.2f", it);
    TString I0S1HistTitle = Form("I0S1");
    CF_I0S1[iCutOff] = new TH1F(I0S1HistName.Data(), I0S1HistTitle.Data(),
                                momBins, kMin, kMax);
//    I0S1HistName+="WF";
//    CF_I0S1WaveFunction[iCutOff] = new TH2F(I0S1HistName.Data(), I0S1HistTitle.Data(),
//                                momBins, kMin, kMax,RadBins,rMin,rMax);

    TString I1S0HistName = Form("hI1S0_%.2f", it);
    TString I1S0HistTitle = Form("I1S0");
    CF_I1S0[iCutOff] = new TH1F(I1S0HistName.Data(), I1S0HistTitle.Data(),
                                momBins, kMin, kMax);
//    I1S0HistName+="WF";
//    CF_I1S0WaveFunction[iCutOff] = new TH2F(I1S0HistName.Data(), I1S0HistTitle.Data(),
//                                momBins, kMin, kMax,RadBins,rMin,rMax);

    TString I1S1HistName = Form("hI1S1_%.2f", it);
    TString I1S1HistTitle = Form("I1S1");
    CF_I1S1[iCutOff] = new TH1F(I1S1HistName.Data(), I1S1HistTitle.Data(),
                                momBins, kMin, kMax);
//    I1S1HistName+="WF";
//    CF_I1S1WaveFunction[iCutOff] = new TH2F(I1S1HistName.Data(), I1S1HistTitle.Data(),
//                                momBins, kMin, kMax,RadBins,rMin,rMax);

    double GaussSourceSize = 1.4;

    CATS AB_pXim_I0S0;
    double Pars_pXiI0S0[4] = { 0, 0, 0, GaussSourceSize };
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
    AB_pXim_I0S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S0);

    AB_pXim_I0S0.KillTheCat();

    CATS AB_pXim_I0S1;
    double Pars_pXiI0S1[4] = { 0, 0, 0, GaussSourceSize };
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
    AB_pXim_I0S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI0S1);

    AB_pXim_I0S1.KillTheCat();

    CATS AB_pXim_I1S0;
    double Pars_pXiI1S0[4] = { 0, 0, 0, GaussSourceSize };
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
    AB_pXim_I1S0.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S0);

    AB_pXim_I1S0.KillTheCat();

    CATS AB_pXim_I1S1;
    double Pars_pXiI1S1[4] = { 0, 0, 0, GaussSourceSize };
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
    AB_pXim_I1S1.SetShortRangePotential(0, 0, fDlmPot, pXimPotParsI1S1);

    AB_pXim_I1S1.KillTheCat();

    for (int ikS = 0; ikS < momBins; ++ikS) {
//      std::cout << "=======" << ikS << "======= \n";
//      std::cout << "GetCorrFunction: " << AB_pXim_I0S0.GetCorrFun(ikS + 1)
//                << std::endl;
//      std::cout
//          << "EvalBinnedFun: "
//          << AB_pXim_I0S0.EvalCorrFun(CF_I0S0[iCutOff]->GetBinCenter(ikS + 1))
//          << std::endl;
      CF_I0S0[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I0S0.GetCorrFun(ikS));
      CF_I0S1[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I0S1.GetCorrFun(ikS));
      CF_I1S0[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I1S0.GetCorrFun(ikS));
      CF_I1S1[iCutOff]->SetBinContent(ikS + 1, AB_pXim_I1S1.GetCorrFun(ikS));

//      for (int iRad = 0; iRad < RadBins; ++iRad) {
//        complex<double> I0S0WF = AB_pXim_I0S0.EvalRadialWaveFunction(ikS,0,0,dRad*iRad,false);
//        CF_I0S0WaveFunction[iCutOff]->SetBinContent(ikS+1,iRad+1,norm(I0S0WF));
//        complex<double> I0S1WF = AB_pXim_I0S1.EvalRadialWaveFunction(ikS,0,0,dRad*iRad,false);
//        CF_I0S1WaveFunction[iCutOff]->SetBinContent(ikS+1,iRad+1,norm(I0S1WF));
//        complex<double> I1S0WF = AB_pXim_I1S0.EvalRadialWaveFunction(ikS,0,0,dRad*iRad,false);
//        CF_I1S0WaveFunction[iCutOff]->SetBinContent(ikS+1,iRad+1,norm(I1S0WF));
//        complex<double> I1S1WF = AB_pXim_I1S1.EvalRadialWaveFunction(ikS,0,0,dRad*iRad,false);
//        CF_I1S1WaveFunction[iCutOff]->SetBinContent(ikS+1,iRad+1,norm(I1S1WF));
//      }
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
    outfile->cd();
    ListI0S0->Add(CF_I0S0[iCutOff]);
//    ListI0S0->Add(CF_I0S0WaveFunction[iCutOff]);

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

    outfile->cd();
    ListI0S1->Add(CF_I0S1[iCutOff]);
//    ListI0S1->Add(CF_I0S1WaveFunction[iCutOff]);

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
    outfile->cd();
    ListI1S0->Add(CF_I1S0[iCutOff]);
//    ListI1S0->Add(CF_I1S0WaveFunction[iCutOff]);

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
    outfile->cd();
    ListI1S1->Add(CF_I1S1[iCutOff]);
//    ListI1S1->Add(CF_I1S1WaveFunction[iCutOff]);

    leg.AddEntry(CF_I0S1[iCutOff], TString::Format("Cut off %.2f", it),
                 "lp");

    iCutOff++;
  }
  outfile->cd();
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
  outfile->cd();
  ListI0S0->Write("I0S0", 1);
  ListI0S1->Write("I0S1", 1);
  ListI1S0->Write("I1S0", 1);
  ListI1S1->Write("I1S1", 1);
  outfile->Close();
  return;
}

int main(int argc, char *argv[]) {

  PlotPotentials();
  return 0;
}
