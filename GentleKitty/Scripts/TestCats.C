#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TidyCats.h"
#include "TDatabasePDG.h"

void testCats() {

  double PotPars1S0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  double PotPars1D2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 2, 2 };  //the last 3 digits are s,l,j
  double PotPars3P0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  double PotPars3P1[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
  double PotPars3P2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };
  double PotPars3P2Coupled[8] = { NN_AV18, v18_SingleChannelMagic, 1, 1, 1, 1, 1, 2 };  // magically accounts for the coupling to F something

  CATSparameters *cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars1S0->SetParameters(PotPars1S0);
  CATSparameters *cPotPars1D2 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars1D2->SetParameters(PotPars1D2);
  CATSparameters *cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P0->SetParameters(PotPars3P0);
  CATSparameters *cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P1->SetParameters(PotPars3P1);
  CATSparameters *cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P2->SetParameters(PotPars3P2);
  CATSparameters *cPotPars3P2Coupled = new CATSparameters(CATSparameters::tPotential,
                                                   8, true);
  cPotPars3P2Coupled->SetParameters(PotPars3P2Coupled);
  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  CATSparameters* cPars1S0;
  cPars1S0 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars1S0->SetParameter(0, 1.2);

  CATSparameters* cPars1D2;
  cPars1D2 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars1D2->SetParameter(0, 1.2);

  CATSparameters* cPars3P0;
  cPars3P0 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars3P0->SetParameter(0, 1.2);

  CATSparameters* cPars3P1;
  cPars3P1 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars3P1->SetParameter(0, 1.2);

  CATSparameters* cPars3P2;
  cPars3P2 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars3P2->SetParameter(0, 1.2);

  CATSparameters* cPars3P2Coupled;
  cPars3P2Coupled = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars3P2Coupled->SetParameter(0, 1.2);

  CATS* AB_pp1S0 = new CATS();
  AB_pp1S0->SetAnaSource(GaussSource, *cPars3P2);
  AB_pp1S0->SetAnaSource(0, 1.2);
  AB_pp1S0->SetAnaSource(1, 2.0);
  AB_pp1S0->SetEpsilonConv(1e-8);
  AB_pp1S0->SetEpsilonProp(1e-8);
  AB_pp1S0->SetUseAnalyticSource(true);
  AB_pp1S0->SetMomentumDependentSource(false);
  AB_pp1S0->SetThetaDependentSource(false);
  AB_pp1S0->SetExcludeFailedBins(false);
  AB_pp1S0->SetMomBins(200, 0, 400);
  AB_pp1S0->SetQ1Q2(1);
  AB_pp1S0->SetPdgId(2212, 2212);
  AB_pp1S0->SetRedMass(0.5 * massProton);
  AB_pp1S0->SetNumChannels(1);
  AB_pp1S0->SetNumPW(0, 3);
  AB_pp1S0->SetSpin(0, 0);
  AB_pp1S0->SetChannelWeight(0, 1);
  AB_pp1S0->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp1S0->KillTheCat();

  CATS* AB_pp1D2 = new CATS();
  AB_pp1D2->SetAnaSource(GaussSource, *cPars3P2);
  AB_pp1D2->SetAnaSource(0, 1.2);
  AB_pp1D2->SetAnaSource(1, 2.0);
  AB_pp1D2->SetEpsilonConv(1e-8);
  AB_pp1D2->SetEpsilonProp(1e-8);
  AB_pp1D2->SetUseAnalyticSource(true);
  AB_pp1D2->SetMomentumDependentSource(false);
  AB_pp1D2->SetThetaDependentSource(false);
  AB_pp1D2->SetExcludeFailedBins(false);
  AB_pp1D2->SetMomBins(200, 0, 400);
  AB_pp1D2->SetQ1Q2(1);
  AB_pp1D2->SetPdgId(2212, 2212);
  AB_pp1D2->SetRedMass(0.5 * massProton);
  AB_pp1D2->SetNumChannels(1);
  AB_pp1D2->SetNumPW(0, 3);
  AB_pp1D2->SetSpin(0, 0);
  AB_pp1D2->SetChannelWeight(0, 1);
  AB_pp1D2->SetShortRangePotential(0, 2, fDlmPot, *cPotPars3P2);
  AB_pp1D2->KillTheCat();

  CATS* AB_pp3P0 = new CATS();
  AB_pp3P0->SetAnaSource(GaussSource, *cPars3P2);
  AB_pp3P0->SetAnaSource(0, 1.2);
  AB_pp3P0->SetAnaSource(1, 2.0);
  AB_pp3P0->SetEpsilonConv(1e-8);
  AB_pp3P0->SetEpsilonProp(1e-8);
  AB_pp3P0->SetUseAnalyticSource(true);
  AB_pp3P0->SetMomentumDependentSource(false);
  AB_pp3P0->SetThetaDependentSource(false);
  AB_pp3P0->SetExcludeFailedBins(false);
  AB_pp3P0->SetMomBins(200, 0, 400);
  AB_pp3P0->SetQ1Q2(1);
  AB_pp3P0->SetPdgId(2212, 2212);
  AB_pp3P0->SetRedMass(0.5 * massProton);
  AB_pp3P0->SetNumChannels(1);
  AB_pp3P0->SetNumPW(0, 3);
  AB_pp3P0->SetSpin(0, 1);
  AB_pp3P0->SetChannelWeight(0, 1);
  AB_pp3P0->SetShortRangePotential(0, 1, fDlmPot, *cPotPars3P0);
  AB_pp3P0->KillTheCat();

  CATS* AB_pp3P1 = new CATS();
  AB_pp3P1->SetAnaSource(GaussSource, *cPars3P2);
  AB_pp3P1->SetAnaSource(0, 1.2);
  AB_pp3P1->SetAnaSource(1, 2.0);
  AB_pp3P1->SetEpsilonConv(1e-8);
  AB_pp3P1->SetEpsilonProp(1e-8);
  AB_pp3P1->SetUseAnalyticSource(true);
  AB_pp3P1->SetMomentumDependentSource(false);
  AB_pp3P1->SetThetaDependentSource(false);
  AB_pp3P1->SetExcludeFailedBins(false);
  AB_pp3P1->SetMomBins(200, 0, 400);
  AB_pp3P1->SetQ1Q2(1);
  AB_pp3P1->SetPdgId(2212, 2212);
  AB_pp3P1->SetRedMass(0.5 * massProton);
  AB_pp3P1->SetNumChannels(1);
  AB_pp3P1->SetNumPW(0, 3);
  AB_pp3P1->SetSpin(0, 1);
  AB_pp3P1->SetChannelWeight(0, 1);
  AB_pp3P1->SetShortRangePotential(0, 1, fDlmPot, *cPotPars3P1);
  AB_pp3P1->KillTheCat();

  CATS* AB_pp3P2 = new CATS();
  AB_pp3P2->SetAnaSource(GaussSource, *cPars3P2);
  AB_pp3P2->SetAnaSource(0, 1.2);
  AB_pp3P2->SetAnaSource(1, 2.0);
  AB_pp3P2->SetEpsilonConv(1e-8);
  AB_pp3P2->SetEpsilonProp(1e-8);
  AB_pp3P2->SetUseAnalyticSource(true);
  AB_pp3P2->SetMomentumDependentSource(false);
  AB_pp3P2->SetThetaDependentSource(false);
  AB_pp3P2->SetExcludeFailedBins(false);
  AB_pp3P2->SetMomBins(200, 0, 400);
  AB_pp3P2->SetQ1Q2(1);
  AB_pp3P2->SetPdgId(2212, 2212);
  AB_pp3P2->SetRedMass(0.5 * massProton);
  AB_pp3P2->SetNumChannels(1);
  AB_pp3P2->SetNumPW(0, 3);
  AB_pp3P2->SetSpin(0, 1);
  AB_pp3P2->SetChannelWeight(0, 1);
  AB_pp3P2->SetShortRangePotential(0, 1, fDlmPot, *cPotPars3P2);
  AB_pp3P2->KillTheCat();

  CATS* AB_pp3P2Coupled = new CATS();
  AB_pp3P2Coupled->SetAnaSource(GaussSource, *cPars3P2Coupled);
  AB_pp3P2Coupled->SetAnaSource(0, 1.2);
  AB_pp3P2Coupled->SetAnaSource(1, 2.0);
  AB_pp3P2Coupled->SetEpsilonConv(1e-8);
  AB_pp3P2Coupled->SetEpsilonProp(1e-8);
  AB_pp3P2Coupled->SetUseAnalyticSource(true);
  AB_pp3P2Coupled->SetMomentumDependentSource(false);
  AB_pp3P2Coupled->SetThetaDependentSource(false);
  AB_pp3P2Coupled->SetExcludeFailedBins(false);
  AB_pp3P2Coupled->SetMomBins(200, 0, 400);
  AB_pp3P2Coupled->SetQ1Q2(1);
  AB_pp3P2Coupled->SetPdgId(2212, 2212);
  AB_pp3P2Coupled->SetRedMass(0.5 * massProton);
  AB_pp3P2Coupled->SetNumChannels(1);
  AB_pp3P2Coupled->SetNumPW(0, 3);
  AB_pp3P2Coupled->SetSpin(0, 1);
  AB_pp3P2Coupled->SetChannelWeight(0, 1);
  AB_pp3P2Coupled->SetShortRangePotential(0, 1, fDlmPot, *cPotPars3P2Coupled);
  AB_pp3P2Coupled->KillTheCat();

  TGraph Graph1S0;
  Graph1S0.SetLineWidth(3);
  Graph1S0.SetLineColor(1);
  
  TGraph Graph1D2;
  Graph1D2.SetLineWidth(3);
  Graph1D2.SetLineColor(2);
  
  TGraph Graph3P0;
  Graph3P0.SetLineWidth(3);
  Graph3P0.SetLineColor(3);
  
  TGraph Graph3P1;
  Graph3P1.SetLineWidth(3);
  Graph3P1.SetLineColor(4);
  
  TGraph Graph3P2;
  Graph3P2.SetLineWidth(3);
  Graph3P2.SetLineColor(5);
  
  TGraph Graph3P2Coupled;
  Graph3P2Coupled.SetLineWidth(3);
  Graph3P2Coupled.SetLineColor(6);
  
  int dkStar = 2;
  for (int ik = 0; ik < 200; ik++) {
    Graph1S0.SetPoint(ik,dkStar*ik,AB_pp1S0->GetCorrFun(ik));
    Graph1D2.SetPoint(ik,dkStar*ik,AB_pp1D2->GetCorrFun(ik));
    Graph3P0.SetPoint(ik,dkStar*ik,AB_pp3P0->GetCorrFun(ik));
    Graph3P1.SetPoint(ik,dkStar*ik,AB_pp3P1->GetCorrFun(ik));
    Graph3P2.SetPoint(ik,dkStar*ik,AB_pp3P2->GetCorrFun(ik));
    Graph3P2Coupled.SetPoint(ik,dkStar*ik,AB_pp3P2Coupled->GetCorrFun(ik));
  }

  auto can = new TCanvas();
  Graph1S0.Draw("APL");
  Graph1D2.Draw("PLSame");
  Graph3P0.Draw("PLSame");
  Graph3P1.Draw("PLSame");
  Graph3P2.Draw("PLSame");
  Graph3P2Coupled.Draw("PLSame");
  TFile* out = new TFile("out.root","recreate");
  out->cd();
  Graph1S0.SetName("1S0");
  Graph1S0.Write();

  Graph1D2.SetName("1D2");
  Graph1D2.Write();

  Graph3P0.SetName("3P0");
  Graph3P0.Write();

  Graph3P1.SetName("3P1");
  Graph3P1.Write();

  Graph3P2.SetName("3P2");
  Graph3P2.Write();

  Graph3P2Coupled.SetName("3P2Coupled");
  Graph3P2Coupled.Write();

  can->Write();
  out->Write();
  out->Close();
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}
