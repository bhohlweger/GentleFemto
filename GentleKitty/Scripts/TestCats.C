#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TidyCats.h"
#include "TDatabasePDG.h"

void testCats() {

  TGraph Graph1S03PNoCoupling;
  Graph1S03PNoCoupling.SetLineWidth(3);
  Graph1S03PNoCoupling.SetLineColor(7);

  TGraph Graph1S03PCoupling;
  Graph1S03PCoupling.SetLineWidth(3);
  Graph1S03PCoupling.SetLineColor(8);

  TGraph Graph1S01D23PNoCoupling;
  Graph1S01D23PNoCoupling.SetLineWidth(3);
  Graph1S01D23PNoCoupling.SetLineColor(9);

  TGraph Graph1S01D23PCoupling;
  Graph1S01D23PCoupling.SetLineWidth(3);
  Graph1S01D23PCoupling.SetLineColor(11);

  const double Weight1S0 = 3. / 12.;  // also the weight for 1D2
  const double Weight3P0 = 1. / 12.;
  const double Weight3P1 = 3. / 12.;
  const double Weight3P2 = 5. / 12.;

  double PotPars1S0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
  double PotPars1D2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 2, 2 };  //the last 3 digits are s,l,j
  double PotPars3P0[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
  double PotPars3P1[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
  double PotPars3P2[8] = { NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };
  double PotPars3P2Coupled[8] = { NN_AV18, v18_SingleChannelMagic, 1, 1, 1, 1,
      1, 2 };  // magically accounts for the coupling to F something

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
  
  CATSparameters *cPotPars3P2Coupled = new CATSparameters(
      CATSparameters::tPotential, 8, true);
  cPotPars3P2Coupled->SetParameters(PotPars3P2Coupled);

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;

  CATSparameters* cPars1S03PNoCoupling;
  cPars1S03PNoCoupling = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars1S03PNoCoupling->SetParameter(0, 1.2);

  CATSparameters* cPars1S03PCoupling;
  cPars1S03PCoupling = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars1S03PCoupling->SetParameter(0, 1.2);

  CATSparameters* cPars1S01D23PNoCoupling;
  cPars1S01D23PNoCoupling = new CATSparameters(CATSparameters::tSource, 1,
                                               true);
  cPars1S01D23PNoCoupling->SetParameter(0, 1.2);

  CATSparameters* cPars1S01D23PCoupling;
  cPars1S01D23PCoupling = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars1S01D23PCoupling->SetParameter(0, 1.2);
  int dkStar = 4;
  int kStarBins = 300;
  int kStarMax = dkStar*kStarBins; 
  
  CATS* AB_pp1S03PNoCoupling = new CATS();
  AB_pp1S03PNoCoupling->SetAnaSource(GaussSource, *cPars1S03PNoCoupling);
  AB_pp1S03PNoCoupling->SetAnaSource(0, 1.2);
  AB_pp1S03PNoCoupling->SetAnaSource(1, 2.0);
  AB_pp1S03PNoCoupling->SetEpsilonConv(1e-7);
  AB_pp1S03PNoCoupling->SetEpsilonProp(1e-7);
  AB_pp1S03PNoCoupling->SetUseAnalyticSource(true);
  AB_pp1S03PNoCoupling->SetMomentumDependentSource(false);
  AB_pp1S03PNoCoupling->SetThetaDependentSource(false);
  AB_pp1S03PNoCoupling->SetExcludeFailedBins(false);
  AB_pp1S03PNoCoupling->SetMomBins(kStarBins, 0, kStarMax);
  AB_pp1S03PNoCoupling->SetQ1Q2(1);
  AB_pp1S03PNoCoupling->SetPdgId(2212, 2212);
  AB_pp1S03PNoCoupling->SetRedMass(0.5 * massProton);
  AB_pp1S03PNoCoupling->SetNumChannels(4);
  AB_pp1S03PNoCoupling->SetNumPW(0, 2);
  AB_pp1S03PNoCoupling->SetNumPW(1, 2);
  AB_pp1S03PNoCoupling->SetNumPW(2, 2);
  AB_pp1S03PNoCoupling->SetNumPW(3, 2);
  AB_pp1S03PNoCoupling->SetSpin(0, 0);
  AB_pp1S03PNoCoupling->SetSpin(1, 1);
  AB_pp1S03PNoCoupling->SetSpin(2, 1);
  AB_pp1S03PNoCoupling->SetSpin(3, 1);
  AB_pp1S03PNoCoupling->SetChannelWeight(0, Weight1S0);
  AB_pp1S03PNoCoupling->SetChannelWeight(1, Weight3P0);
  AB_pp1S03PNoCoupling->SetChannelWeight(2, Weight3P1);
  AB_pp1S03PNoCoupling->SetChannelWeight(3, Weight3P2);
  AB_pp1S03PNoCoupling->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp1S03PNoCoupling->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  AB_pp1S03PNoCoupling->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  AB_pp1S03PNoCoupling->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
  AB_pp1S03PNoCoupling->KillTheCat();

  for (int ik = 0; ik < kStarBins; ik++) {
    Graph1S03PNoCoupling.SetPoint(ik, dkStar * ik,
                                  AB_pp1S03PNoCoupling->GetCorrFun(ik));
  }

  delete AB_pp1S03PNoCoupling;
  if (cPars1S03PNoCoupling)
    delete cPars1S03PNoCoupling;

  CATS* AB_pp1S03PCoupling = new CATS();
  AB_pp1S03PCoupling->SetAnaSource(GaussSource, *cPars1S03PCoupling);
  AB_pp1S03PCoupling->SetAnaSource(0, 1.2);
  AB_pp1S03PCoupling->SetAnaSource(1, 2.0);
  AB_pp1S03PCoupling->SetEpsilonConv(1e-7);
  AB_pp1S03PCoupling->SetEpsilonProp(1e-7);
  AB_pp1S03PCoupling->SetUseAnalyticSource(true);
  AB_pp1S03PCoupling->SetMomentumDependentSource(false);
  AB_pp1S03PCoupling->SetThetaDependentSource(false);
  AB_pp1S03PCoupling->SetExcludeFailedBins(false);
  AB_pp1S03PCoupling->SetMomBins(kStarBins, 0, kStarMax);
  AB_pp1S03PCoupling->SetQ1Q2(1);
  AB_pp1S03PCoupling->SetPdgId(2212, 2212);
  AB_pp1S03PCoupling->SetRedMass(0.5 * massProton);
  AB_pp1S03PCoupling->SetNumChannels(4);
  AB_pp1S03PCoupling->SetNumPW(0, 2);
  AB_pp1S03PCoupling->SetNumPW(1, 2);
  AB_pp1S03PCoupling->SetNumPW(2, 2);
  AB_pp1S03PCoupling->SetNumPW(3, 2);
  AB_pp1S03PCoupling->SetSpin(0, 0);
  AB_pp1S03PCoupling->SetSpin(1, 1);
  AB_pp1S03PCoupling->SetSpin(2, 1);
  AB_pp1S03PCoupling->SetSpin(3, 1);
  AB_pp1S03PCoupling->SetChannelWeight(0, Weight1S0);
  AB_pp1S03PCoupling->SetChannelWeight(1, Weight3P0);
  AB_pp1S03PCoupling->SetChannelWeight(2, Weight3P1);
  AB_pp1S03PCoupling->SetChannelWeight(3, Weight3P2);
  AB_pp1S03PCoupling->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp1S03PCoupling->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  AB_pp1S03PCoupling->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  AB_pp1S03PCoupling->SetShortRangePotential(3, 1, fDlmPot,
                                             *cPotPars3P2Coupled);
  AB_pp1S03PCoupling->KillTheCat();

  for (int ik = 0; ik < kStarBins; ik++) {
    Graph1S03PCoupling.SetPoint(ik, dkStar * ik,
                                AB_pp1S03PCoupling->GetCorrFun(ik));
  }

  delete AB_pp1S03PCoupling;
  if (cPars1S03PCoupling)
    delete cPars1S03PCoupling;

  CATS* AB_pp1S01D23PNoCoupling = new CATS();
  AB_pp1S01D23PNoCoupling->SetAnaSource(GaussSource, *cPars1S01D23PNoCoupling);
  AB_pp1S01D23PNoCoupling->SetAnaSource(0, 1.2);
  AB_pp1S01D23PNoCoupling->SetAnaSource(1, 2.0);
  AB_pp1S01D23PNoCoupling->SetEpsilonConv(1e-7);
  AB_pp1S01D23PNoCoupling->SetEpsilonProp(1e-7);
  AB_pp1S01D23PNoCoupling->SetUseAnalyticSource(true);
  AB_pp1S01D23PNoCoupling->SetMomentumDependentSource(false);
  AB_pp1S01D23PNoCoupling->SetThetaDependentSource(false);
  AB_pp1S01D23PNoCoupling->SetExcludeFailedBins(false);
  AB_pp1S01D23PNoCoupling->SetMomBins(kStarBins, 0, kStarMax);
  AB_pp1S01D23PNoCoupling->SetQ1Q2(1);
  AB_pp1S01D23PNoCoupling->SetPdgId(2212, 2212);
  AB_pp1S01D23PNoCoupling->SetRedMass(0.5 * massProton);
  AB_pp1S01D23PNoCoupling->SetNumChannels(4);
  AB_pp1S01D23PNoCoupling->SetNumPW(0, 3);
  AB_pp1S01D23PNoCoupling->SetNumPW(1, 3);
  AB_pp1S01D23PNoCoupling->SetNumPW(2, 3);
  AB_pp1S01D23PNoCoupling->SetNumPW(3, 3);
  AB_pp1S01D23PNoCoupling->SetSpin(0, 0);
  AB_pp1S01D23PNoCoupling->SetSpin(1, 1);
  AB_pp1S01D23PNoCoupling->SetSpin(2, 1);
  AB_pp1S01D23PNoCoupling->SetSpin(3, 1);
  AB_pp1S01D23PNoCoupling->SetChannelWeight(0, Weight1S0);
  AB_pp1S01D23PNoCoupling->SetChannelWeight(1, Weight3P0);
  AB_pp1S01D23PNoCoupling->SetChannelWeight(2, Weight3P1);
  AB_pp1S01D23PNoCoupling->SetChannelWeight(3, Weight3P2);
  AB_pp1S01D23PNoCoupling->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp1S01D23PNoCoupling->SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
  AB_pp1S01D23PNoCoupling->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  AB_pp1S01D23PNoCoupling->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  AB_pp1S01D23PNoCoupling->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
  AB_pp1S01D23PNoCoupling->KillTheCat();

  for (int ik = 0; ik < kStarBins; ik++) {
    Graph1S01D23PNoCoupling.SetPoint(ik, dkStar * ik,
                                     AB_pp1S01D23PNoCoupling->GetCorrFun(ik));
  }

  delete AB_pp1S01D23PNoCoupling;
  if (cPars1S01D23PNoCoupling)
    delete cPars1S01D23PNoCoupling;

  CATS* AB_pp1S01D23PCoupling = new CATS();
  AB_pp1S01D23PCoupling->SetAnaSource(GaussSource, *cPars1S01D23PCoupling);
  AB_pp1S01D23PCoupling->SetAnaSource(0, 1.2);
  AB_pp1S01D23PCoupling->SetAnaSource(1, 2.0);
  AB_pp1S01D23PCoupling->SetEpsilonConv(1e-7);
  AB_pp1S01D23PCoupling->SetEpsilonProp(1e-7);
  AB_pp1S01D23PCoupling->SetUseAnalyticSource(true);
  AB_pp1S01D23PCoupling->SetMomentumDependentSource(false);
  AB_pp1S01D23PCoupling->SetThetaDependentSource(false);
  AB_pp1S01D23PCoupling->SetExcludeFailedBins(false);
  AB_pp1S01D23PCoupling->SetMomBins(kStarBins, 0, kStarMax);
  AB_pp1S01D23PCoupling->SetQ1Q2(1);
  AB_pp1S01D23PCoupling->SetPdgId(2212, 2212);
  AB_pp1S01D23PCoupling->SetRedMass(0.5 * massProton);
  AB_pp1S01D23PCoupling->SetNumChannels(4);
  AB_pp1S01D23PCoupling->SetNumPW(0, 3);
  AB_pp1S01D23PCoupling->SetNumPW(1, 3);
  AB_pp1S01D23PCoupling->SetNumPW(2, 3);
  AB_pp1S01D23PCoupling->SetNumPW(3, 3);
  AB_pp1S01D23PCoupling->SetSpin(0, 0);
  AB_pp1S01D23PCoupling->SetSpin(1, 1);
  AB_pp1S01D23PCoupling->SetSpin(2, 1);
  AB_pp1S01D23PCoupling->SetSpin(3, 1);
  AB_pp1S01D23PCoupling->SetChannelWeight(0, Weight1S0);
  AB_pp1S01D23PCoupling->SetChannelWeight(1, Weight3P0);
  AB_pp1S01D23PCoupling->SetChannelWeight(2, Weight3P1);
  AB_pp1S01D23PCoupling->SetChannelWeight(3, Weight3P2);
  AB_pp1S01D23PCoupling->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
  AB_pp1S01D23PCoupling->SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
  AB_pp1S01D23PCoupling->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
  AB_pp1S01D23PCoupling->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
  AB_pp1S01D23PCoupling->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2Coupled);
  AB_pp1S01D23PCoupling->KillTheCat();

  for (int ik = 0; ik < kStarBins; ik++) {
    Graph1S01D23PCoupling.SetPoint(ik, dkStar * ik,
                                     AB_pp1S01D23PCoupling->GetCorrFun(ik));
  }

  delete AB_pp1S01D23PCoupling;
  if (cPars1S01D23PCoupling)
    delete cPars1S01D23PCoupling;
  
  TFile* out = new TFile("out.root", "recreate");
  out->cd();

  Graph1S03PNoCoupling.SetName("1S03P2NoCoupling");
  Graph1S03PNoCoupling.Write();

  Graph1S03PCoupling.SetName("1S03P2Coupling");
  Graph1S03PCoupling.Write();

  Graph1S01D23PNoCoupling.SetName("1S01D23P2NoCoupling");
  Graph1S01D23PNoCoupling.Write();

  Graph1S01D23PCoupling.SetName("1S01D23P2Coupling");
  Graph1S01D23PCoupling.Write();

  out->Write();
  out->Close();
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}
