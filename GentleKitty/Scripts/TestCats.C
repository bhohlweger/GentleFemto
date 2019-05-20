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
  CATSparameters* cPars;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, 1.2);

  CATSparameters* cPars2;
  cPars2 = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars2->SetParameter(0, 1.2);


  CATS* AB_pp3P2 = new CATS();
  AB_pp3P2->SetAnaSource(GaussSource, *cPars);
//  AB_pp3P2->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
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
  AB_pp3P2->SetNumPW(0, 1);
  AB_pp3P2->SetSpin(0, 1);
  AB_pp3P2->SetChannelWeight(0, 1);
  AB_pp3P2->SetShortRangePotential(0, 0, fDlmPot, *cPotPars3P2);
  AB_pp3P2->KillTheCat();
  CATS* AB_pp3P2Coupled = new CATS();
  AB_pp3P2Coupled->SetAnaSource(GaussSource, *cPars2);
//  AB_pp3P2Coupled->SetAnaSource(CatsSourceForwarder, fppCleverMcLevy, 2);
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
  AB_pp3P2Coupled->SetNumPW(0, 1);
  AB_pp3P2Coupled->SetSpin(0, 1);
  AB_pp3P2Coupled->SetChannelWeight(0, 1);
  AB_pp3P2Coupled->SetShortRangePotential(0, 0, fDlmPot, *cPotPars3P2Coupled);
  AB_pp3P2Coupled->KillTheCat();
  TGraph Graph3P2;
  TGraph Graph3P2Coupled;
  int dkStar = 2;
  for (int ik = 0; ik < 200; ik++) {
    Graph3P2.SetPoint(ik,dkStar*ik,AB_pp3P2->GetCorrFun(ik));
    Graph3P2Coupled.SetPoint(ik,dkStar*ik,AB_pp3P2Coupled->GetCorrFun(ik));
  }
  auto can = new TCanvas();
  Graph3P2.Draw("APL");
  Graph3P2Coupled.Draw("PLSame");
  TFile* out = new TFile("out.root","recreate");
  out->cd();
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
