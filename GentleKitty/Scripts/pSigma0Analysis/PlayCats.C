#include "TAxis.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "gsl_sf_dawson.h"
#include "DLM_Histo.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "TLatex.h"
#include "TApplication.h"
#include "TGraph.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "DreamPlot.h"
#include "CATSInputSigma0.h"
#include "DLM_CkDecomposition.h"
#include "CATSLambdaParam.h"
#include "TDatabasePDG.h"
#include "TidyCats.h"
#include "TStyle.h"

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
int main(int argc, char* argv[]) {

  DreamPlot::SetStyle();

  TidyCats tidy;

  TApplication *app = new TApplication("app", 0, 0);
  const double radius1 = 1.25;
  const double radius2 = 4;

  auto sourceSmall = new CATSparameters(CATSparameters::tSource, 1, true);
  sourceSmall->SetParameter(0, radius1);

  auto sourceLarge = new CATSparameters(CATSparameters::tSource, 1, true);
  sourceLarge->SetParameter(0, radius2);

  const int NumMomBins = 600;
  const double kMin = 0;
  const double kMax = 300;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            kMax);
  DreamPlot::SetStyleHisto(dummyHist, 20, kWhite);
  dummyHist->GetYaxis()->SetRangeUser(0., 3);
  dummyHist->GetYaxis()->SetTitleOffset(1.15);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);

  const int lineWidth = 3;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Identical Boson
  auto grIdenticalBosonRadius1 = new TGraph();
  grIdenticalBosonRadius1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grIdenticalBosonRadius1->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grIdenticalBosonRadius1, 20, kBlue + 3, 0.8);

  auto grIdenticalBosonRadius2 = new TGraph();
  grIdenticalBosonRadius2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grIdenticalBosonRadius2->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grIdenticalBosonRadius2, 20, kBlue + 3, 0.8);
  grIdenticalBosonRadius2->SetLineStyle(2);

  CATS IdenticalBosonKitty;
  IdenticalBosonKitty.SetUseAnalyticSource(true);
  IdenticalBosonKitty.SetMomBins(NumMomBins, kMin, kMax);
  IdenticalBosonKitty.SetNumChannels(1);
  IdenticalBosonKitty.SetNumPW(0, 0);
  IdenticalBosonKitty.SetSpin(0, 0);
  IdenticalBosonKitty.SetChannelWeight(0, 1);
  IdenticalBosonKitty.SetQ1Q2(0);
  IdenticalBosonKitty.SetPdgId(211, 211);
  IdenticalBosonKitty.SetQuantumStatistics(true);
  IdenticalBosonKitty.SetRedMass(0.5 * TDatabasePDG::Instance()->GetParticle(211)->Mass()* 1000);

  IdenticalBosonKitty.SetAnaSource(GaussSource, *sourceSmall);
  IdenticalBosonKitty.KillTheCat();
  FillWaveGraph(IdenticalBosonKitty, grIdenticalBosonRadius1);

  IdenticalBosonKitty.SetAnaSource(GaussSource, *sourceLarge);
  IdenticalBosonKitty.KillTheCat();
  FillWaveGraph(IdenticalBosonKitty, grIdenticalBosonRadius2);

  auto c = new TCanvas("IdenticalBoson", "IdenticalBoson", 0, 0, 650, 550);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.025);
  dummyHist->Draw();
  grIdenticalBosonRadius1->Draw("L3");
  grIdenticalBosonRadius2->Draw("L3");

  auto leg1 = new TLegend(0.55, 0.7, 0.75, 0.91);
  leg1->SetBorderSize(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg1->SetHeader("Identical bosons");
  leg1->AddEntry(grIdenticalBosonRadius1, Form("#it{r}_{0} = %.2f fm", radius1), "l");
  leg1->AddEntry(grIdenticalBosonRadius2, Form("#it{r}_{0} = %.2f fm", radius2),
                 "l");
  leg1->Draw("same");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Identical fermion
  auto grIdenticalFermionRadius1 = new TGraph();
  grIdenticalFermionRadius1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grIdenticalFermionRadius1->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grIdenticalFermionRadius1, 20, kGreen + 3, 0.8);

  auto grIdenticalFermionRadius2 = new TGraph();
  grIdenticalFermionRadius2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grIdenticalFermionRadius2->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grIdenticalFermionRadius2, 20, kGreen + 3, 0.8);
  grIdenticalFermionRadius2->SetLineStyle(2);

  CATS IdenticalFermionKitty;
  IdenticalFermionKitty.SetMomBins(NumMomBins, kMin, kMax);
  IdenticalFermionKitty.SetUseAnalyticSource(true);
  IdenticalFermionKitty.SetQ1Q2(0);
  IdenticalFermionKitty.SetPdgId(2212, 2212);
  IdenticalFermionKitty.SetRedMass(
      0.5 * TDatabasePDG::Instance()->GetParticle(2212)->Mass()* 1000);
  IdenticalFermionKitty.SetNumChannels(4);
  IdenticalFermionKitty.SetNumPW(0, 3);
  IdenticalFermionKitty.SetNumPW(1, 3);
  IdenticalFermionKitty.SetNumPW(2, 3);
  IdenticalFermionKitty.SetNumPW(3, 3);
  IdenticalFermionKitty.SetSpin(0, 0);
  IdenticalFermionKitty.SetSpin(1, 1);
  IdenticalFermionKitty.SetSpin(2, 1);
  IdenticalFermionKitty.SetSpin(3, 1);
  IdenticalFermionKitty.SetChannelWeight(0, 3. / 12.);
  IdenticalFermionKitty.SetChannelWeight(1, 1. / 12.);
  IdenticalFermionKitty.SetChannelWeight(2, 3. / 12.);
  IdenticalFermionKitty.SetChannelWeight(3, 5. / 12.);
  IdenticalFermionKitty.SetQuantumStatistics(true);

  IdenticalFermionKitty.SetAnaSource(GaussSource, *sourceSmall);
  IdenticalFermionKitty.KillTheCat();
  FillWaveGraph(IdenticalFermionKitty, grIdenticalFermionRadius1);

  IdenticalFermionKitty.SetAnaSource(GaussSource, *sourceLarge);
  IdenticalFermionKitty.KillTheCat();
  FillWaveGraph(IdenticalFermionKitty, grIdenticalFermionRadius2);

  auto d = new TCanvas("IdenticalFermion", "IdenticalFermion", 0, 0, 650, 550);
  d->SetRightMargin(0.04);
  d->SetTopMargin(0.025);
  dummyHist->Draw();
  grIdenticalFermionRadius1->Draw("L3");
  grIdenticalFermionRadius2->Draw("L3");

  auto leg2 = new TLegend(0.55, 0.7, 0.75, 0.91);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg2->SetHeader("Identical fermions");
  leg2->AddEntry(grIdenticalFermionRadius1,
                 Form("#it{r}_{0} = %.2f fm", radius1), "l");
  leg2->AddEntry(grIdenticalFermionRadius2,
                 Form("#it{r}_{0} = %.2f fm", radius2), "l");
  leg2->Draw("same");


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Attractive Coulomb only

  auto grCoulombAttractionRadius1 = new TGraph();
  grCoulombAttractionRadius1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grCoulombAttractionRadius1->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grCoulombAttractionRadius1, 20, kAzure, 0.8);

  auto grCoulombAttractionRadius2 = new TGraph();
  grCoulombAttractionRadius2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grCoulombAttractionRadius2->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grCoulombAttractionRadius2, 20, kAzure, 0.8);
  grCoulombAttractionRadius2->SetLineStyle(2);

  CATS CoulombKitty;
  CoulombKitty.SetUseAnalyticSource(true);
  CoulombKitty.SetMomBins(NumMomBins, kMin, kMax);
  CoulombKitty.SetNumChannels(1);
  CoulombKitty.SetNumPW(0, 0);
  CoulombKitty.SetSpin(0, 0);
  CoulombKitty.SetChannelWeight(0, 1);
  CoulombKitty.SetQ1Q2(-1);
  CoulombKitty.SetPdgId(211, 211);
  CoulombKitty.SetQuantumStatistics(false);
  CoulombKitty.SetRedMass(0.5 * TDatabasePDG::Instance()->GetParticle(211)->Mass()* 1000);

  CoulombKitty.SetAnaSource(GaussSource, *sourceSmall);
  CoulombKitty.KillTheCat();
  FillWaveGraph(CoulombKitty, grCoulombAttractionRadius1);

  CoulombKitty.SetAnaSource(GaussSource, *sourceLarge);
  CoulombKitty.KillTheCat();
  FillWaveGraph(CoulombKitty, grCoulombAttractionRadius2);

  auto e = new TCanvas("Attractive Coulomb", "Attractive Coulomb", 0, 0, 650, 550);
  e->SetRightMargin(0.04);
  e->SetTopMargin(0.025);
  dummyHist->Draw();
  grCoulombAttractionRadius1->Draw("L3");
  grCoulombAttractionRadius2->Draw("L3");

  auto leg3 = new TLegend(0.55, 0.7, 0.75, 0.91);
  leg3->SetBorderSize(0);
  leg3->SetTextFont(42);
  leg3->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg3->SetHeader("Attractive Coulomb");
  leg3->AddEntry(grCoulombAttractionRadius1,
                 Form("#it{r}_{0} = %.2f fm", radius1), "l");
  leg3->AddEntry(grCoulombAttractionRadius2,
                 Form("#it{r}_{0} = %.2f fm", radius2), "l");
  leg3->Draw("same");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Repulsive Coulomb only

  auto grCoulombRepulsionRadius1 = new TGraph();
  grCoulombRepulsionRadius1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grCoulombRepulsionRadius1->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grCoulombRepulsionRadius1, 20, kOrange + 2, 0.8);

  auto grCoulombRepulsionRadius2 = new TGraph();
  grCoulombRepulsionRadius2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grCoulombRepulsionRadius2->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grCoulombRepulsionRadius2, 20, kOrange + 2, 0.8);
  grCoulombRepulsionRadius2->SetLineStyle(2);

  CoulombKitty.SetQ1Q2(1);

  CoulombKitty.SetAnaSource(GaussSource, *sourceSmall);
  CoulombKitty.KillTheCat();
  FillWaveGraph(CoulombKitty, grCoulombRepulsionRadius1);

  CoulombKitty.SetAnaSource(GaussSource, *sourceLarge);
  CoulombKitty.KillTheCat();
  FillWaveGraph(CoulombKitty, grCoulombRepulsionRadius2);

  auto f = new TCanvas("Repulsive Coulomb", "Repulsive Coulomb", 0, 0, 650, 550);
  f->SetRightMargin(0.04);
  f->SetTopMargin(0.025);
  dummyHist->Draw();
  grCoulombRepulsionRadius1->Draw("L3");
  grCoulombRepulsionRadius2->Draw("L3");

  auto leg4 = new TLegend(0.55, 0.7, 0.75, 0.91);
  leg4->SetBorderSize(0);
  leg4->SetTextFont(42);
  leg4->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg4->SetHeader("Repulsive Coulomb");
  leg4->AddEntry(grCoulombRepulsionRadius1,
                 Form("#it{r}_{0} = %.2f fm", radius1), "l");
  leg4->AddEntry(grCoulombRepulsionRadius2,
                 Form("#it{r}_{0} = %.2f fm", radius2), "l");
  leg4->Draw("same");


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Summary plots
  auto g = new TCanvas("Identical Particles", "Identical Particles", 0, 0, 650, 550);
  g->SetRightMargin(0.04);
  g->SetTopMargin(0.025);
  dummyHist->Draw();
  grIdenticalBosonRadius1->Draw("L3");
  grIdenticalBosonRadius2->Draw("L3");
  grIdenticalFermionRadius1->Draw("L3");
  grIdenticalFermionRadius2->Draw("L3");

  auto legg1 = new TLegend(0.45, 0.63, 0.65, 0.84);
  legg1->SetBorderSize(0);
  legg1->SetTextFont(42);
  legg1->SetTextSize(gStyle->GetTextSize() * 0.9);
  legg1->SetHeader(Form("#it{r}_{0} = %.2f fm", radius1));
  legg1->AddEntry(grIdenticalBosonRadius1, "Boson", "l");
  legg1->AddEntry(grIdenticalFermionRadius1, "Fermion", "l");
  legg1->Draw("same");

  auto legg2 = new TLegend(0.7, 0.63, 0.9, 0.84);
  legg2->SetBorderSize(0);
  legg2->SetTextFont(42);
  legg2->SetTextSize(gStyle->GetTextSize() * 0.9);
  legg2->SetHeader(Form("#it{r}_{0} = %.2f fm", radius2));
  legg2->AddEntry(grIdenticalBosonRadius2, "Boson", "l");
  legg2->AddEntry(grIdenticalFermionRadius2, "Fermion", "l");
  legg2->Draw("same");

  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(gStyle->GetTextSize() * 0.9);
  text.DrawLatex(0.45, 0.855, "#bf{Quantum Statistics}");

  g->Print("Identicals.pdf");

  auto h = new TCanvas("Coulomb", "Coulomb", 0, 0, 650, 550);
  h->SetRightMargin(0.04);
  h->SetTopMargin(0.025);
  dummyHist->Draw();
  grCoulombAttractionRadius1->Draw("L3");
  grCoulombAttractionRadius2->Draw("L3");
  grCoulombRepulsionRadius1->Draw("L3");
  grCoulombRepulsionRadius2->Draw("L3");

  auto legh1 = new TLegend(0.45, 0.63, 0.65, 0.84);
  legh1->SetBorderSize(0);
  legh1->SetTextFont(42);
  legh1->SetTextSize(gStyle->GetTextSize() * 0.9);
  legh1->SetHeader(Form("#it{r}_{0} = %.2f fm", radius1));
  legh1->AddEntry(grCoulombAttractionRadius1, "Attraction", "l");
  legh1->AddEntry(grCoulombRepulsionRadius1, "Repulsion", "l");
  legh1->Draw("same");

  auto legh2 = new TLegend(0.7, 0.63, 0.9, 0.84);
  legh2->SetBorderSize(0);
  legh2->SetTextFont(42);
  legh2->SetTextSize(gStyle->GetTextSize() * 0.9);
  legh2->SetHeader(Form("#it{r}_{0} = %.2f fm", radius2));
  legh2->AddEntry(grCoulombAttractionRadius2, "Attraction", "l");
  legh2->AddEntry(grCoulombRepulsionRadius2, "Repulsion", "l");
  legh2->Draw("same");

  text.DrawLatex(0.45, 0.855, "#bf{Coulomb}");

  h->Print("Coulomb.pdf");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Proton proton

  auto grPPFull = new TGraph();
  grPPFull->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPPFull->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPPFull, 20, kTeal + 2, 0.8);

  auto grPPQS = new TGraph();
  grPPQS->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPPQS->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPPQS, 20, kTeal + 2, 0.8);
  grPPQS->SetLineStyle(2);

  auto grPPCoulomb = new TGraph();
  grPPCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPPCoulomb->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPPCoulomb, 20, kTeal + 2, 0.8);
  grPPCoulomb->SetLineStyle(3);

  auto grPPStrong = new TGraph();
  grPPStrong->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPPStrong->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPPStrong, 20, kTeal + 2, 0.8);
  grPPStrong->SetLineStyle(4);

  CATS ppKitty;
  tidy.GetCatsProtonProton(&ppKitty, NumMomBins, kMin, kMax, TidyCats::sGaussian);
  ppKitty.SetAnaSource(GaussSource, *sourceSmall);
  ppKitty.KillTheCat();
  FillWaveGraph(ppKitty, grPPFull);

  ppKitty.SetQ1Q2(0);
  ppKitty.SetQuantumStatistics(false);
  ppKitty.KillTheCat();
  FillWaveGraph(ppKitty, grPPStrong);

  IdenticalFermionKitty.SetAnaSource(GaussSource, *sourceSmall);
  IdenticalFermionKitty.KillTheCat();
  FillWaveGraph(IdenticalFermionKitty, grPPQS);

  IdenticalFermionKitty.SetQ1Q2(1);
  IdenticalFermionKitty.SetQuantumStatistics(false);
  IdenticalFermionKitty.KillTheCat();
  FillWaveGraph(IdenticalFermionKitty, grPPCoulomb);

  auto pp = new TCanvas("pp", "pp", 0, 0, 650, 550);
  pp->SetRightMargin(0.04);
  pp->SetTopMargin(0.025);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0., 5);
  dummyHist->GetYaxis()->SetNdivisions(505);
  grPPFull->Draw("L3");
  grPPQS->Draw("L3");
  grPPCoulomb->Draw("L3");
  grPPStrong->Draw("L3");

  auto legpp = new TLegend(0.45, 0.49, 0.75, 0.84);
  legpp->SetBorderSize(0);
  legpp->SetTextFont(42);
  legpp->SetTextSize(gStyle->GetTextSize() * 0.9);
  legpp->SetHeader(Form("#it{r}_{0} = %.2f fm", radius1));
  legpp->AddEntry(grPPFull, "Full", "l");
  legpp->AddEntry(grPPQS, "Quantum Statistics", "l");
  legpp->AddEntry(grPPCoulomb, "Coulomb Interaction", "l");
  legpp->AddEntry(grPPStrong, "Strong Interaction", "l");
  legpp->Draw("same");

  text.DrawLatex(0.45, 0.855, "p#minus#kern[-0.95]{ }p");

  pp->Print("ppFull.pdf");


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Pion Pion

  auto grPipiFull = new TGraph();
  grPipiFull->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPipiFull->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPipiFull, 20, kCyan -1, 0.8);

  auto grPipiQS = new TGraph();
  grPipiQS->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPipiQS->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPipiQS, 20, kCyan -1, 0.8);
  grPipiQS->SetLineStyle(2);

  auto grPipiCoulomb = new TGraph();
  grPipiCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grPipiCoulomb->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grPipiCoulomb, 20, kCyan -1, 0.8);
  grPipiCoulomb->SetLineStyle(3);

  IdenticalBosonKitty.SetAnaSource(GaussSource, *sourceSmall);
  IdenticalBosonKitty.SetQ1Q2(1);
  IdenticalBosonKitty.KillTheCat();
  FillWaveGraph(IdenticalBosonKitty, grPipiFull);

  IdenticalBosonKitty.SetQ1Q2(0);
  IdenticalBosonKitty.KillTheCat();
  FillWaveGraph(IdenticalBosonKitty, grPipiQS);

  IdenticalBosonKitty.SetQ1Q2(1);
  IdenticalBosonKitty.SetQuantumStatistics(false);
  IdenticalBosonKitty.KillTheCat();
  FillWaveGraph(IdenticalBosonKitty, grPipiCoulomb);

  auto pipi = new TCanvas("pipi", "pipi", 0, 0, 650, 550);
  pipi->SetRightMargin(0.04);
  pipi->SetTopMargin(0.025);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0., 3);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grPipiFull->Draw("L3");
  grPipiQS->Draw("L3");
  grPipiCoulomb->Draw("L3");

  auto legpipi = new TLegend(0.45, 0.56, 0.75, 0.84);
  legpipi->SetBorderSize(0);
  legpipi->SetTextFont(42);
  legpipi->SetTextSize(gStyle->GetTextSize() * 0.9);
  legpipi->SetHeader(Form("#it{r}_{0} = %.2f fm", radius1));
  legpipi->AddEntry(grPipiFull, "Full", "l");
  legpipi->AddEntry(grPipiQS, "Quantum Statistics", "l");
  legpipi->AddEntry(grPipiCoulomb, "Coulomb Interaction", "l");
  legpipi->Draw("same");

  text.DrawLatex(0.45, 0.855, "#pi^{+}#minus#kern[-0.5]{ }#pi^{+}");

  pipi->Print("pipiFull.pdf");

  app->Run();
}
