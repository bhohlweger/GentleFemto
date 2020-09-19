#include "TAxis.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "gsl_sf_dawson.h"
#include "DLM_Histo.h"
#include "CATS.h"
#include "DLM_Source.h"
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
void FillWaveGraph(CATS &kitty, TGraph *gr, std::vector<float> weights) {
  int channel = 0;
  for (auto &it : weights) {
    kitty.SetChannelWeight(channel++, it);

  }
  kitty.KillTheCat();
  FillWaveGraph(kitty, gr);
}

/// =====================================================================================
int main(int argc, char* argv[]) {
  DreamPlot::SetStyle();
  double* radius = new double[1];
  radius[0] = 1.25;
  int momBins = 50;
  int kmin = -4.99;
  int kmax = 495.01;

  const int colorchiEFT = kAzure;
  const float alphaEFT = 0.6;
  const int colorNSC97f = kOrange - 3;
  const float alphaNSC97f = 0.8;
  const int lineWidth = 3;

  auto grSigma0chiEFT = new TGraph();
  grSigma0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grSigma0chiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grSigma0chiEFT, 20, colorchiEFT, alphaEFT);

  auto grpS0chiEFT = new TGraph();
  grpS0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpS0chiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpS0chiEFT, 20, colorchiEFT, alphaEFT);
  grpS0chiEFT->SetLineStyle(3);

  auto grpS0nSpluschiEFT = new TGraph();
  grpS0nSpluschiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpS0nSpluschiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpS0nSpluschiEFT, 20, colorchiEFT, alphaEFT);
  grpS0nSpluschiEFT->SetLineStyle(7);


  auto grSigma0NSC97f = new TGraph();
  grSigma0NSC97f->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grSigma0NSC97f->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grSigma0NSC97f, 20, colorNSC97f, alphaNSC97f);

  auto grpS0NSC97f = new TGraph();
  grpS0NSC97f->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpS0NSC97f->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpS0NSC97f, 20, colorNSC97f, alphaNSC97f);
  grpS0NSC97f->SetLineStyle(3);

  auto grpS0nSplusNSC97f = new TGraph();
  grpS0nSplusNSC97f->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpS0nSplusNSC97f->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpS0nSplusNSC97f, 20, colorNSC97f, alphaNSC97f);
  grpS0nSplusNSC97f->SetLineStyle(7);

  /// chiEFT
  TidyCats *tidy = new TidyCats();
  CATS chiEFTKitty;
  tidy->GetCatsProtonSigma0(&chiEFTKitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0Haidenbauer);
  chiEFTKitty.KillTheCat();
  auto Sigma0chiEFT = new DLM_Ck(1, 0, chiEFTKitty);
  Sigma0chiEFT->SetSourcePar(0, radius[0]);
  Sigma0chiEFT->Update();
  FillWaveGraph(chiEFTKitty, grSigma0chiEFT);
  FillWaveGraph(chiEFTKitty, grpS0chiEFT, { 0.25, 0.75, 0, 0, 0, 0 });
  FillWaveGraph(chiEFTKitty, grpS0nSpluschiEFT, { 0.25, 0.75, 0, 0, 0.25, 0.75 });

  /// NSC97f
  CATS NSC97fKitty;
  tidy->GetCatsProtonSigma0(&NSC97fKitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
  NSC97fKitty.KillTheCat();
  auto Sigma0NSC97f = new DLM_Ck(1, 0, NSC97fKitty);
  Sigma0NSC97f->SetSourcePar(0, radius[0]);
  Sigma0NSC97f->Update();
  FillWaveGraph(NSC97fKitty, grSigma0NSC97f);
  FillWaveGraph(NSC97fKitty, grpS0NSC97f, { 0.25, 0.75, 0, 0, 0, 0 });
  FillWaveGraph(NSC97fKitty, grpS0nSplusNSC97f, { 0.25, 0.75, 0, 0, 0.25, 0.75 });

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

  const float kmaxdraw = 200;
  const float right = 0.04;
  const float top = 0.025;

  float yminSigma = 0.6;
  const float ymaxSigma = 2.5;

  auto dummyHist = new TH1F("dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, kmaxdraw);
  DreamPlot::SetStyleHisto(dummyHist, 20, colorchiEFT);

  auto c = new TCanvas("CFpSigma_chiEFT_smear", "CFpSigma_chiEFT_smear", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);
  grSigma0chiEFT->Draw("L3same");
  grpS0chiEFT->Draw("L3same");
  grpS0nSpluschiEFT->Draw("L3same");
  float xmin = 0.4;
  float xwidth = 0.45;
  float ymin = 0.65;
  float ywidth = 0.275;
  auto leg = new TLegend(xmin, ymin, xmin+xwidth, ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg->SetHeader("#chiEFT (NLO)");
  leg->AddEntry(grpS0chiEFT, "Genuine N#minus#kern[-0.95]{ }#Sigma (#it{I} = 1/2, 3/2)", "l");
  leg->AddEntry(grpS0nSpluschiEFT, "incl. n#minus#kern[-0.95]{ }#Sigma^{+} #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg->AddEntry(grSigma0chiEFT, "incl. p#minus#kern[-0.95]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg->Draw("same");
  c->Print("Sigma0_chiEFT_coupled.pdf");

  auto d = new TCanvas("CFpSigma_NSC97f_smear", "CFpSigma_NSC97f_smear", 0, 0, 650, 550);
  d->SetRightMargin(right);
  d->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);
  grSigma0NSC97f->Draw("L3same");
  grpS0NSC97f->Draw("L3same");
  grpS0nSplusNSC97f->Draw("L3same");
  auto leg2 = new TLegend(xmin, ymin, xmin+xwidth, ymin+ywidth);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg2->SetHeader("NSC97f");
  leg2->AddEntry(grpS0NSC97f, "Genuine N#minus#kern[-0.95]{ }#Sigma (#it{I} = 1/2, 3/2)", "l");
  leg2->AddEntry(grpS0nSplusNSC97f, "incl. n#minus#kern[-0.95]{ }#Sigma^{+} #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg2->AddEntry(grSigma0NSC97f, "incl. p#minus#kern[-0.95]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg2->Draw("same");
  d->Print("Sigma0_NSC97f_coupled.pdf");

  return 1;
}
