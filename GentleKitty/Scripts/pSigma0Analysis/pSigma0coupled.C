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
  const float alphaEFT = 0.5;
  const int lineWidth = 3;

  auto grSigma0chiEFT = new TGraph();
  grSigma0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0chiEFT, 20, colorchiEFT, alphaEFT);

  auto grpS0chiEFT = new TGraph();
  grpS0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grpS0chiEFT, 20, colorchiEFT, alphaEFT);
  grpS0chiEFT->SetLineStyle(3);

  auto grpS0nSpluschiEFT = new TGraph();
  grpS0nSpluschiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grpS0nSpluschiEFT, 20, colorchiEFT, alphaEFT);
  grpS0nSpluschiEFT->SetLineStyle(7);


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

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

  const float kmaxdraw = 200;
  const float right = 0.04;
  const float top = 0.025;

  float yminSigma = 0.6;
  const float ymaxSigma = 2.;

  auto dummyHist = new TH1F("dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, kmaxdraw);
  DreamPlot::SetStyleHisto(dummyHist, 20, colorchiEFT);

  auto c = new TCanvas("CFpSigma_smear", "CFpSigma_smear", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grSigma0chiEFT->Draw("L3same");
  grpS0chiEFT->Draw("L3same");
  grpS0nSpluschiEFT->Draw("L3same");
  float xmin = 0.49;
  float xwidth = 0.45;
  float ymin = 0.725;
  float ywidth = 0.225;
  auto leg = new TLegend(xmin, ymin, xmin+xwidth, ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg->AddEntry(grpS0chiEFT, "Genuine p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg->AddEntry(grpS0nSpluschiEFT, "incl. n#minus#kern[-0.95]{ }#Sigma^{+} #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg->AddEntry(grSigma0chiEFT, "incl. p#minus#kern[-0.95]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}", "l");
  leg->Draw("same");


  c->Print("Sigma0_coupled.pdf");

  return 1;
}
