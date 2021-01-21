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
int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  double *radius = new double[1];
  radius[0] = 1.;
  int momBins = 50;
  int kmin = -4.99;
  int kmax = 495.01;

  const int colorchiEFT = kAzure;
  const float alphaEFT = 0.6;
  const int lineWidth = 3;

  auto grpDminuschiEFT = new TGraph();
  grpDminuschiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpDminuschiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpDminuschiEFT, 20, colorchiEFT, alphaEFT);

  auto grpDminusElastchiEFT = new TGraph();
  grpDminusElastchiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpDminusElastchiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpDminusElastchiEFT, 20, colorchiEFT, alphaEFT);
  grpDminusElastchiEFT->SetLineStyle(3);

  auto grpDminusnDzerochiEFT = new TGraph();
  grpDminusnDzerochiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grpDminusnDzerochiEFT->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grpDminusnDzerochiEFT, 20, colorchiEFT, alphaEFT);
  grpDminusnDzerochiEFT->SetLineStyle(7);

  /// chiEFT
  TidyCats *tidy = new TidyCats();
  CATS chiEFTKitty;
  tidy->GetCatsProtonDminus(&chiEFTKitty, momBins, kmin, kmax,
                            TidyCats::pDminusHaidenbauer, TidyCats::sGaussian);
  chiEFTKitty.SetAnaSource(0, radius[0]);
  chiEFTKitty.KillTheCat();

  FillWaveGraph(chiEFTKitty, grpDminuschiEFT);
  FillWaveGraph(chiEFTKitty, grpDminusElastchiEFT, { 1, 0 });
  FillWaveGraph(chiEFTKitty, grpDminusnDzerochiEFT, { 0, 1 });

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

  const float kmaxdraw = 250;
  const float right = 0.04;
  const float top = 0.025;

  float yminSigma = 0.6;
  const float ymaxSigma = 1.5;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            kmaxdraw);
  DreamPlot::SetStyleHisto(dummyHist, 20, colorchiEFT);

  auto c = new TCanvas("CFpSigma_chiEFT_smear", "CFpSigma_chiEFT_smear", 0, 0,
                       650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);
  grpDminuschiEFT->Draw("L3same");
  grpDminusElastchiEFT->Draw("L3same");
  grpDminusnDzerochiEFT->Draw("L3same");
  float xmin = 0.4;
  float xwidth = 0.45;
  float ymin = 0.65;
  float ywidth = 0.18;
  auto leg = new TLegend(xmin, ymin, xmin + xwidth, ymin + ywidth);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  //leg->SetHeader("#chiEFT (NLO)");
  leg->SetHeader(Form("#it{r}_{0} = %.1f fm", radius[0]));
  leg->AddEntry(grpDminusElastchiEFT, "p#minus#kern[-0.95]{ }D^{-}",
                "l");
  leg->AddEntry(
      grpDminuschiEFT,
      "incl. n#minus#kern[-0.95]{ }D^{0} #rightarrow p#minus#kern[-0.95]{ }D^{-}",
      "l");
  leg->Draw("same");
  c->Print("pDminus_chiEFT_coupled.pdf");

  return 1;
}
