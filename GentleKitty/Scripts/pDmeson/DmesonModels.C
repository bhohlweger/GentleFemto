#include "TSystem.h"
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
int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  const double radius = 0.9;
  int momBins = 150;
  int kmin = 1;
  int kmax = 301;

  const int colHaide = kOrange;
  const int colCoulomb = kGreen + 2;
  const int colModel1 = kBlue + 3;
  const int colModel2 = kOrange;
  const int colModel3 = kRed + 3;
  const int colModel4 = kTeal;
  const int lineWidth = 3;

  auto grHaidenbauer = new TGraph();
  grHaidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grHaidenbauer->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grHaidenbauer, 20, colHaide, 0.8);

  auto grCoulomb = new TGraph();
  grCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grCoulomb->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grCoulomb, 20, colCoulomb, 0.8);

  TString HomeDir = gSystem->GetHomeDirectory().c_str();
  auto grModel1 = new TGraph(
      TString::Format(
          "%s/CERNHome/D-mesons/Analysis/Models/corr_model1_0.9fm_wC.dat",
          HomeDir.Data()));
  grModel1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grModel1->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grModel1, 20, colModel1, 0.8);

  auto grModel2 = new TGraph(
      TString::Format(
          "%s/CERNHome/D-mesons/Analysis/Models/corr_model2_0.9fm_wC.dat",
          HomeDir.Data()));
  grModel2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grModel2->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grModel2, 20, colModel2, 0.8);

  auto grModel3 = new TGraph(
      TString::Format(
          "%s/CERNHome/D-mesons/Analysis/Models/corr_model3_0.9fm_wC.dat",
          HomeDir.Data()));
  grModel3->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grModel3->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(grModel3, 20, colModel3, 0.8);

  auto grModel4 = new TGraphErrors();
  double x, y, x2, y2;
  auto ingraph1 = new TGraph(
      TString::Format(
          "%s/CERNHome/D-mesons/Analysis/Models/corr_model4_1_0.9fm_wC.dat",
          HomeDir.Data()));
  auto ingraph2 = new TGraph(
      TString::Format(
          "%s/CERNHome/D-mesons/Analysis/Models/corr_model4_2_0.9fm_wC.dat",
          HomeDir.Data()));
  for (int i = 0; i < ingraph1->GetN(); ++i) {
    ingraph1->GetPoint(i, x, y);
    ingraph2->GetPoint(i, x2, y2);
    if (std::abs(x - x2) > 0.01)
      continue;
    grModel4->SetPoint(i, x, 0.5f * (y + y2));
    grModel4->SetPointError(i, 0., 0.5 * (y - y2));
  }
  DreamPlot::SetStyleGraph(grModel4, 20, colModel4, 0.8);
  grModel4->SetFillColorAlpha(colModel4, 0.8);
  auto grModel4Dummy = (TGraphErrors*) grModel4->Clone("grModelDummy");
  grModel4Dummy->SetLineWidth(3);
  grModel4->SetLineWidth(0);

  TidyCats *tidy = new TidyCats();
  CATS chiEFTKitty, CoulombKitty;
  tidy->GetCatsProtonDminus(&chiEFTKitty, momBins, kmin, kmax,
                            TidyCats::pDminusHaidenbauer, TidyCats::sGaussian);
  chiEFTKitty.SetAnaSource(0, radius);
  chiEFTKitty.KillTheCat();
  FillWaveGraph(chiEFTKitty, grHaidenbauer);

  tidy->GetCatsProtonDminus(&CoulombKitty, momBins, kmin, kmax,
                            TidyCats::pDCoulombOnly, TidyCats::sGaussian);
  CoulombKitty.SetAnaSource(0, radius);
  CoulombKitty.KillTheCat();
  FillWaveGraph(CoulombKitty, grCoulomb);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

  const float kmaxdraw = 300;
  const float right = 0.04;
  const float top = 0.025;

  float yminSigma = 0.8;
  const float ymaxSigma = 3.5;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            kmaxdraw);
  DreamPlot::SetStyleHisto(dummyHist, 20, kBlack);

  auto c = new TCanvas("CFpSigma_chiEFT_smear", "CFpSigma_chiEFT_smear", 0, 0,
                       650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grHaidenbauer->Draw("L3same");
  grCoulomb->Draw("L3same");
  grModel1->Draw("L3same");
  grModel3->Draw("l3same");
  grModel4->Draw("l3same");

  float xmin = 0.5;
  float xwidth = 0.45;
  float ymin = 0.65;
  float ywidth = 0.3;
  auto leg = new TLegend(xmin, ymin, xmin + xwidth, ymin + ywidth);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.8);
  leg->AddEntry(grCoulomb, "Coulomb", "l");
  leg->AddEntry(grHaidenbauer, "J. Haidenbauer #it{et al.}", "l");
  leg->AddEntry(grModel1, "J. Hofmann and M. Lutz", "l");
  leg->AddEntry(grModel3, "Y. Yamaguchi #it{et al.}", "l");
  leg->AddEntry(grModel4Dummy, "C. Fontoura #it{et al.}", "l");
  leg->SetHeader(Form("#it{r}_{0} = %.1f fm", radius));
  leg->Draw("same");
  c->Print("models.pdf");

  auto outfile = new TFile("models.root", "RECREATE");
  grHaidenbauer->Write("Haidenbauer");
  grCoulomb->Write("Coulomb");
  grModel1->Write("Model1");
  grModel2->Write("Model2");
  grModel3->Write("Model3");
  grModel4->Write("Model4");
  outfile->Close();
  return 1;
}
