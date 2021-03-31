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
double integration(TGraph* pw, TF1* source, const double stepsize) {
  double x, y;
  double val = 0;
  for (int i=0; i<pw->GetN(); ++i) {
    pw->GetPoint(i, x, y);
    val += std::abs(y) * std::abs(y) * source->Eval(x) * stepsize;
  }
  return 1 + val;
}
/// =====================================================================================
double integration(TGraph* pwReal, TGraph* pwImag, TF1* source, const double stepsize) {
  double xReal, yReal, xImag, yImag;
  double val = 0;
  for (int i=0; i<pwReal->GetN(); ++i) {
    pwReal->GetPoint(i, xReal, yReal);
    pwImag->GetPoint(i, xImag, yImag);
    if ( xReal > 10 ) break;
    std::complex<double> wf(yReal, yImag);
    std::cout << xReal << " " << yReal << " " << yImag << " " << std::abs(wf) << " " << val << " " << std::abs(wf) * std::abs(wf) * source->Eval(xReal) << "\n";
    val += std::abs(wf) * std::abs(wf) * source->Eval(xReal) * stepsize;
  }
  return 1 + val;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  double *radius = new double[1];
  radius[0] = 1.;
  int momBins = 100;
  int kmin = 1;
  int kmax = 201;

  auto gaussFit =
      new TF1(
          "gaus",
          [&](double *x, double *p) {return
	    4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
            std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
          },
          0, 10, 1);
  gaussFit->SetParameter(0, radius[0]);
  gaussFit->SetNpx(1000);

  const int colorchiEFT = kAzure;
  const float alphaEFT = 1;
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
  CATS chiEFTKitty, CoulombKitty;
  tidy->GetCatsProtonDminus(&chiEFTKitty, momBins, kmin, kmax,
                            TidyCats::pDminusHaidenbauer, TidyCats::sGaussian);
  chiEFTKitty.SetAnaSource(0, radius[0]);
  chiEFTKitty.KillTheCat();

  tidy->GetCatsProtonDminus(&CoulombKitty, momBins, kmin, kmax,
                            TidyCats::pDCoulombOnly, TidyCats::sGaussian);
  CoulombKitty.SetAnaSource(0, radius[0]);
  CoulombKitty.KillTheCat();

  auto grCoulomb = (TGraph*)grpDminuschiEFT->Clone("Coulomb");
  FillWaveGraph(CoulombKitty, grCoulomb);
  
  FillWaveGraph(chiEFTKitty, grpDminusElastchiEFT, { 1, 0 });
  FillWaveGraph(chiEFTKitty, grpDminusnDzerochiEFT, { 0, 1 });
  FillWaveGraph(chiEFTKitty, grpDminuschiEFT, {1, 1});

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

  // evaluate the wave functions
  bool divideByR = true;
  auto outputFile = new TFile("wavefunctions.root", "recreate");
  c->Write("cfs");
  grCoulomb->Write("CF_Coulomb");
  grpDminuschiEFT->Write("CF_tot");
  grpDminusElastchiEFT->Write("CF_elast");
  grpDminusnDzerochiEFT->Write("CF_cc");

  auto grIntPsi = new TGraph();
  grIntPsi->SetName("grIntPsi");
  auto grIntChi = new TGraph();
  grIntChi->SetName("grIntChi");
  
  for (unsigned int i = 0; i < chiEFTKitty.GetNumMomBins(); ++i) {
    const double mom = chiEFTKitty.GetMomentum(i);

    auto grWF2 = new TGraph();
    
    auto grCurpDReal = new TGraph();
    grCurpDReal->SetName(Form("pDReal_%.3f", mom));
    auto grCurpDImag = new TGraph();
    grCurpDImag->SetName(Form("pDImag_%.3f", mom));

    auto grCurnDReal = new TGraph();
    grCurnDReal->SetName(Form("nDReal_%.3f", mom));
    auto grCurnDImag = new TGraph();
    grCurnDImag->SetName(Form("nDImag_%.3f", mom));

    auto grRefpDReal = new TGraph();
    grRefpDReal->SetName(Form("pDrefReal_%.3f", mom));
    
    double stepsize = 0.01;
    double rad = stepsize;
    int counter = 0;
    while (rad < 20) {
      auto pDVal = chiEFTKitty.EvalRadialWaveFunction(i, 0, 0, rad, divideByR);
      grCurpDReal->SetPoint(counter, rad, pDVal.real());
      grCurpDImag->SetPoint(counter, rad, pDVal.imag());

      auto nDVal = chiEFTKitty.EvalRadialWaveFunction(i, 1, 0, rad, divideByR);
      grCurnDReal->SetPoint(counter, rad, nDVal.real());
      grCurnDImag->SetPoint(counter, rad, nDVal.imag());

      auto ref = chiEFTKitty.EvalReferenceRadialWF(i, 0, rad, divideByR);
      grRefpDReal->SetPoint(counter, rad, ref.real());

      rad += stepsize;
      counter++;
    }
    grIntPsi->SetPoint(i, mom, integration(grRefpDReal, gaussFit, stepsize));
    grIntChi->SetPoint(i, mom, integration(grCurpDReal, grCurpDImag, gaussFit, stepsize));
    
    grCurpDReal->Write();
    grCurpDImag->Write();
    grCurnDReal->Write();
    grCurnDImag->Write();
    grRefpDReal->Write();

    delete grCurpDReal;
    delete grCurpDImag;
    delete grCurnDReal;
    delete grCurnDImag;
    delete grRefpDReal;
  }
  gaussFit->Write("source");
  grIntPsi->Write();
  grIntChi->Write();
  outputFile->Close();
  return 1;
}
