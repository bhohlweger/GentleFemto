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
#include "TLatex.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}

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
  double* radius = new double[1];
  radius[0] = 1.154;
  int momBins = 80;
  int kmin = -2.49;
  int kmax = 397.51;

  const int colorfss2 = kRed + 1;
  const int colorchiEFT = kAzure;
  const int colorESC16 = kGreen + 2;
  const int colorNSC97f = kOrange - 3;
  const int lineWidth = 3;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Proton-Lambda correlation functions
  auto grLambdafss2 = new TGraph();
  grLambdafss2->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grLambdafss2, 20, colorfss2);
  grLambdafss2->SetLineWidth(lineWidth);
  auto grLambdachiEFT = new TGraph();
  grLambdachiEFT->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grLambdachiEFT, 20, colorchiEFT);
  grLambdachiEFT->SetLineWidth(lineWidth);
  auto grLambdaESC16 = new TGraph();
  grLambdaESC16->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grLambdaESC16, 20, colorESC16);
  grLambdaESC16->SetLineWidth(lineWidth);
  auto grLambdaNSC97f = new TGraph();
  grLambdaNSC97f->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grLambdaNSC97f, 20, colorNSC97f);
  grLambdaNSC97f->SetLineWidth(lineWidth);

  /// chiEFT
  auto LambdaChiEFTHaidenbauer = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                            Lednicky_SingletTriplet);
  LambdaChiEFTHaidenbauer->SetPotPar(0, 2.91);
  LambdaChiEFTHaidenbauer->SetPotPar(1, 2.78);
  LambdaChiEFTHaidenbauer->SetPotPar(2, 1.54);
  LambdaChiEFTHaidenbauer->SetPotPar(3, 2.72);
  LambdaChiEFTHaidenbauer->SetSourcePar(0, radius[0]);
  LambdaChiEFTHaidenbauer->Update();
  FillCkGraph(LambdaChiEFTHaidenbauer, grLambdachiEFT);

  /// fss2
  auto Lambdafss2 = new DLM_Ck(1, 4, momBins, kmin, kmax,
                               Lednicky_SingletTriplet);
  Lambdafss2->SetPotPar(0, 2.59);
  Lambdafss2->SetPotPar(1, 2.83);
  Lambdafss2->SetPotPar(2, 1.60);
  Lambdafss2->SetPotPar(3, 3.01);
  Lambdafss2->SetSourcePar(0, radius[0]);
  Lambdafss2->Update();
  FillCkGraph(Lambdafss2, grLambdafss2);

  /// ESC16
  auto LambdaESC16 = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                Lednicky_SingletTriplet);
  LambdaESC16->SetPotPar(0, 1.88);
  LambdaESC16->SetPotPar(1, 3.58);
  LambdaESC16->SetPotPar(2, 1.86);
  LambdaESC16->SetPotPar(3, 3.37);
  LambdaESC16->SetSourcePar(0, radius[0]);
  LambdaESC16->Update();
  FillCkGraph(LambdaESC16, grLambdaESC16);

  /// NSC97f
  auto LambdaNSC97f = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                 Lednicky_SingletTriplet);
  LambdaNSC97f->SetPotPar(0, 2.51);
  LambdaNSC97f->SetPotPar(1, 3.03);
  LambdaNSC97f->SetPotPar(2, 1.75);
  LambdaNSC97f->SetPotPar(3, 3.32);
  LambdaNSC97f->SetSourcePar(0, radius[0]);
  LambdaNSC97f->Update();
  FillCkGraph(LambdaNSC97f, grLambdaNSC97f);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Proton-Sigma0 correlation functions
  auto grSigma0fss2 = new TGraph();
  grSigma0fss2->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0fss2, 20, colorfss2);
  grSigma0fss2->SetLineWidth(lineWidth);
  auto grSigma0chiEFT = new TGraph();
  grSigma0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0chiEFT, 20, colorchiEFT);
  grSigma0chiEFT->SetLineWidth(lineWidth);
  auto grSigma0ESC16 = new TGraph();
  grSigma0ESC16->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0ESC16, 20, colorESC16);
  grSigma0ESC16->SetLineWidth(lineWidth);
  auto grSigma0NSC97f = new TGraph();
  grSigma0NSC97f->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0NSC97f, 20, colorNSC97f);
  grSigma0NSC97f->SetLineWidth(lineWidth);

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

  /// ESC16
  CATS ESC16Kitty;
  tidy->GetCatsProtonSigma0(&ESC16Kitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0ESC16);
  ESC16Kitty.KillTheCat();
  auto Sigma0ESC16 = new DLM_Ck(1, 0, ESC16Kitty);
  Sigma0ESC16->SetSourcePar(0, radius[0]);
  Sigma0ESC16->Update();
  FillWaveGraph(ESC16Kitty, grSigma0ESC16);

  /// NSC97f
  CATS NSC97fKitty;
  tidy->GetCatsProtonSigma0(&NSC97fKitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
  NSC97fKitty.KillTheCat();
  auto Sigma0NSC97f = new DLM_Ck(1, 0, NSC97fKitty);
  Sigma0NSC97f->SetSourcePar(0, radius[0]);
  Sigma0NSC97f->Update();
  FillWaveGraph(NSC97fKitty, grSigma0NSC97f);

  /// fss2
  DLM_Ck* Sigma0fss2 = new DLM_Ck(1, 0, momBins, kmin, kmax,
                                  Lednicky_gauss_Sigma0);
  Sigma0fss2->SetSourcePar(0, radius[0]);
  Sigma0fss2->Update();
  FillCkGraph(Sigma0fss2, grSigma0fss2);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting
  gStyle->SetHatchesSpacing(0.75);
  const float kmaxdraw = 300;
  const float right = 0.03;
  const float top = 0.025;


  auto d = new TCanvas("CFpLambda", "CFpLambda", 0, 0, 650, 550);
  d->SetRightMargin(right);
  d->SetTopMargin(top);
  grLambdachiEFT->Draw("AL3");
  grLambdachiEFT->GetXaxis()->SetRangeUser(0, kmaxdraw);
  grLambdachiEFT->GetYaxis()->SetRangeUser(0.8, 3.5);
  grLambdaESC16->Draw("L3same");
  grLambdaNSC97f->Draw("L3same");
  grLambdafss2->Draw("L3same");
  float xmin = 0.49;
  float xwidth = 0.45;
  float ymin = 0.725;
  float ywidth = 0.225;
  auto leg = new TLegend(xmin, ymin, xmin+xwidth, ymin+ywidth);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetNColumns(2);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg->AddEntry(grSigma0fss2, "fss2", "l");
  leg->AddEntry(grSigma0chiEFT, "#chiEFT (NLO)", "l");
  leg->AddEntry(grSigma0ESC16, "ESC16", "l");
  leg->AddEntry(grSigma0NSC97f, "NSC97f", "l");
  leg->Draw("same");
  d->Print("Lambda.pdf");

  auto c = new TCanvas("CFpSigma", "CFpSigma", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  grSigma0chiEFT->Draw("AL3");
  grSigma0chiEFT->GetXaxis()->SetRangeUser(0, kmaxdraw);
  grSigma0chiEFT->GetYaxis()->SetRangeUser(0.75, 2.75);
  grSigma0ESC16->Draw("L3same");
  grSigma0NSC97f->Draw("L3same");
  grSigma0fss2->Draw("L3same");
  leg->Draw("same");
  c->Print("Sigma0.pdf");
  xmin = 0.33;
  xwidth = 0.6;
  ymin = 0.31;
  ywidth = 0.39;
  TPad *inset_pad = new TPad("insert", "insertPad", xmin, ymin, xmin+xwidth, ymin+ywidth);
  inset_pad->SetTopMargin(0.01);
  inset_pad->SetRightMargin(0.05);
  inset_pad->SetBottomMargin(0.28);
  inset_pad->SetLeftMargin(0.28);
  inset_pad->SetFillStyle(4000);
  inset_pad->Draw();
  inset_pad->cd();
  grLambdachiEFT->Draw("AL3");
  grLambdachiEFT->GetXaxis()->SetRangeUser(0, kmaxdraw);
  grLambdachiEFT->GetYaxis()->SetRangeUser(0.8, 3.5);
  grLambdachiEFT->GetYaxis()->SetNdivisions(203);
  grLambdachiEFT->GetXaxis()->SetNdivisions(204);
  grLambdachiEFT->GetXaxis()->SetTitleOffset(2.5);
  grLambdachiEFT->GetYaxis()->SetTitleOffset(1.1);
  const float scaleSize = 0.85;
  grLambdachiEFT->GetXaxis()->SetTitleSize(scaleSize * grLambdachiEFT->GetXaxis()->GetTitleSize());
  grLambdachiEFT->GetYaxis()->SetTitleSize(scaleSize * grLambdachiEFT->GetYaxis()->GetTitleSize());
  grLambdachiEFT->GetXaxis()->SetLabelSize(scaleSize * grLambdachiEFT->GetXaxis()->GetLabelSize());
  grLambdachiEFT->GetYaxis()->SetLabelSize(scaleSize * grLambdachiEFT->GetYaxis()->GetLabelSize());
  grLambdaESC16->Draw("L3same");
  grLambdaNSC97f->Draw("L3same");
  grLambdafss2->Draw("L3same");
  c->cd();
  TLatex pair;
  pair.SetTextSize(gStyle->GetTextSize() * 1.25);
  pair.SetNDC(kTRUE);
  pair.DrawLatex(0.275, 0.885, "p#minus#kern[-0.65]{ }#Sigma^{0}");
  pair.DrawLatex(0.775, 0.625, "p#minus#kern[-0.65]{ }#Lambda");
  c->Print("Sigma0-Lambda.pdf");
}
