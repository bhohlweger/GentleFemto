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
  radius[0] = 1.25;
  int momBins = 50;
  int kmin = -4.99;
  int kmax = 495.01;

  const int colorfss2 = kRed + 2;
  const float alphafss2 = 0.5;
  const int colorchiEFT = kAzure;
  const float alphaEFT = 0.5;
  const int colorESC16 = kGreen + 2;
  const float alphaESC16 = 0.6;
  const int colorNSC97f = kOrange - 3;
  const float alphaNSC97f = 0.5;
  const int lineWidth = 3;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Proton-Sigma0 correlation functions
  auto grSigma0fss2 = new TGraph();
  grSigma0fss2->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0fss2, 20, colorfss2, alphafss2);
  grSigma0fss2->SetLineWidth(lineWidth);
  auto grSigma0chiEFT = new TGraph();
  grSigma0chiEFT->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0chiEFT, 20, colorchiEFT, alphaEFT);
  grSigma0chiEFT->SetLineWidth(lineWidth);
  auto grSigma0ESC16 = new TGraph();
  grSigma0ESC16->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0ESC16, 20, colorESC16, alphaESC16);
  grSigma0ESC16->SetLineWidth(lineWidth);
  auto grSigma0NSC97f = new TGraph();
  grSigma0NSC97f->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grSigma0NSC97f, 20, colorNSC97f, alphaNSC97f);
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

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Smearing
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();

  auto grSigma0fss2_smear = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0fss2_smear, 20, colorfss2, alphafss2);
  grSigma0fss2_smear->SetLineWidth(lineWidth);
  grSigma0fss2_smear->SetLineStyle(2);
  auto grSigma0chiEFT_smear = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0chiEFT_smear, 20, colorchiEFT, alphaEFT);
  grSigma0chiEFT_smear->SetLineWidth(lineWidth);
  grSigma0chiEFT_smear->SetLineStyle(2);
  auto grSigma0ESC16_smear = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0ESC16_smear, 20, colorESC16, alphaESC16);
  grSigma0ESC16_smear->SetLineWidth(lineWidth);
  grSigma0ESC16_smear->SetLineStyle(2);
  auto grSigma0NSC97f_smear = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0NSC97f_smear, 20, colorNSC97f, alphaNSC97f);
  grSigma0NSC97f_smear->SetLineWidth(lineWidth);
  grSigma0NSC97f_smear->SetLineStyle(2);

  DLM_CkDecomposition CkDec_Sigma0chiEFT_smear("CkDec_Sigma0chiEFT_smear", 1,
                                               *Sigma0chiEFT,
                                               CATSinput->GetSigmaFile(1));
  DLM_CkDecomposition CkDec_Sigma0ESC16_smear("CkDec_Sigma0ESC16_smear", 1,
                                              *Sigma0ESC16,
                                              CATSinput->GetSigmaFile(1));
  DLM_CkDecomposition CkDec_Sigma0NSC97f_smear("CkDec_Sigma0NSC97f_smear", 1,
                                               *Sigma0NSC97f,
                                               CATSinput->GetSigmaFile(1));
  DLM_CkDecomposition CkDec_Sigma0fss2_smear("CkDec_Sigma0fss2_smear", 1,
                                             *Sigma0fss2,
                                             CATSinput->GetSigmaFile(1));

  for (unsigned int i = 0; i < Sigma0chiEFT->GetNbins(); ++i) {
    const float mom = Sigma0chiEFT->GetBinCenter(0, i);
    grSigma0chiEFT_smear->SetPoint(i, mom,
                                   CkDec_Sigma0chiEFT_smear.EvalCk(mom));
    grSigma0ESC16_smear->SetPoint(i, mom, CkDec_Sigma0ESC16_smear.EvalCk(mom));
    grSigma0NSC97f_smear->SetPoint(i, mom,
                                   CkDec_Sigma0NSC97f_smear.EvalCk(mom));
    grSigma0fss2_smear->SetPoint(i, mom, CkDec_Sigma0fss2_smear.EvalCk(mom));
  }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Lambda Param

  const double lambdaPrimary = 0.22;

  auto grSigma0fss2_lambda = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0fss2_lambda, 20, colorfss2, alphafss2);
  grSigma0fss2_lambda->SetLineWidth(lineWidth);
  grSigma0fss2_lambda->SetLineStyle(2);
  auto grSigma0chiEFT_lambda = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0chiEFT_lambda, 20, colorchiEFT, alphaEFT);
  grSigma0chiEFT_lambda->SetLineWidth(lineWidth);
  grSigma0chiEFT_lambda->SetLineStyle(2);
  auto grSigma0ESC16_lambda = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0ESC16_lambda, 20, colorESC16, alphaESC16);
  grSigma0ESC16_lambda->SetLineWidth(lineWidth);
  grSigma0ESC16_lambda->SetLineStyle(2);
  auto grSigma0NSC97f_lambda = new TGraph();
  DreamPlot::SetStyleGraph(grSigma0NSC97f_lambda, 20, colorNSC97f, alphaNSC97f);
  grSigma0NSC97f_lambda->SetLineWidth(lineWidth);
  grSigma0NSC97f_lambda->SetLineStyle(2);


  DLM_CkDecomposition CkDec_Sigma0chiEFT_lambda("CkDec_Sigma0chiEFT_lambda", 1,
                                                *Sigma0chiEFT, nullptr);
  DLM_CkDecomposition CkDec_Sigma0ESC16_lambda("CkDec_Sigma0ESC16_lambda", 1,
                                               *Sigma0ESC16, nullptr);
  DLM_CkDecomposition CkDec_Sigma0NSC97f_lambda("CkDec_Sigma0NSC97f_lambda", 1,
                                                *Sigma0NSC97f, nullptr);
  DLM_CkDecomposition CkDec_Sigma0fss2_lambda("CkDec_Sigma0fss2_lambda", 1,
                                              *Sigma0fss2, nullptr);

  CkDec_Sigma0chiEFT_lambda.AddContribution(0, 1.f - lambdaPrimary, DLM_CkDecomposition::cFake);
  CkDec_Sigma0chiEFT_lambda.Update();
  CkDec_Sigma0ESC16_lambda.AddContribution(0, 1.f - lambdaPrimary, DLM_CkDecomposition::cFake);
  CkDec_Sigma0ESC16_lambda.Update();
  CkDec_Sigma0NSC97f_lambda.AddContribution(0, 1.f - lambdaPrimary, DLM_CkDecomposition::cFake);
  CkDec_Sigma0NSC97f_lambda.Update();
  CkDec_Sigma0fss2_lambda.AddContribution(0, 1.f - lambdaPrimary, DLM_CkDecomposition::cFake);
  CkDec_Sigma0fss2_lambda.Update();

  for (unsigned int i = 0; i < Sigma0chiEFT->GetNbins(); ++i) {
    const float mom = Sigma0chiEFT->GetBinCenter(0, i);
    grSigma0chiEFT_lambda->SetPoint(i, mom,
                                   CkDec_Sigma0chiEFT_lambda.EvalCk(mom));
    grSigma0ESC16_lambda->SetPoint(i, mom, CkDec_Sigma0ESC16_lambda.EvalCk(mom));
    grSigma0NSC97f_lambda->SetPoint(i, mom,
                                   CkDec_Sigma0NSC97f_lambda.EvalCk(mom));
    grSigma0fss2_lambda->SetPoint(i, mom, CkDec_Sigma0fss2_lambda.EvalCk(mom));
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

  const float kmaxdraw = 200;
  const float right = 0.04;
  const float top = 0.025;

  float yminSigma = 0.25;
  const float ymaxSigma = 2.5;

  auto dummyHist = new TH1F("dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, kmaxdraw);
  DreamPlot::SetStyleHisto(dummyHist, 20, colorNSC97f);

  auto dummyGraph = new TGraph();
  DreamPlot::SetStyleGraph(dummyGraph, 20, kBlack);
  auto dummyGraphSmeared = new TGraph();
  dummyGraph->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(dummyGraphSmeared, 20, kBlack);
  dummyGraphSmeared->SetLineWidth(lineWidth);
  dummyGraphSmeared->SetLineStyle(2);

  auto c = new TCanvas("CFpSigma_smear", "CFpSigma_smear", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grSigma0chiEFT->Draw("L3same");
  grSigma0ESC16->Draw("L3same");
  grSigma0NSC97f->Draw("L3same");
  grSigma0fss2->Draw("L3same");

  grSigma0chiEFT_smear->Draw("L3same");
  grSigma0ESC16_smear->Draw("L3same");
  grSigma0NSC97f_smear->Draw("L3same");
  grSigma0fss2_smear->Draw("L3same");

  auto leg = new TLegend(0.49, 0.725, 0.49+0.45, 0.725+0.225);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetNColumns(2);
  leg->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg->AddEntry(grSigma0fss2, "fss2", "l");
  leg->AddEntry(grSigma0chiEFT, "#chiEFT (NLO)", "l");
  leg->AddEntry(grSigma0ESC16, "ESC16", "l");
  leg->AddEntry(grSigma0NSC97f, "NSC97f", "l");
  leg->AddEntry(dummyGraph, "Genuine", "l");
  leg->AddEntry(dummyGraphSmeared, "Smeared", "l");
  leg->Draw("same");
  leg->Draw("same");
  c->Print("Sigma0_smeared.pdf");

  auto d = new TCanvas("CFpSigma_lambda", "CFpSigma_lambda", 0, 0, 650, 550);
  d->SetRightMargin(right);
  d->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminSigma, ymaxSigma);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grSigma0chiEFT->Draw("L3same");
  grSigma0ESC16->Draw("L3same");
  grSigma0NSC97f->Draw("L3same");
  grSigma0fss2->Draw("L3same");

  grSigma0chiEFT_lambda->Draw("L3same");
  grSigma0ESC16_lambda->Draw("L3same");
  grSigma0NSC97f_lambda->Draw("L3same");
  grSigma0fss2_lambda->Draw("L3same");

  auto leg2 = new TLegend(0.49, 0.725, 0.49+0.45, 0.725+0.225);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetNColumns(2);
  leg2->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg2->AddEntry(grSigma0fss2, "fss2", "l");
  leg2->AddEntry(grSigma0chiEFT, "#chiEFT (NLO)", "l");
  leg2->AddEntry(grSigma0ESC16, "ESC16", "l");
  leg2->AddEntry(grSigma0NSC97f, "NSC97f", "l");
  leg2->AddEntry(dummyGraph, "Genuine", "l");
  leg2->AddEntry(dummyGraphSmeared, "#lambda = 0.22", "l");
  leg2->Draw("same");
  leg2->Draw("same");
  d->Print("Sigma0_lambda.pdf");

}
