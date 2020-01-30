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


  // Repeat the exercise for pp
  auto grPPFull = new TGraph();
  grPPFull->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grPPFull, 20, kTeal + 2, 0.8);
  grPPFull->SetLineWidth(lineWidth);
  auto grPPSmear = new TGraph();
  DreamPlot::SetStyleGraph(grPPSmear, 20, kTeal + 2, 0.8);
  grPPSmear->SetLineWidth(lineWidth);
  grPPSmear->SetLineStyle(2);
  auto grPPLambda = new TGraph();
  DreamPlot::SetStyleGraph(grPPLambda, 20,  kTeal + 2, 0.8);
  grPPLambda->SetLineWidth(lineWidth);
  grPPLambda->SetLineStyle(2);

  CATS ppKitty;
  tidy->GetCatsProtonProton(&ppKitty, 150, 0, 300, TidyCats::sGaussian);
  ppKitty.KillTheCat();
  auto ppCk = new DLM_Ck(1, 0, ppKitty);
  ppCk->SetSourcePar(0, radius[0]);
  ppCk->Update();
  FillWaveGraph(ppKitty, grPPFull);

  const float lambdaProtonPrim = 0.67;
  DLM_CkDecomposition CkDec_pp_smear("CkDec_pp_smear", 1, *ppCk,
                                     CATSinput->GetSigmaFile(0));
  DLM_CkDecomposition CkDec_pp_lambda("CkDec_pp_lambda", 1, *ppCk, nullptr);
  CkDec_pp_lambda.AddContribution(0, 1.f - lambdaProtonPrim, DLM_CkDecomposition::cFake);
  CkDec_pp_lambda.Update();

  for (unsigned int i = 0; i < ppCk->GetNbins(); ++i) {
    const float mom = ppCk->GetBinCenter(0, i);
    grPPSmear->SetPoint(i, mom, CkDec_pp_smear.EvalCk(mom));
    grPPLambda->SetPoint(i, mom, CkDec_pp_lambda.EvalCk(mom));
  }

  float yminpp = 0.;
  const float ymaxpp = 5;

  auto c2 = new TCanvas("CFpp_smear", "CFpp_smear", 0, 0, 650, 550);
  c2->SetRightMargin(right);
  c2->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminpp, ymaxpp);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grPPFull->Draw("L3same");
  grPPSmear->Draw("L3same");

  auto leg4 = new TLegend(0.49, 0.725, 0.49+0.45, 0.725+0.225);
  leg4->SetBorderSize(0);
  leg4->SetTextFont(42);
  leg4->SetHeader("p#minus#kern[-0.95]{ }p Coulomb + Argonne #nu_{18}");
  leg4->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg4->AddEntry(grPPFull, "Genuine", "l");
  leg4->AddEntry(grPPSmear, "Smeared", "l");
  leg4->Draw("same");
  c2->Print("pp_smeared.pdf");

  auto c3 = new TCanvas("CFpp_lambda", "CFpp_lambda", 0, 0, 650, 550);
  c3->SetRightMargin(right);
  c3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(yminpp, ymaxpp);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grPPFull->Draw("L3same");
  grPPLambda->Draw("L3same");

  auto leg5 = new TLegend(0.49, 0.725, 0.49+0.45, 0.725+0.225);
  leg5->SetBorderSize(0);
  leg5->SetTextFont(42);
  leg5->SetHeader("p#minus#kern[-0.95]{ }p Coulomb + Argonne #nu_{18}");
  leg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg5->AddEntry(grPPFull, "Genuine", "l");
  leg5->AddEntry(grPPLambda, "#lambda = 0.67", "l");
  leg5->Draw("same");
  c3->Print("pp_lambda.pdf");

  // Now the Lambda feeding to pp
  auto grPLFull = new TGraph();
  grPLFull->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grPLFull, 20, kBlue + 3, 0.7);
  grPLFull->SetLineWidth(lineWidth);
  auto grPLSmear = new TGraph();
  DreamPlot::SetStyleGraph(grPLSmear, 20, kBlue + 3, 0.7);
  grPLSmear->SetLineWidth(lineWidth);
  grPLSmear->SetLineStyle(2);
  auto grPLLambda = new TGraph();
  DreamPlot::SetStyleGraph(grPLLambda, 20,  kBlue + 3, 0.7);
  grPLLambda->SetLineWidth(lineWidth);
  grPLLambda->SetLineStyle(3);

  const float PurityLambda = 0.961;
  const float PrimLambdaAndSigma = 0.785;  //fraction of primary Lambdas + Sigma 0
  const float SecLambda = 1 - PrimLambdaAndSigma;  //fraction of weak decay Lambdas


  const float lambdaLambdaProton = 0.203;
  CATS AB_pL;
  tidy->GetCatsProtonLambda(&AB_pL, 31, -9.99, 300, TidyCats::sGaussian,
                            TidyCats::pNLOWF);
  AB_pL.KillTheCat();
  DLM_Ck* Ck_pL = new DLM_Ck(1, 0, AB_pL);
  Ck_pL->SetSourcePar(0, radius[0]);
  Ck_pL->Update();
  FillWaveGraph(AB_pL, grPLFull);

  DLM_CkDecomposition CkDec_pL("pLambda", 0, *Ck_pL,
                               CATSinput->GetSigmaFile(1));

  DLM_Ck* Ck_pL_flat = new DLM_Ck(1, 0, 150, 0, 300, Flat_Residual);
  Ck_pL_flat->SetSourcePar(0, radius[0]);
  Ck_pL_flat->Update();
  DLM_CkDecomposition CkDec_pL_lambdasmear("pLambda2", 1, *Ck_pL_flat, nullptr);
  CkDec_pL_lambdasmear.AddContribution(0, 0.999, DLM_CkDecomposition::cFeedDown,
                                       &CkDec_pL, CATSinput->GetResFile(0));

  DLM_CkDecomposition CkDec_pL_lambdasmearlambda("pLambda3", 1, *Ck_pL_flat, nullptr);
  CkDec_pL_lambdasmearlambda.AddContribution(0, lambdaLambdaProton,
                                             DLM_CkDecomposition::cFeedDown,
                                             &CkDec_pL,
                                             CATSinput->GetResFile(0));
  CkDec_pL_lambdasmear.Update();
  CkDec_pL_lambdasmearlambda.Update();

  for (unsigned int i = 0; i < Ck_pL->GetNbins(); ++i) {
    const float mom = Ck_pL->GetBinCenter(0, i);
    grPLSmear->SetPoint(i, mom, CkDec_pL_lambdasmear.EvalCk(mom));
    grPLLambda->SetPoint(i, mom, CkDec_pL_lambdasmearlambda.EvalCk(mom));
  }

  auto gc3 = new TCanvas("CFpL_lambda", "CFpL_lambda", 0, 0, 650, 550);
  gc3->SetRightMargin(right);
  gc3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 3.5);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grPLFull->Draw("L3same");
  grPLSmear->Draw("L3same");
  grPLLambda->Draw("L3same");

  auto gleg5 = new TLegend(0.45, 0.65, 0.45 + 0.45, 0.95);
  gleg5->SetBorderSize(0);
  gleg5->SetTextFont(42);
  gleg5->SetHeader("p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)");
  gleg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gleg5->AddEntry(grPLFull, "Genuine", "l");
  gleg5->AddEntry(
      grPLSmear,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }p", "l");
  gleg5->AddEntry(
      grPLLambda,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }p (#lambda = 0.203)",
      "l");
  gleg5->Draw("same");
  gc3->Print("pL_feed.pdf");
}
