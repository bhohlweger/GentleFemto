#include "TidyCats.h"
#include "TRandom3.h"
#include "SidebandSigma.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "TDatabasePDG.h"
#include <iostream>
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TStopwatch.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TStyle.h"

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
TGraph *GetSmearedCF(TGraph* CF, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF);
  }
  return smearedCF;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  TString trigger = argv[2];
  TString suffix = argv[3];

  const float lambdaSideband = 0.731;

  auto filename = TString::Format("%s/SherlockSideband.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");

  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  auto side = new SidebandSigma();
  side->SetRebin(10);
  side->SetSideBandFile(InputDir.Data(), trigger.Data(), suffix.Data());
  const double sidebandNormDown = 250;
  const double sidebandNormUp = 400;
  side->SetNormalizationRange(sidebandNormDown, sidebandNormUp);

  side->SideBandCFs();
  auto SBmerge = side->GetSideBandGraph(5);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Get the smearing matrix
  auto filenameSmear = TString::Format("%s/SmearSideband.root", InputDir.Data());
  auto infile = TFile::Open(filenameSmear);
  if (!infile) {
    std::cout
        << "No smearing matrix found - start the sideband computation task!\n";
    return 0;
  }
  auto histSmear = (TH2D*) infile->Get("histSmear");
  auto histSmearSideband = (TH2D*) infile->Get("histSmearSideband");
  infile->Close();

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                       suffix.Data());
  CATSinput->ObtainCFs(10, 250, 400);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  // pp radius
  const double ppRadius = 1.154;

  CATS AB_pXim;
  tidy->GetCatsProtonXiMinus(&AB_pXim, 41, -4.9, 400, TidyCats::sGaussian,
                             TidyCats::pHALQCD, 12);
  AB_pXim.SetAnaSource(0, ppRadius);
  AB_pXim.KillTheCat();
  DLM_Ck* Ck_pXim = new DLM_Ck(1, 0, AB_pXim);
  Ck_pXim->SetSourcePar(0, ppRadius);
  Ck_pXim->Update();

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, 41, -4.9, 400,
                                 TidyCats::sGaussian);
  AB_pXim1530.SetAnaSource(0, ppRadius);
  AB_pXim1530.KillTheCat();
  DLM_Ck* Ck_pXim1530 = new DLM_Ck(1, 0, AB_pXim1530);
  Ck_pXim1530->SetSourcePar(0, ppRadius);
  Ck_pXim1530->Update();

  CATS AB_pSigma0;
  tidy->GetCatsProtonSigma0(&AB_pSigma0, 21, -9.99, 400, TidyCats::sGaussian,
                            TidyCats::pSigma0ESC16);
  AB_pSigma0.KillTheCat();
  DLM_Ck* Ck_pSigma0;
  Ck_pSigma0 = new DLM_Ck(1, 0, AB_pSigma0);
  Ck_pSigma0->SetSourcePar(0, ppRadius);
  Ck_pSigma0->Update();


  CATS AB_pL;
  tidy->GetCatsProtonLambda(&AB_pL, 41, -9.99, 400, TidyCats::sGaussian,
                            TidyCats::pNLOWF);
  AB_pL.KillTheCat();
  DLM_Ck* Ck_pL = new DLM_Ck(1, 0, AB_pL);
  Ck_pL->SetSourcePar(0, ppRadius);
  Ck_pL->Update();

  auto grLambdaGenuine = new TGraph();
  FillWaveGraph(AB_pL, grLambdaGenuine);

  const double protonPurity = 0.9943;
  const double protonPrimary = 0.823;
  const double protonLambda = 0.125;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  // some parameters, stolen from PLOnly.C
  const double lambdaPurity = 0.946;
  const double lambdaPrimaryAndSigma0 = 0.76; // includes primary Lambda and feed-down from Sigma0
  double PrimLambda = 3. / 4. * lambdaPrimaryAndSigma0; // isospin 3 : 1
  double SecSigLambda = 1. / 4. * lambdaPrimaryAndSigma0;  // decay probability = 100%!
  double SecXiLambda = (1. - lambdaPrimaryAndSigma0) / 2.; // same for Xi0 and Xim
  const Particle lambda(lambdaPurity, PrimLambda, { SecSigLambda, SecXiLambda,
                            SecXiLambda });

  // some parameters, stolen from PLOnly.C
  const double xiPurity = 0.915;
  const double Xi01530XimProdFraction = 1 / 2.;
  const double Xim1530XimProdFraction = 1 / 2.;
  const double Xi01530Xim_BR = 2 / 3.;
  const double Xim1530Xim_BR = 1 / 3.;
  const double OmegamXimProdFraction = 1 / 10.;
  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %
  double XiNormalization = 1 + OmegamXimProdFraction * OmegamXim_BR
      + Xi01530XimProdFraction * Xi01530Xim_BR
      + Xim1530XimProdFraction * Xim1530Xim_BR;
  double SecOmegaXim = OmegamXimProdFraction * OmegamXim_BR
      / (double) XiNormalization;
  double SecXi01530Xim = Xi01530XimProdFraction * Xi01530Xim_BR
      / (double) XiNormalization;
  double SecXim1530Xim = Xim1530XimProdFraction * Xim1530Xim_BR
      / (double) XiNormalization;
  double PrimXim = 1. / (double) XiNormalization;
  const Particle xi(xiPurity, PrimXim, { SecOmegaXim, SecXi01530Xim,
                        SecXim1530Xim });


  const CATSLambdaParam lambdaParamPL(proton, lambda);
  const CATSLambdaParam lambdaParamPXi(proton, xi);

  const double lam_pL = lambdaParamPL.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::Primary);
  const double lam_pL_fake = lambdaParamPL.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::Fake)
      + lambdaParamPL.GetLambdaParam(CATSLambdaParam::Fake,
                                   CATSLambdaParam::Primary)
      + lambdaParamPL.GetLambdaParam(CATSLambdaParam::Fake,
                                   CATSLambdaParam::Fake);
  const double lam_pL_pS0 = lambdaParamPL.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0);
  const double lam_pL_pXm = lambdaParamPL.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1);

  const double lam_pXim = lambdaParamPXi.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::Primary);
  const double lam_pXim_pXim1530 = lambdaParamPXi.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
  const double lam_pXim_fake = lambdaParamPXi.GetLambdaParam(
      CATSLambdaParam::Primary, CATSLambdaParam::Fake)
      + lambdaParamPXi.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Primary)
      + lambdaParamPXi.GetLambdaParam(CATSLambdaParam::Fake,
                                      CATSLambdaParam::Fake);

  DLM_CkDecomposition CkDec_pL("pLambda", 4, *Ck_pL, nullptr);  // we apply the momentum smearing later on, doesn't make a difference
  DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0, nullptr);
  DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim, nullptr);
  DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530, nullptr);

  CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
                             DLM_CkDecomposition::cFeedDown, &CkDec_pXim1530,
                             CATSinput->GetResFile(3));  //from Xi-(1530)
  CkDec_pXim.AddContribution(1,
                             1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
                             DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
  CkDec_pXim.AddContribution(2, lam_pXim_fake, DLM_CkDecomposition::cFake);
  CkDec_pXim.Update();

  CkDec_pL.AddContribution(0, lam_pL_pS0, DLM_CkDecomposition::cFeedDown,
                           &CkDec_pSigma0, CATSinput->GetResFile(1));
  CkDec_pL.AddContribution(1, lam_pL_pXm, DLM_CkDecomposition::cFeedDown,
                           &CkDec_pXim, CATSinput->GetResFile(2));
  CkDec_pL.AddContribution(2,
                           1. - lam_pL - lam_pL_pS0 - lam_pL_pXm - lam_pL_fake,
                           DLM_CkDecomposition::cFeedDown);
  CkDec_pL.AddContribution(3, lam_pL_fake, DLM_CkDecomposition::cFake);
  CkDec_pL.Update();

  auto grLambdaDecomposition = new TGraph();
  for (int i = 0; i < Ck_pL->GetNbins(); ++i) {
    const double mom = Ck_pL->GetBinCenter(0, i);
    grLambdaDecomposition->SetPoint(i, mom, CkDec_pL.EvalCk(mom));
  }

  const float lambdaplambda = lambdaParamPL.GetLambdaParam(
      CATSLambdaParam::Primary);

  std::cout << "considering a primary fraction for the p-L CF of "
            << lambdaplambda << "\n";

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // apply the smearing with the photon, and the momentum smearing

  auto grGenuineSmearedSigma = GetSmearedCF(GetSmearedCF(grLambdaGenuine, histSmear),
                                 CATSinput->GetSigmaFile(1));
  auto grDecompSmearedSigma = GetSmearedCF(GetSmearedCF(grLambdaDecomposition, histSmear),
                                       CATSinput->GetSigmaFile(1));

  auto grGenuineSmearedSideband = GetSmearedCF(
      GetSmearedCF(grLambdaGenuine, histSmearSideband),
      CATSinput->GetSigmaFile(1));
  auto grDecompSmearedSideband = GetSmearedCF(
      GetSmearedCF(grLambdaDecomposition, histSmearSideband),
      CATSinput->GetSigmaFile(1));

  auto grSmearedSigmaForSigma = new TGraph();
  auto grSmearedSidebandForSigma= new TGraph();

  double x, y;
  for (int i = 0; i < 1000; ++i) {
    grDecompSmearedSigma->GetPoint(i, x, y);
    grSmearedSigmaForSigma->SetPoint(i, x, 1 + (y - 1) * lambdaSideband);
  }

  for (int i = 0; i < 1000; ++i) {
    grDecompSmearedSideband->GetPoint(i, x, y);
    grSmearedSidebandForSigma->SetPoint(i, x, 1 + (y - 1) * lambdaSideband);
  }

  auto grDeviationSigma = new TGraph();
  double x1, x2, y1, y2;
  for (int i = 0; i < 1000; ++i) {
    grSmearedSigmaForSigma->GetPoint(i, x1, y1);
    grSmearedSidebandForSigma->GetPoint(i, x2, y2);
    grDeviationSigma->SetPoint(i, x1, (y1-y2)/y1);
  }


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Plotting

  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            350);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  DreamPlot::SetStyleHisto(histSmear);
  auto c = new TCanvas("c", "c", 650, 550);
  histSmear->Draw("col");
  histSmear->GetXaxis()->SetNdivisions(505);
  histSmear->GetYaxis()->SetNdivisions(505);
  histSmear->SetTitle(
      "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})");
  c->Print("pLambdaSmearingMatrixSigma.pdf");

  DreamPlot::SetStyleHisto(histSmearSideband);
  auto d = new TCanvas("d", "d", 650, 550);
  histSmearSideband->Draw("col");
  histSmearSideband->GetXaxis()->SetNdivisions(505);
  histSmearSideband->GetYaxis()->SetNdivisions(505);
  histSmearSideband->SetTitle(
      "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus(#Lambda#gamma)} (MeV/#it{c})");
  d->Print("pLambdaSmearingMatrixSideband.pdf");

  int lineWidth = 2;
  DreamPlot::SetStyleGraph(grLambdaGenuine, 20, kGreen + 3, 0.7);
  grLambdaGenuine->SetLineWidth(lineWidth);
  grLambdaGenuine->SetLineStyle(3);
  DreamPlot::SetStyleGraph(grGenuineSmearedSigma, 20, kGreen + 3, 0.7);
  grGenuineSmearedSigma->SetLineWidth(lineWidth);
  grGenuineSmearedSigma->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grDecompSmearedSigma, 20, kGreen + 3, 0.7);
  grDecompSmearedSigma->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(SBmerge, kOpenCircle, kBlue + 3);

  auto gc3 = new TCanvas("CFpL_lambda", "CFpL_lambda", 0, 0, 650, 550);
  gc3->SetRightMargin(right);
  gc3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 3.5);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grLambdaGenuine->Draw("L3same");
  grGenuineSmearedSigma->Draw("L3same");
  grDecompSmearedSigma->Draw("L3same");
  SBmerge->Draw("pezsame");
  auto gleg5 = new TLegend(0.4, 0.6, 0.4 + 0.45, 0.9);
  gleg5->SetBorderSize(0);
  gleg5->SetTextFont(42);
  gleg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gleg5->AddEntry(grLambdaGenuine,
                  "Genuine p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)", "l");
  gleg5->AddEntry(
      grGenuineSmearedSigma,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}",
      "l");
  gleg5->AddEntry(
      grDecompSmearedSigma,
      Form(
          "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0} (#lambda = %.3f)",
          lambdaplambda),
      "l");
  gleg5->AddEntry(
      SBmerge,
      "p#minus#kern[-0.65]{ }(#Lambda#gamma) #oplus #bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma)",
      "pez");
  gleg5->Draw("same");
  gc3->Print("pLambdaSmearingSigma.pdf");

  DreamPlot::SetStyleGraph(grGenuineSmearedSideband, 20, kGreen + 3, 0.7);
  grGenuineSmearedSideband->SetLineWidth(lineWidth);
  grGenuineSmearedSideband->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grDecompSmearedSideband, 20, kGreen + 3, 0.7);
  grDecompSmearedSideband->SetLineWidth(lineWidth);

  auto g2c3 = new TCanvas("CFpL_lambda_SB", "CFpL_lambda_SB", 0, 0, 650, 550);
  g2c3->SetRightMargin(right);
  g2c3->SetTopMargin(top);
  dummyHist->Draw();
  grLambdaGenuine->Draw("L3same");
  grGenuineSmearedSideband->Draw("L3same");
  grDecompSmearedSideband->Draw("L3same");
  SBmerge->Draw("pezsame");
  auto gle1g5 = new TLegend(0.4, 0.6, 0.4 + 0.45, 0.9);
  gle1g5->SetBorderSize(0);
  gle1g5->SetTextFont(42);
  gle1g5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gle1g5->AddEntry(grLambdaGenuine,
                   "Genuine p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)", "l");
  gle1g5->AddEntry(
      grGenuineSmearedSideband,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }(#Lambda#gamma)",
      "l");
  gle1g5->AddEntry(
      grDecompSmearedSideband,
      Form(
          "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }(#Lambda#gamma) (#lambda = %.3f)",
          lambdaplambda),
      "l");
  gle1g5->AddEntry(
      SBmerge,
      "p#minus#kern[-0.65]{ }(#Lambda#gamma) #oplus #bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma)",
      "pez");
  gle1g5->Draw("same");
  g2c3->Print("pLambdaSmearingSideband.pdf");


  outfile->cd();
  histSmear->Write();
  histSmearSideband->Write();
  grLambdaGenuine->Write("Genuine p-Lambda");
  grGenuineSmearedSigma->Write("Genuine p-Lambda smeared Sigma");
  grGenuineSmearedSideband->Write("Genuine p-Lambda smeared sideband");

  grLambdaDecomposition->Write("Decomp p-Lambda");
  grDecompSmearedSigma->Write("Decomp p-Lambda smeared Sigma");
  grDecompSmearedSideband->Write("Decomp p-Lambda smeared sideband");

  grSmearedSigmaForSigma->Write("p-Lambda in Sigma");
  grSmearedSidebandForSigma->Write("p-Lambda in Sigma sideband");
  grDeviationSigma->Write("Deviation");
  SBmerge->Write();
  CATSinput->GetSigmaFile(1)->Write();
  outfile->Close();
}
