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
  int nParticles = atoi(argv[4]);

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
  auto histSmear = (TH2D*) infile->Get("histSmearSigma");
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

  CATS AB_pL;
  tidy->GetCatsProtonLambda(&AB_pL, 41, -9.99, 400, TidyCats::sGaussian,
                            TidyCats::pNLOWF);
  AB_pL.KillTheCat();
  DLM_Ck* Ck_pL = new DLM_Ck(1, 0, AB_pL);
  Ck_pL->SetSourcePar(0, ppRadius);
  Ck_pL->Update();
  auto grSidebandRaw = new TGraph();
  FillWaveGraph(AB_pL, grSidebandRaw);

  const double protonPurity = 0.9943;
  const double protonPrimary = 0.823;
  const double protonLambda = 0.125;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  const double lambdaPurity = 0.946;
  const double lambdaPrimary = 0.595;
  const Particle lambda(lambdaPurity, lambdaPrimary,
                        { { 1.f - lambdaPrimary } });

  const CATSLambdaParam lambdaParam(proton, lambda);
  const float lambdaplambda = lambdaParam.GetLambdaParam(
      CATSLambdaParam::Primary);

  std::cout << "considering a primary fraction for the p-L CF of "
            << lambdaplambda << "\n";

  // apply the smearing with the photon, and the momentum smearing
  auto grSideband = GetSmearedCF(GetSmearedCF(grSidebandRaw, histSmear),
                                 CATSinput->GetSigmaFile(1));
  auto grSidebandLambda = new TGraph();

  auto grSidebandSideband = GetSmearedCF(
      GetSmearedCF(grSidebandRaw, histSmearSideband),
      CATSinput->GetSigmaFile(1));
  auto grSidebandSidebandLambda = new TGraph();

  double x, y;
  for (int i = 0; i < 1000; ++i) {
    grSideband->GetPoint(i, x, y);
    grSidebandLambda->SetPoint(i, x, 1 + (y - 1) * lambdaplambda);
  }

  for (int i = 0; i < 1000; ++i) {
    grSidebandSideband->GetPoint(i, x, y);
    grSidebandSidebandLambda->SetPoint(i, x, 1 + (y - 1) * lambdaplambda);
  }

  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            350);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  DreamPlot::SetStyleHisto(histSmear);
  auto c = new TCanvas("c", "c", 650, 550);
  histSmear->Draw("col");

  histSmear->SetTitle(
      "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})");
  c->Print("pLambdaSmearingMatrixSigma.pdf");

  DreamPlot::SetStyleHisto(histSmearSideband);
  auto d = new TCanvas("d", "d", 650, 550);
  histSmearSideband->Draw("col");
  histSmearSideband->SetTitle(
      "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus(#Lambda#gamma)} (MeV/#it{c})");
  d->Print("pLambdaSmearingMatrixSideband.pdf");

  int lineWidth = 2;
  DreamPlot::SetStyleGraph(grSidebandRaw, 20, kGreen + 3, 0.7);
  grSidebandRaw->SetLineWidth(lineWidth);
  grSidebandRaw->SetLineStyle(3);
  DreamPlot::SetStyleGraph(grSideband, 20, kGreen + 3, 0.7);
  grSideband->SetLineWidth(lineWidth);
  grSideband->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSidebandLambda, 20, kGreen + 3, 0.7);
  grSidebandLambda->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(SBmerge, kOpenCircle, kBlue + 3);

  auto gc3 = new TCanvas("CFpL_lambda", "CFpL_lambda", 0, 0, 650, 550);
  gc3->SetRightMargin(right);
  gc3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 3.5);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grSidebandRaw->Draw("L3same");
  grSideband->Draw("L3same");
  grSidebandLambda->Draw("L3same");
  SBmerge->Draw("pezsame");
  auto gleg5 = new TLegend(0.4, 0.6, 0.4 + 0.45, 0.9);
  gleg5->SetBorderSize(0);
  gleg5->SetTextFont(42);
  gleg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gleg5->AddEntry(grSidebandRaw,
                  "Genuine p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)", "l");
  gleg5->AddEntry(
      grSideband,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}",
      "l");
  gleg5->AddEntry(
      grSidebandLambda,
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

  DreamPlot::SetStyleGraph(grSidebandSideband, 20, kGreen + 3, 0.7);
  grSidebandSideband->SetLineWidth(lineWidth);
  grSidebandSideband->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSidebandSidebandLambda, 20, kGreen + 3, 0.7);
  grSidebandSidebandLambda->SetLineWidth(lineWidth);

  auto g2c3 = new TCanvas("CFpL_lambda_SB", "CFpL_lambda_SB", 0, 0, 650, 550);
  g2c3->SetRightMargin(right);
  g2c3->SetTopMargin(top);
  dummyHist->Draw();
  grSidebandRaw->Draw("L3same");
  grSidebandSideband->Draw("L3same");
  grSidebandSidebandLambda->Draw("L3same");
  SBmerge->Draw("pezsame");
  auto gle1g5 = new TLegend(0.4, 0.6, 0.4 + 0.45, 0.9);
  gle1g5->SetBorderSize(0);
  gle1g5->SetTextFont(42);
  gle1g5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gle1g5->AddEntry(grSidebandRaw,
                   "Genuine p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)", "l");
  gle1g5->AddEntry(
      grSideband,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }(#Lambda#gamma)",
      "l");
  gle1g5->AddEntry(
      grSidebandLambda,
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
  grSidebandRaw->Write("p-Lambda raw");
  grSideband->Write("p-Lambda smeared");
  grSidebandLambda->Write("p-Lambda smeared lambda");
  grSidebandSideband->Write("p-Lambda smeared sideband");
  grSidebandSidebandLambda->Write("p-Lambda smeared lambda");
  SBmerge->Write();
  CATSinput->GetSigmaFile(1)->Write();
  outfile->Close();
}
