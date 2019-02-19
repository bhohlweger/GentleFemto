#include "TidyCats.h"
#include "TRandom3.h"
#include "SidebandSigma.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "TCanvas.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include <iostream>

int main(int argc, char *argv[]) {

  TString InputDir = argv[1];
  TString appendix = argv[2];

  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  auto side = new SidebandSigma();
  side->SetRebin(10);
  side->SetSideBandFile(InputDir.Data(), appendix.Data());
  const double sidebandNormDown = 340;
  const double sidebandNormUp = 440;
  side->SetNormalizationRange(sidebandNormDown, sidebandNormUp);

  side->SideBandCFs();
  auto SBmerge = side->GetSideBands(5);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Get the smearing matrix

  auto file = TFile::Open(Form("%s/AnalysisResults.root", InputDir.Data()));

  TString name = "Sigma0_Femto_";
  name += appendix;
  TDirectory *dir = file->GetDirectory(name);
  name = "femto_";
  name += appendix;
  auto histoList = (TList *) dir->Get(name);
  auto histLambdaGamma = (TH2F*) histoList->FindObject(
      "fHistCorrelationPSigmaPLambda");
  histLambdaGamma->Add(
      (TH2F*) histoList->FindObject(
          "fHistCorrelationAntiPAntiSigmaAntiPAntiLambda"));

  /// Convert to MeV/c
  auto histLambdaGammaMev = new TH2F(Form("%s_MeV", histLambdaGamma->GetName()),
                                     histLambdaGamma->GetTitle(),
                                     histLambdaGamma->GetNbinsX(), 0, 3000,
                                     histLambdaGamma->GetNbinsY(), 0, 3000);

  for (int i = 0; i < histLambdaGamma->GetNbinsX(); ++i) {
    for (int j = 0; j < histLambdaGamma->GetNbinsY(); ++j) {
      histLambdaGammaMev->SetBinContent(j, i,
                                        histLambdaGamma->GetBinContent(i, j));
    }
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendix.Data());
  CATSinput->ObtainCFs(10, 340, 440);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.
  const int binwidth = 8;
  const unsigned NumMomBins = int(800 / binwidth);
  double kMin = SBmerge->GetBinCenter(1) - binwidth / 2.f;
  const double kMax = kMin + binwidth * NumMomBins;

  std::cout << "kMin_pSigma: " << kMin << std::endl;
  std::cout << "kMax_pSigma: " << kMax << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pSigma: " << NumMomBins << std::endl;

  // pp radius systematic variations
  const double ppRadius = 1.32;

  DLM_Ck* Ck_pL = new DLM_Ck(1, 4, NumMomBins, kMin, kMax,
                             Lednicky_SingletTriplet);

  // NLO
  Ck_pL->SetPotPar(0, 2.91);
  Ck_pL->SetPotPar(1, 2.78);
  Ck_pL->SetPotPar(2, 1.54);
  Ck_pL->SetPotPar(3, 2.72);

  // LO
//  Ck_pL->SetPotPar(0, 1.91);
//  Ck_pL->SetPotPar(1, 1.4);
//  Ck_pL->SetPotPar(2, 1.23);
//  Ck_pL->SetPotPar(3, 2.13);

  Ck_pL->SetSourcePar(0, ppRadius);
  Ck_pL->Update();

  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  const double lambdaPurity = 0.9;
  const double lambdaPrimary = 0.619493;
  const Particle lambda(lambdaPurity, lambdaPrimary,
                        { { 1.f - lambdaPrimary } });

  const CATSLambdaParam lambdaParam(proton, lambda);

  DLM_CkDecomposition CkDec_pL("pL", 1, *Ck_pL, histLambdaGammaMev);
  CkDec_pL.AddContribution(
      0, 1. - lambdaParam.GetLambdaParam(CATSLambdaParam::Primary),
      DLM_CkDecomposition::cFake);
  CkDec_pL.Update();

  auto grSideband = new TGraph();
  for (unsigned int i = 0; i < Ck_pL->GetNbins(); ++i) {
    grSideband->SetPoint(i, Ck_pL->GetBinCenter(i),
                         CkDec_pL.EvalCk(Ck_pL->GetBinCenter(i)));
  }

  auto c = new TCanvas();
  SBmerge->Draw();
  SBmerge->SetMaximum(1.7);
  SBmerge->GetXaxis()->SetRangeUser(0, 350);
  grSideband->Draw("L3 same");
  c->Print("Sideband_fit.pdf");

  auto d = new TCanvas();
  histLambdaGamma->Draw("colz");
  histLambdaGamma->GetXaxis()->SetRangeUser(0, 0.5);
  histLambdaGamma->GetYaxis()->SetRangeUser(0, 0.5);
  d->Print("momRes.pdf");

}
