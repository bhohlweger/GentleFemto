#include "LambdaGami.h"
#include "CATSLambdaParam.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "DLM_CkDecomposition.h"

static double cutOff = 1000;  // at this value the calculation and doing of the cf stops

double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
                       double Xi1530Var);
TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF);
TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF);
int main(int argc, char *argv[]) {
  LambdaGami* XiGami = new LambdaGami();
  //read histogram
  TH1F* CFMeasured = nullptr;

  //Get rid of the bassline or vaseline

  //Setup the Lambda Parameters
  double lamGenuine = SetupLambdaPars(XiGami, 1., 1., 1.);
  //Get rid of the sidebands
  TH1F* unfoldedSideBand = XimSideband(XiGami, CFMeasured);
  //Unfold for momentum resolution
  TH1F* unfoldedMomRes;
  //Get rid of the p-Xim1530 smeared
  TH1F* unfoldedFeedDown = Xim1530FeedDown(XiGami,unfoldedMomRes);
  //Unfold to the genuine CF
  TH1F* unfoldedGenuine = XiGami->UnfoldGenuine(unfoldedFeedDown, lamGenuine);
  return 0;
}

TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF) {
  return XiGami->UnfoldResidual(dataCF, nullptr, XiGami->GetLamdaPar(0));
}

TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF) {
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  TString SigmaFileName = "Sample6_MeV_compact.root";

  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName(SigmaFileName.Data());
  CATSinput->ReadSigmaFile();

  TidyCats* tidy = new TidyCats();

  unsigned int nkBin = 0;
  double kMin, kMax;
  kMin = dataCF->GetXaxis()->GetXmin();
  if (cutOff > dataCF->GetXaxis()->GetXmax()) {
    nkBin = dataCF->FindBin(cutOff);
    kMax = kMin + nkBin * dataCF->GetBinWidth(1);  //assumes a constant binning
  } else {
    kMax = dataCF->GetXaxis()->GetXmax();
    nkBin = dataCF->GetNbinsX();
  }

  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, nkBin, kMin, kMax,
                                 TidyCats::sGaussian);

  AB_pXim1530.SetAnaSource(0, 0.92);  //for now 1.2 fm .. ADJUST!
  AB_pXim1530.KillTheCat();
  DLM_Ck* ck = new DLM_Ck(AB_pXim1530.GetNumSourcePars(), 0, AB_pXim1530);
  DLM_CkDecomposition dec = DLM_CkDecomposition("dummy", 0, *ck, nullptr);
  DLM_ResponseMatrix* resp = new DLM_ResponseMatrix(AB_pXim1530, NULL,
                                                    CATSinput->GetResFile(3),
                                                    false);
  DLM_Histo<double>* smeared = new DLM_Histo<double>(*ck);
  tidy->Smear(ck, resp, smeared);
  TH1F* pXim1530Converted = tidy->Convert2LesserOf2Evils(smeared, dataCF);
  pXim1530Converted->SetName("pXim1530Converted");
  pXim1530Converted->SetTitle("pXim1530Converted");

  delete ck;
  delete resp;
  return XiGami->UnfoldResidual(dataCF, pXim1530Converted,
                                XiGami->GetLamdaPar(1));
}

double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
                       double Xi1530Var) {
  //Pro Var: Varies the ratio between the composition of lambdas and Sigma + feeding
  //to protons
  double PurityProton = 0.9943;

  double PrimProton = 0.822;
  double SecLamProton = 0.124;  //Fraction of Lambdas

  double SecFracSigma = 1. - PrimProton - ProVar * SecLamProton;
  Particle Proton = Particle(PurityProton, PrimProton, { ProVar * SecLamProton,
                                 SecFracSigma });
  double PurityXi = 0.915;
  // Xim Production xseciton: dN/dy = 5.3e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // Xi1530 Production xsection: dN/dy = 2.48e-3 (https://doi.org/10.1140/epjc/s10052-017-4943-1)
  const double Xi01530XimProdFraction = 1 / 2.;  //Production ratio
  const double Xim1530XimProdFraction = 1 / 2.;  //Same production ratio

  // 2/3 of Xi0(1530) decays via Xi- + pi+ (Isospin considerations)
  const double Xi01530Xim_BR = 2 / 3.;
  // 1/3 of Xi-(1530) decays via Xi- + pi0 (Isospin considerations)
  const double Xim1530Xim_BR = 1 / 3.;

  // Omega production xsection: dN/dy = 0.67e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
  // -> Production Ratio ~ 1/10
  const double OmegamXimProdFraction = 1 / 10.;
  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %

  // Produce N Xi's -> Produce:
  // 1 ) N* 1/10 Omegas -> See N* 1/10 * 8.6% more Xi's
  // 2)  N* 1/2 Xi0_1530 -> See N*1/2*2/3 = N* 1/3 more Xi's
  // 3)  N* 1/2 Xim_1530 -> See N*1/2*1/3 = N* 1/6 more Xi's
  // Total Sample:  N(1+0.0086+1/3+1/6) ->
  // Primary Fraction = N / N(1+0.0086+1/3+1/6)
  // Secondary Omegas = N*0.0086  / N(1+0.0086+1/3+1/6)
  // etc.

  double XiNormalization = 1 + OmegaVar * OmegamXimProdFraction * OmegamXim_BR
      + Xi1530Var * Xi01530XimProdFraction * Xi01530Xim_BR
      + Xi1530Var * Xim1530XimProdFraction * Xim1530Xim_BR;
  double SecOmegaXim = OmegaVar * OmegamXimProdFraction * OmegamXim_BR
      / (double) XiNormalization;
  double SecXi01530Xim = Xi1530Var * Xi01530XimProdFraction * Xi01530Xim_BR
      / (double) XiNormalization;
  double SecXim1530Xim = Xi1530Var * Xim1530XimProdFraction * Xim1530Xim_BR
      / (double) XiNormalization;
  double PrimXim = 1. / (double) XiNormalization;
  Particle Xi = Particle(PurityXi, PrimXim, { SecOmegaXim, SecXi01530Xim,
                             SecXim1530Xim });

  CATSLambdaParam lamPar(Proton, Xi);
  lamPar.PrintLambdaParams();
  //Sideband first
  XiGami->SetLambdaPar(
      lamPar.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
          + lamPar.GetLambdaParam(CATSLambdaParam::Fake,
                                  CATSLambdaParam::Primary)
          + lamPar.GetLambdaParam(CATSLambdaParam::Fake,
                                  CATSLambdaParam::Fake));
  //p-Xi1530 second
  XiGami->SetLambdaPar(
      lamPar.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown,
                            0, 2));
  return lamPar.GetLambdaParam(CATSLambdaParam::Primary,
                               CATSLambdaParam::Primary);
}
