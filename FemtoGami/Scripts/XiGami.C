#include "LambdaGami.h"
#include "CATSLambdaParam.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "SideBandFit.h"
#include "AnalyseProXi.h"
#include "TF1.h"
#include "TSystem.h"

#include "DLM_CkDecomposition.h"

static double cutOff = 1000;  // at this value the calculation and doing of the cf stops

int main(int argc, char *argv[]) {
  const char* fileName = argv[1];
  const char* prefix = argv[2];
  AnalyseProXi* ana = new AnalyseProXi(cutOff, 0.5);
  ana->SetAnalysisFile(fileName, prefix);
  ana->Default();
  ana->GetVariation(0);

  return 0;
}

//TFile* QAOutput = nullptr;

//double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
//                       double Xi1530Var);
//TH1F* BaseLine(TH1F* dataCF);
//TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF, unsigned int varSideNorm);
//TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF);
//void StoreModels(TH1F* unfoldedGenuine, TFile* QAOutput);
//
//  LambdaGami* XiGami = new LambdaGami();
//  MomentumGami* MomGami = new MomentumGami(0.5);
//
//  auto CATSinput = new CATSInput();
//  CATSinput->SetNormalization(0.500, 0.700);
//  CATSinput->SetFixedkStarMinBin(true, 0.);
//  const int rebin = 5;
//
//  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
//  TString SigmaFileName = "Sample6_MeV_compact.root";
//  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
//  CATSinput->SetSigmaFileName(SigmaFileName.Data());
//  CATSinput->ReadSigmaFile();
//  TH2F* momReso = (TH2F*)CATSinput->GetSigmaFile(3)->Clone("MomResolution");
//  momReso->Rebin2D(2);
//  MomGami->SetResolution(momReso, 1000);
//  CATSinput->SetMomentumGami(MomGami);
//
//  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
//  DreamFile->SetAnalysisFile(fileName, prefix, "0");
//
//  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
//  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
//  DreamCF* CFpXiDef = CATSinput->ObtainCFSyst(rebin, "pXiVar0", pXi, ApAXi);
//  TH1F* CFMeasured = CFpXiDef->FindCorrelationFunction(
//      "hCk_UnfoldedpXiVar0MeV_0");
//  CFMeasured = (TH1F*)CFMeasured->Clone("InputCF");
//
//  CFpXiDef->WriteOutput("CFinput.root");
//
//  QAOutput = TFile::Open(
//      TString::Format("%s/debug.root", gSystem->pwd()).Data(), "recreate");
////  CFMeasured->SetDirectory(0);
//  XiGami->StoreStatErr(CFMeasured);
//  QAOutput->cd();
//  CFMeasured->Write();
//
//  double lamGenuine = SetupLambdaPars(XiGami, 1., 1., 1.);
//  //Loop over variations, but set pointer for default to be passed on ?
//
//  //Get rid of the sidebands
//  TH1F* unfoldedSideBand = XimSideband(XiGami, CFMeasured, 1);
//  QAOutput->cd();
//  XiGami->AddStatErr(unfoldedSideBand);
//  unfoldedSideBand->Write();
//
//  //Get rid of the bassline or vaseline
//  TH1F* CFBLFree = BaseLine(unfoldedSideBand);
//  QAOutput->cd();
//  XiGami->AddStatErr(CFBLFree);
//  CFBLFree->Write();
//
//  //Setup the Lambda Parameters
//
//  //Get rid of the p-Xim1530 smeared
//  TH1F* unfoldedFeedDown = Xim1530FeedDown(XiGami, CFBLFree);
//  QAOutput->cd();
//  XiGami->AddStatErr(unfoldedFeedDown);
//  unfoldedFeedDown->Write();
//
//  //Unfold to the genuine CF
//  TH1F* unfoldedGenuine = XiGami->UnfoldGenuine(unfoldedFeedDown, lamGenuine);
//  QAOutput->cd();
//  XiGami->AddStatErr(unfoldedGenuine);
//  unfoldedGenuine->Write();
//
////  StoreModels(unfoldedGenuine, QAOutput);
//
//  QAOutput->Close();


//TH1F* BaseLine(TH1F* dataCF) {
//  TString bslName = TString::Format("%s_woBL", dataCF->GetName()).Data();
//  TH1F* BaseLineFree = (TH1F*) dataCF->Clone(bslName.Data());
//  TF1* funct_1 = new TF1("myPol1", "pol1", 0, cutOff);
//  dataCF->Fit(funct_1, "SN", "", 400, 900);
//  BaseLineFree->Divide(funct_1);
//  QAOutput->cd();
//  funct_1->Write();
//  return BaseLineFree;
//}
//
//TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF, unsigned int varSideNorm) {
//  double normvarCont[3][2];
//  normvarCont[0][0] = 450;
//  normvarCont[0][1] = 650;
//  normvarCont[1][0] = 400;
//  normvarCont[1][1] = 600;
//  normvarCont[2][0] = 500;
//  normvarCont[2][1] = 700;
//
//  SideBandFit* side = new SideBandFit();
//  side->SetSideBandFile("~/cernbox/HM13TeV/AnalysisData/latestSystematic", "HM",
//                        "103", "104");
//  side->SetNormalizationRange(normvarCont[varSideNorm][0],
//                              normvarCont[varSideNorm][1]);
//  side->SideBandCFs(false);
//  TH1F* fitme = side->GetSideBands(5);
//  double SideBandPars[4];
//  side->FitSideBands(fitme, SideBandPars);
//
//  unsigned int nkBin = 0;
//  double kMin, kMax;
//  kMin = dataCF->GetXaxis()->GetXmin();
//  if (cutOff > dataCF->GetXaxis()->GetXmax()) {
//    nkBin = dataCF->FindBin(cutOff);
//    kMax = kMin + nkBin * dataCF->GetBinWidth(1);  //assumes a constant binning
//  } else {
//    kMax = dataCF->GetXaxis()->GetXmax();
//    nkBin = dataCF->GetNbinsX();
//  }
//  TString histName = TString::Format("%sSideBand", dataCF->GetName());
//  TH1F* sideBand = new TH1F(histName.Data(), histName.Data(), nkBin, kMin,
//                            kMax);
//  double* mom = new double(0);
//  for (int iBims = 1; iBims < nkBin + 1; ++iBims) {
//    *mom = sideBand->GetBinCenter(iBims);
//    sideBand->SetBinContent(
//        iBims, SideBandFit::ParameterizationROOT(mom, SideBandPars));
//  }
//  delete mom;
//  QAOutput->cd();
//  sideBand->Write();
//  return XiGami->UnfoldResidual(dataCF, sideBand, XiGami->GetLamdaPar(0));
//}
//
//TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF) {
//  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
//  TString SigmaFileName = "Sample6_MeV_compact.root";
//
//  CATSInput *CATSinput = new CATSInput();
//  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
//  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
//  CATSinput->ReadResFile();
//  CATSinput->SetSigmaFileName(SigmaFileName.Data());
//  CATSinput->ReadSigmaFile();
//
//  TidyCats* tidy = new TidyCats();
//
//  unsigned int nkBin = 0;
//  double kMin, kMax;
//  kMin = dataCF->GetXaxis()->GetXmin();
//  if (cutOff > dataCF->GetXaxis()->GetXmax()) {
//    nkBin = dataCF->FindBin(cutOff);
//    std::cout << "dataCF->GetBinWidth(1): " << dataCF->GetBinWidth(1)
//              << std::endl;
//    kMax = kMin + nkBin * dataCF->GetBinWidth(1);  //assumes a constant binning
//  } else {
//    kMax = dataCF->GetXaxis()->GetXmax();
//    nkBin = dataCF->GetNbinsX();
//  }
//  CATS AB_pXim1530;
//  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, nkBin, kMin, kMax,
//                                 TidyCats::sGaussian);
//
//  AB_pXim1530.SetAnaSource(0, 0.92);  //for now 1.2 fm .. ADJUST!
//  AB_pXim1530.KillTheCat();
//  DLM_Ck* ck = new DLM_Ck(AB_pXim1530.GetNumSourcePars(), 0, AB_pXim1530);
//  DLM_CkDecomposition dec = DLM_CkDecomposition("dummy", 0, *ck, nullptr);
//  DLM_ResponseMatrix* resp = new DLM_ResponseMatrix(*ck, NULL,
//                                                    CATSinput->GetResFile(3),
//                                                    false);
////  smeared = nullptr;
//  DLM_Histo<double>*CkSmeared = new DLM_Histo<double>();
//  CkSmeared->SetUp(1);
//  CkSmeared->SetUp(0, ck->GetNbins(0), ck->GetLowEdge(0), ck->GetUpEdge(0));
//  CkSmeared->Initialize();
//  CkSmeared->SetBinContentAll(0);
//  CkSmeared->SetBinErrorAll(0);
//
//  tidy->Smear(ck, resp, CkSmeared);
//
//  TH1F* pXim1530Converted = tidy->Convert2LesserOf2Evils(CkSmeared, dataCF);
//  pXim1530Converted->SetName("pXim1530Converted");
//  pXim1530Converted->SetTitle("pXim1530Converted");
//  for (int iBin = pXim1530Converted->FindBin(500);
//      iBin < pXim1530Converted->GetNbinsX() + 1; ++iBin) {
//    pXim1530Converted->SetBinContent(iBin, 1);
//  }
//  delete ck;
//  delete resp;
//  QAOutput->cd();
//  pXim1530Converted->Write();
//  return XiGami->UnfoldResidual(dataCF, pXim1530Converted,
//                                XiGami->GetLamdaPar(1));
//}
//
//double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
//                       double Xi1530Var) {
//  //Pro Var: Varies the ratio between the composition of lambdas and Sigma + feeding
//  //to protons
//  double PurityProton = 0.9943;
//
//  double PrimProton = 0.822;
//  double SecLamProton = 0.124;  //Fraction of Lambdas
//
//  double SecFracSigma = 1. - PrimProton - ProVar * SecLamProton;
//  Particle Proton = Particle(PurityProton, PrimProton, { ProVar * SecLamProton,
//                                 SecFracSigma });
//  double PurityXi = 0.915;
//  // Xim Production xseciton: dN/dy = 5.3e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
//  // Xi1530 Production xsection: dN/dy = 2.48e-3 (https://doi.org/10.1140/epjc/s10052-017-4943-1)
//  const double Xi01530XimProdFraction = 1 / 2.;  //Production ratio
//  const double Xim1530XimProdFraction = 1 / 2.;  //Same production ratio
//
//  // 2/3 of Xi0(1530) decays via Xi- + pi+ (Isospin considerations)
//  const double Xi01530Xim_BR = 2 / 3.;
//  // 1/3 of Xi-(1530) decays via Xi- + pi0 (Isospin considerations)
//  const double Xim1530Xim_BR = 1 / 3.;
//
//  // Omega production xsection: dN/dy = 0.67e-3 (https://www.sciencedirect.com/science/article/pii/S037026931200528X)
//  // -> Production Ratio ~ 1/10
//  const double OmegamXimProdFraction = 1 / 10.;
//  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %
//
//  // Produce N Xi's -> Produce:
//  // 1 ) N* 1/10 Omegas -> See N* 1/10 * 8.6% more Xi's
//  // 2)  N* 1/2 Xi0_1530 -> See N*1/2*2/3 = N* 1/3 more Xi's
//  // 3)  N* 1/2 Xim_1530 -> See N*1/2*1/3 = N* 1/6 more Xi's
//  // Total Sample:  N(1+0.0086+1/3+1/6) ->
//  // Primary Fraction = N / N(1+0.0086+1/3+1/6)
//  // Secondary Omegas = N*0.0086  / N(1+0.0086+1/3+1/6)
//  // etc.
//
//  double XiNormalization = 1 + OmegaVar * OmegamXimProdFraction * OmegamXim_BR
//      + Xi1530Var * Xi01530XimProdFraction * Xi01530Xim_BR
//      + Xi1530Var * Xim1530XimProdFraction * Xim1530Xim_BR;
//  double SecOmegaXim = OmegaVar * OmegamXimProdFraction * OmegamXim_BR
//      / (double) XiNormalization;
//  double SecXi01530Xim = Xi1530Var * Xi01530XimProdFraction * Xi01530Xim_BR
//      / (double) XiNormalization;
//  double SecXim1530Xim = Xi1530Var * Xim1530XimProdFraction * Xim1530Xim_BR
//      / (double) XiNormalization;
//  double PrimXim = 1. / (double) XiNormalization;
//  Particle Xi = Particle(PurityXi, PrimXim, { SecOmegaXim, SecXi01530Xim,
//                             SecXim1530Xim });
//
//  CATSLambdaParam lamPar(Proton, Xi);
//  lamPar.PrintLambdaParams();
//  //Sideband first
//  XiGami->SetLambdaPar(
//      lamPar.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::Fake)
//          + lamPar.GetLambdaParam(CATSLambdaParam::Fake,
//                                  CATSLambdaParam::Primary)
//          + lamPar.GetLambdaParam(CATSLambdaParam::Fake,
//                                  CATSLambdaParam::Fake));
//  //p-Xi1530 second
//  XiGami->SetLambdaPar(
//      lamPar.GetLambdaParam(CATSLambdaParam::Primary, CATSLambdaParam::FeedDown,
//                            0, 2));
//  return lamPar.GetLambdaParam(CATSLambdaParam::Primary,
//                               CATSLambdaParam::Primary);
//}
//
//void StoreModels(TH1F* unfoldedGenuine, TFile* QAOutput) {
//  int nBins = unfoldedGenuine->FindBin(520);
//  double xmin = 0;
//  double xmax = 520;
//  TidyCats* tidy = new TidyCats();
//  TidyCats::Sources TheSource = TidyCats::sGaussian;
//  float ppRadii[3];
//  ppRadii[0] = 0.76;
//  ppRadii[1] = 0.79;
//  ppRadii[2] = 0.82;
//
//  CATS Coulomb;
//
//  tidy->GetCatsProtonXiMinus(&Coulomb, nBins, xmin, xmax, TheSource,
//                             TidyCats::pCoulomb, 0);
//  Coulomb.SetAnaSource(0, ppRadii[1]);
//  Coulomb.KillTheCat();
//
//  CATS HAL;
//  tidy->GetCatsProtonXiMinus(&HAL, nBins, xmin, xmax, TheSource,
//                             TidyCats::pHALQCD, 12);
//  HAL.SetAnaSource(0, ppRadii[1]);
//  HAL.KillTheCat();
//
//  CATS ESC16;
//  tidy->GetCatsProtonXiMinus(&ESC16, nBins, xmin, xmax, TheSource,
//                             TidyCats::pRikkenPot, 0);
//  ESC16.SetAnaSource(0, ppRadii[1]);
//  ESC16.KillTheCat();
//
//  TGraph* cou = new TGraph();
//  TGraph* hal = new TGraph();
//  TGraph* esc = new TGraph();
//
//  for (auto it = 0 ; it < unfoldedGenuine->FindBin(xmax); ++it) {
//    double kStar = unfoldedGenuine->GetBinCenter(it+1);
//    cou->SetPoint(it , kStar, Coulomb.GetCorrFun(it));
//    hal->SetPoint(it , kStar, HAL.GetCorrFun(it));
//    esc->SetPoint(it , kStar, ESC16.GetCorrFun(it));
//  }
//  QAOutput->cd();
//  cou->SetName("coulomb");
//  cou->Write();
//  hal->SetName("halqcd");
//  hal->Write();
//  esc->SetName("esc");
//  esc->Write();
//  return;
//}
