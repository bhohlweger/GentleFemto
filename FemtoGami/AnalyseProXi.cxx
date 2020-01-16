/*
 * AnalyseProXi.cxx
 *
 *  Created on: 29 Oct 2019
 *      Author: bernhardhohlweger
 */

#include "AnalyseProXi.h"
#include "CATSLambdaParam.h"
#include "TidyCats.h"
#include "CATSInput.h"
#include "SideBandFit.h"

#include "TF1.h"
#include "TSystem.h"
#include "TFile.h"

#include "DLM_CkDecomposition.h"

AnalyseProXi::AnalyseProXi(double cutoff, double smearing)
    : fcutOff(cutoff),
      fLimitCFRange(false),
      fCFLimit(1000),
      fFilename(""),
      fPrefix(""),
      fSuffix(""),
      fQAOutput(nullptr),
      fXiGami(new LambdaGami()),
      fMomGami(new MomentumGami(smearing)),
      fNormVar(1),
      fSideNormVar(1),
      fBaselineVar(1),
      fLamVarProton(1),
      fLamVarOmega(1),
      fLamVarXim1530(1),
      fRadVarXim1530(1),
      fnormalizationLeft(0.4),
      fnormalizationRight(0.6) {
  fNormMin = {1000, 500, 500, 500, 600, 600, 600, 700, 700, 800};
  fNormMax = {1000, 900, 800, 700, 1000, 900, 800, 1000, 900, 1000};
  fSidebandNormMin = {450,400,500};
  fSidebandNormMax = {650,600,700};
  fBLfunct = {"pol0","pol1"};
  fLamVars = {0.8,1.0,1.2};
  fXim1530Rad = {0.95,1.00,1.05};

  TString CalibBaseDir = "/home/schmollweger/cernbox/HM13TeV/AnalysisData/336_pXiMCNano/ResolutionpXi.root";
  TFile* inFile = TFile::Open(CalibBaseDir.Data(), "read");
  if (!inFile) {
    std::cout << "No Infile set, no Momentum resolution set, RIP \n";
  } else {
    TH2F* momReso = (TH2F*) inFile->Get("FiveMeV");
    fMomGami->SetResolution(momReso, 1);
  }
}

AnalyseProXi::~AnalyseProXi() {
  delete fXiGami;
  delete fMomGami;
  // TODO Auto-generated destructor stub
}

TH1F* AnalyseProXi::GetVariation(int varnumber, bool getModels) {
  std::cout << "Getting the CF \n";
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(fFilename, fPrefix, fSuffix);

  DreamDist* pXi = DreamFile->GetPairDistributions(0, 2, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 3, "");
  if (fLimitCFRange) {
    std::cout << " Limiting the range to " << fCFLimit << std::endl;
    ResetLimits(pXi);
    ResetLimits(ApAXi);
  }
  DreamCF* CFpXiDef = ObtainCorrFunction("pXiVar0", pXi, ApAXi);
  TH1F* CFMeasured = CFpXiDef->FindCorrelationFunction(
      "hCk_RebinnedpXiVar0MeV_0");
  if (!CFMeasured) {
    return nullptr;
  }

  CFMeasured = (TH1F*) CFMeasured->Clone("InputCF");
  CFpXiDef->WriteOutput(
      TString::Format("%s/CFinput_Var%u.root", gSystem->pwd(), varnumber).Data());

  fQAOutput = TFile::Open(
      TString::Format("%s/debug_Var%u.root", gSystem->pwd(), varnumber).Data(),
      "recreate");
//  CFMeasured->SetDirectory(0);
  fXiGami->UnSetLambdaPar();
  fXiGami->StoreStatErr(CFMeasured);
  fQAOutput->cd();
  fMomGami->GetQAList()->Write("MomGami", 1);
  CFMeasured->Write();
  std::cout << "Calculating Lambda parameter \n";
  double lamGenuine = SetupLambdaPars(fXiGami, fLamVars[fLamVarProton],
                                      fLamVars[fLamVarOmega],
                                      fLamVars[fLamVarXim1530]);
  //Loop over variations, but set pointer for default to be passed on ?

  std::cout << "Get rid of the sidebands \n";
  TH1F* unfoldedSideBand = XimSideband(fXiGami, CFMeasured);
  fQAOutput->cd();
  fXiGami->AddStatErr(unfoldedSideBand);
  unfoldedSideBand->Write();

  std::cout << "Get rid of the bassline or vaseline \n";
  TH1F* CFBLFree = BaseLine(unfoldedSideBand);
  fQAOutput->cd();
  fXiGami->AddStatErr(CFBLFree);
  CFBLFree->Write();

  //Setup the Lambda Parameters

  std::cout << "Get rid of the p-Xim1530 smeared \n";
  TH1F* unfoldedFeedDown = Xim1530FeedDown(fXiGami, CFBLFree);
  fQAOutput->cd();
  fXiGami->AddStatErr(unfoldedFeedDown);
  unfoldedFeedDown->Write();

  std::cout << "Unfold to the genuine CF \n";
  TH1F* unfoldedGenuine = fXiGami->UnfoldGenuine(unfoldedFeedDown, lamGenuine);
  fQAOutput->cd();
  fXiGami->AddStatErr(unfoldedGenuine);
  TH1F* outputHist = (TH1F*) unfoldedGenuine->Clone("outputhist");
  unfoldedGenuine->Write();
  if (getModels) {
    std::cout << "Getting the Models \n";
    GetCoulomb(unfoldedGenuine);
    GetHalQCD(unfoldedGenuine);
//    GetESC16(unfoldedGenuine);
  }
  fQAOutput->Close();

  return outputHist;
}

void AnalyseProXi::ResetLimits(DreamDist* dist) {
  dist->SetSEDist(
      LimitRange(dist->GetSEDist(), fCFLimit,
                 TString::Format("%sLim", dist->GetSEDist()->GetName()).Data()),
      "ted");
  dist->SetSEMultDist(
      LimitRange(
          dist->GetSEMultDist(), fCFLimit,
          TString::Format("%sLim", dist->GetSEMultDist()->GetName()).Data()),
      "ted");
  dist->SetMEDist(
      LimitRange(dist->GetMEDist(), fCFLimit,
                 TString::Format("%sLim", dist->GetMEDist()->GetName()).Data()),
      "ted");
  dist->SetMEMultDist(
      LimitRange(
          dist->GetMEMultDist(), fCFLimit,
          TString::Format("%sLim", dist->GetMEMultDist()->GetName()).Data()),
      "ted");
}

TH1F* AnalyseProXi::LimitRange(TH1F* hist, double limit, const char* name) {
  std::cout << name << std::endl;
  double xMin = hist->GetXaxis()->GetXmin();
  int upperBin = hist->FindBin(limit);
  std::cout << "hist->GetBinCenter(upperBin): "
            << hist->GetXaxis()->GetBinCenter(upperBin)
            << " hist->GetXaxis()->GetBinUpEdge(upperBin): "
            << hist->GetXaxis()->GetBinUpEdge(upperBin) << std::endl;
  int nBins = upperBin;
  double xMax = hist->GetXaxis()->GetBinUpEdge(upperBin);
  TH1F* outHist = new TH1F(name, name, nBins, xMin, xMax);
  for (int iBin = 1; iBin < upperBin + 1; ++iBin) {
    outHist->SetBinContent(iBin, hist->GetBinContent(iBin));
    outHist->SetBinError(iBin, hist->GetBinError(iBin));
  }
  return outHist;
}

TH2F* AnalyseProXi::LimitRange(TH2F* hist, double limit, const char* name) {
  std::cout << name << std::endl;
  double xMin = hist->GetXaxis()->GetXmin();
  int upperBin = hist->GetXaxis()->FindBin(limit);
  std::cout << "hist->GetBinCenter(upperBin): "
            << hist->GetXaxis()->GetBinCenter(upperBin)
            << " hist->GetXaxis()->GetBinUpEdge(upperBin): "
            << hist->GetXaxis()->GetBinUpEdge(upperBin) << std::endl;
  int nBins = upperBin;
  double xMax = hist->GetXaxis()->GetBinUpEdge(upperBin);
  TH2F* outHist = new TH2F(name, name, nBins, xMin, xMax,
                           hist->GetYaxis()->GetNbins(),
                           hist->GetYaxis()->GetXmin(),
                           hist->GetYaxis()->GetXmax());

  for (int iBinX = 1; iBinX < upperBin + 1; ++iBinX) {
    for (int iBinY = 1; iBinY < hist->GetNbinsY() + 1; ++iBinY) {
      outHist->SetBinContent(iBinX, iBinY, hist->GetBinContent(iBinX, iBinY));
      outHist->SetBinError(iBinX, iBinY, hist->GetBinError(iBinX, iBinY));
    }
  }
  return outHist;
}

DreamCF* AnalyseProXi::ObtainCorrFunction(const char* name, DreamDist* partDist,
                                          DreamDist* APartDist) {
  DreamCF* outCF = new DreamCF();
  DreamPair* pp = new DreamPair("Part", fnormalizationLeft,
                                fnormalizationRight);
  DreamPair* ApAp = new DreamPair("AntiPart", fnormalizationLeft,
                                  fnormalizationRight);
  if (fnormalizationLeft == 0 || fnormalizationRight == 0) {
    std::cout << "Normalization is 0! Bad results incoming! \n";
  }

  pp->SetPair(partDist);
  ApAp->SetPair(APartDist);

  pp->ReweightMixedEvent(pp->GetPair(), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPair(), 0.2, 0.9);

  pp->UnfoldMomentum(pp->GetPairReweighted(0), fMomGami);
  ApAp->UnfoldMomentum(ApAp->GetPairReweighted(0), fMomGami);

  pp->FixShift(pp->GetPairUnfolded(0), ApAp->GetPairUnfolded(0), 0.005, true);
  ApAp->FixShift(ApAp->GetPairUnfolded(0), pp->GetPairUnfolded(0), 0.005, true);

  pp->Rebin(pp->GetPairFixShifted(0), 4);
  ApAp->Rebin(ApAp->GetPairFixShifted(0), 4);

  outCF->SetPairs(pp, ApAp);
  outCF->GetCorrelations(name);
  return outCF;
}

TH1F* AnalyseProXi::BaseLine(TH1F* dataCF) {
  TString bslName = TString::Format("%s_woBL", dataCF->GetName()).Data();
  TH1F* BaseLineFree = (TH1F*) dataCF->Clone(bslName.Data());
  TF1* funct_1 = new TF1("myBaseline", fBLfunct[fBaselineVar].Data(), 0,
                         fcutOff);
  dataCF->Fit(funct_1, "SNR", "", fNormMin[fNormVar], fNormMax[fNormVar]);
  BaseLineFree->Divide(funct_1);
  fQAOutput->cd();
  funct_1->Write();
  return BaseLineFree;
}

TH1F* AnalyseProXi::XimSideband(LambdaGami* XiGami, TH1F* dataCF) {

  SideBandFit* side = new SideBandFit();
  side->SetSideBandFile("~/cernbox/HM13TeV/AnalysisData/Systematics_5MeV",
                        "PXi", "103", "104");
  side->SetNormalizationRange(fSidebandNormMin[fSideNormVar],
                              fSidebandNormMax[fSideNormVar]);
  side->SetRebin(4);
  side->SideBandCFs(false);
  TH1F* fitme = side->GetSideBands(5);
  double SideBandPars[4];
  side->FitSideBands(fitme, SideBandPars);

  unsigned int nkBin = 0;
  double kMin, kMax;
  kMin = dataCF->GetXaxis()->GetXmin();
  if (fcutOff > dataCF->GetXaxis()->GetXmax()) {
    nkBin = dataCF->FindBin(fcutOff);
    kMax = kMin + nkBin * dataCF->GetBinWidth(1);  //assumes a constant binning
  } else {
    kMax = dataCF->GetXaxis()->GetXmax();
    nkBin = dataCF->GetNbinsX();
  }
  TString histName = TString::Format("%sSideBand", dataCF->GetName());
  TH1F* sideBand = new TH1F(histName.Data(), histName.Data(), nkBin, kMin,
                            kMax);
  double* mom = new double(0);
  for (int iBims = 1; iBims < nkBin + 1; ++iBims) {
    *mom = sideBand->GetBinCenter(iBims);
    sideBand->SetBinContent(
        iBims, SideBandFit::ParameterizationROOT(mom, SideBandPars));
  }
  delete mom;
  fQAOutput->cd();
  sideBand->Write();
  return XiGami->UnfoldResidual(dataCF, sideBand, XiGami->GetLamdaPar(0));
}

TH1F* AnalyseProXi::Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF) {
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
  if (fcutOff < dataCF->GetXaxis()->GetXmax()) {
    nkBin = dataCF->FindBin(fcutOff);
    std::cout << "dataCF->GetBinWidth(1): " << dataCF->GetBinWidth(1)
              << std::endl;
    kMax = kMin + nkBin * dataCF->GetBinWidth(1);  //assumes a constant binning
  } else {
    kMax = dataCF->GetXaxis()->GetXmax();
    nkBin = dataCF->GetNbinsX();
  }
  std::cout << "kMin: " << kMin << " kMax: " << kMax << " nkBin: " << nkBin
            << std::endl;
  CATS AB_pXim1530;
  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, nkBin, kMin, kMax,
                                 TidyCats::sGaussian);

  AB_pXim1530.SetAnaSource(0, fXim1530Rad.at(fRadVarXim1530));  //for now 1.2 fm .. ADJUST!
  AB_pXim1530.KillTheCat();
  DLM_Ck* ck = new DLM_Ck(AB_pXim1530.GetNumSourcePars(), 0, AB_pXim1530);
  DLM_CkDecomposition dec = DLM_CkDecomposition("dummy", 0, *ck, nullptr);
  DLM_ResponseMatrix* resp = new DLM_ResponseMatrix(*ck, NULL,
                                                    CATSinput->GetResFile(3),
                                                    false);
//  smeared = nullptr;
  DLM_Histo<double>*CkSmeared = new DLM_Histo<double>();
  CkSmeared->SetUp(1);
  CkSmeared->SetUp(0, ck->GetNbins(0), ck->GetLowEdge(0), ck->GetUpEdge(0));
  CkSmeared->Initialize();
  CkSmeared->SetBinContentAll(0);
  CkSmeared->SetBinErrorAll(0);
  tidy->Smear(ck, resp, CkSmeared);

  TH1F* pXim1530Converted = tidy->Convert2LesserOf2Evils(CkSmeared, dataCF);
  pXim1530Converted->SetName("pXim1530Converted");
  pXim1530Converted->SetTitle("pXim1530Converted");
  for (int iBin = pXim1530Converted->FindBin(500);
      iBin < pXim1530Converted->GetNbinsX() + 1; ++iBin) {
    pXim1530Converted->SetBinContent(iBin, 1);
  }
  delete ck;
  delete resp;
  delete tidy;
  delete CATSinput;
  fQAOutput->cd();
  pXim1530Converted->Write();
  return XiGami->UnfoldResidual(dataCF, pXim1530Converted,
                                XiGami->GetLamdaPar(1));
}

double AnalyseProXi::SetupLambdaPars(LambdaGami* XiGami, double ProVar,
                                     double OmegaVar, double Xi1530Var) {
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

TGraphErrors* AnalyseProXi::GetCoulomb(TH1F* unfoldedGenuine) {
  std::cout << "Coulomb\n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 0.932114;
  ppRadii[2] = 0.979453;

  CATS CoulombUp;
  CATS CoulombLow;

  tidy->GetCatsProtonXiMinus(&CoulombUp, nBins, xmin, xmax, TheSource,
                             TidyCats::pCoulomb, 0);
  tidy->GetCatsProtonXiMinus(&CoulombLow, nBins, xmin, xmax, TheSource,
                             TidyCats::pCoulomb, 0);

  CoulombUp.SetAnaSource(0, ppRadii[0]);
  CoulombUp.KillTheCat();

  CoulombLow.SetAnaSource(0, ppRadii[2]);
  CoulombLow.KillTheCat();

  TGraphErrors* cou = new TGraphErrors();
  cou->SetName("Coulomb");
  for (auto it = 0; it < nBins; ++it) {
    double kStar = CoulombUp.GetMomentum(it);
    double mean = 0.5 * (CoulombUp.GetCorrFun(it) + CoulombLow.GetCorrFun(it));
    double Err = TMath::Abs(mean - CoulombUp.GetCorrFun(it));

    cou->SetPoint(it, kStar, mean);
    cou->SetPointError(it, 0, Err);
  }
  fQAOutput->cd();
  cou->Write();
  return cou;
}

TGraphErrors* AnalyseProXi::GetHalQCD(TH1F* unfoldedGenuine) {
  std::cout << "HAL QCD \n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 0.932114;
  ppRadii[2] = 0.979453;

  CATS HalUp;
  CATS HalLow;

  tidy->GetCatsProtonXiMinus(&HalUp, nBins, xmin, xmax, TheSource,
                             TidyCats::pHALQCD, 11);
  tidy->GetCatsProtonXiMinus(&HalLow, nBins, xmin, xmax, TheSource,
                             TidyCats::pHALQCD, 13);

  HalUp.SetAnaSource(0, ppRadii[0]);
  HalUp.KillTheCat();

  HalLow.SetAnaSource(0, ppRadii[2]);
  HalLow.KillTheCat();

  TGraphErrors* hal = new TGraphErrors();
  hal->SetName("HalQCD");
  for (auto it = 0; it < nBins - 10; ++it) {
    double kStar = HalUp.GetMomentum(it);
    double mean = 0.5 * (HalUp.GetCorrFun(it) + HalLow.GetCorrFun(it));
    double Err = TMath::Abs(mean - HalUp.GetCorrFun(it));
    hal->SetPoint(it, kStar, mean);
    hal->SetPointError(it, 0, Err);
  }
  fQAOutput->cd();
  hal->Write();
  return hal;
}

TGraphErrors* AnalyseProXi::GetESC16(TH1F* unfoldedGenuine) {
  std::cout << "ESC16 \n";
  int nBins = 300;
  double xmin = 0.5;
  double xmax = 300.5;
  TidyCats* tidy = new TidyCats();
  TidyCats::Sources TheSource = TidyCats::sResonance;
  float ppRadii[3];
  ppRadii[0] = 0.886647;
  ppRadii[1] = 0.932114;
  ppRadii[2] = 0.979453;

  CATS ESCUp;
  CATS ESCLow;

  tidy->GetCatsProtonXiMinus(&ESCUp, nBins, xmin, xmax, TheSource,
                             TidyCats::pRikkenPot, 0);
  tidy->GetCatsProtonXiMinus(&ESCLow, nBins, xmin, xmax, TheSource,
                             TidyCats::pRikkenPot, 0);

  ESCUp.SetAnaSource(0, ppRadii[0]);
  ESCUp.KillTheCat();

  ESCLow.SetAnaSource(0, ppRadii[2]);
  ESCLow.KillTheCat();

  TGraphErrors* esc = new TGraphErrors();
  esc->SetName("ESC");
  for (auto it = 0; it < nBins; ++it) {
    double kStar = ESCUp.GetMomentum(it);
    double mean = 0.5 * (ESCUp.GetCorrFun(it) + ESCLow.GetCorrFun(it));
    double Err = TMath::Abs(mean - ESCUp.GetCorrFun(it));

    esc->SetPoint(it, kStar, mean);
    esc->SetPointError(it, 0, Err);
  }
  fQAOutput->cd();
  esc->Write();
  return esc;
}
