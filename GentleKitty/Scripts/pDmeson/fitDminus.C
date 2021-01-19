#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "CATSInputSigma0.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include "CATSLambdaParam.h"
#include "SidebandSigma.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"
#include "TNtuple.h"
#include "DreamSystematics.h"
#include "TVirtualFitter.h"
#include "TSpline.h"

int nBins;
double kmin, kmax, binWidth;

/// =====================================================================================
bool getNiceBinning(TGraphAsymmErrors &grCF, std::vector<double> &xNew) {

  const int nMax = 1000000;
  const int nPoints = grCF.GetN();
  const double *x = grCF.GetX();

  // first iteration - place the bin center just as the mean of the difference of the two points
  for (int i = 1; i < nPoints + 1; ++i) {
    xNew[i] = x[i - 1] + (x[i] - x[i - 1]) / 2.;
  }
  // First and last point has to be set manually
  xNew[0] = xNew[1] - x[0];
  xNew[nPoints] = 2. * xNew[nPoints - 1] - xNew[nPoints - 2];

  double totOffset = 100;
  double *xNew2 = new double(nPoints);

  int n = 0;
  while (totOffset > 0.001) {
    if (n > nMax)
      return false;  // avoid infinite loop!
    totOffset = 0.;
    for (int i = 0; i < nPoints; ++i) {
      double binCenter = xNew[i] + (xNew[i + 1] - xNew[i]) / 2.;
      double offset = x[i] - binCenter;
      totOffset += std::abs(offset);
      xNew[i] += offset;
    }
    ++n;
  }

  return true;
}

/// =====================================================================================
TGraphAsymmErrors GetCorrelationGraph(TString filename, TString appendix,
                                      TString suffix, TString graphName,
                                      const double normLower,
                                      const double normUpper, const int rebin) {
  TGraphAsymmErrors outputGraph;

  ReadDreamFile *DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename.Data(), appendix.Data(), suffix.Data());

  DreamCF *CF_pDminus = new DreamCF();
  DreamPair *pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair *apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  pDminus->ShiftForEmpty(pDminus->GetPair());
  apDplus->ShiftForEmpty(apDplus->GetPair());
  pDminus->FixShift(pDminus->GetPairShiftedEmpty(0),
                    apDplus->GetPairShiftedEmpty(0), apDplus->GetFirstBin());
  apDplus->FixShift(apDplus->GetPairShiftedEmpty(0),
                    pDminus->GetPairShiftedEmpty(0), pDminus->GetFirstBin());
  pDminus->Rebin(pDminus->GetPairFixShifted(0), rebin, true);
  apDplus->Rebin(apDplus->GetPairFixShifted(0), rebin, true);
  pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(0), 0.2, 0.9,
                              pDminus->GetPair());
  apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(0), 0.2, 0.9,
                              apDplus->GetPair());
  CF_pDminus->SetPairs(pDminus, apDplus);
  CF_pDminus->GetCorrelations();

  for (auto it : CF_pDminus->GetCorrelationFunctionGraphs()) {
    TString itName = it->GetName();
    if (graphName == itName) {
      std::cout << it->GetName() << std::endl;
      outputGraph = *it;
    }
  }

  delete apDplus;
  delete pDminus;
  //delete CF_pDminus;
  delete DreamFile;
  return outputGraph;
}

/// =====================================================================================
TGraphErrors* WeightedMean(TGraphErrors *gr1, TGraphErrors *gr2,
                           const double weight1) {
  if (gr1->GetN() != gr2->GetN()) {
    std::cout << "Unequal number of data points - aborting!\n";
    return nullptr;
  }
  auto grOut = (TGraphErrors*) gr1->Clone(
      Form("weightedMean_%i", int(gRandom->Uniform() * 10000.f)));
  static double x1, x2, y1, y2;
  for (int i = 0; i < gr1->GetN(); ++i) {
    gr1->GetPoint(i, x1, y1);
    gr2->GetPoint(i, x2, y2);
    if (std::abs(x1 - x2) > 1e-3) {
      std::cout << "x values are not equal - aborting!\n";
      return nullptr;
    }
    grOut->SetPoint(i, x1, weight1 * y1 + (1. - weight1) * y2);
    grOut->SetPointError(
        i,
        gr1->GetErrorX(i),
        std::sqrt(
            weight1 * weight1 * gr1->GetErrorY(i) * gr1->GetErrorY(i)
                + (1. - weight1) * (1. - weight1) * gr2->GetErrorY(i)
                    * gr2->GetErrorY(i)));
  }
  return grOut;
}
/// =====================================================================================
TGraphErrors* WeightedMean(TSpline *gr1, TSpline *gr2, const double weight1) {
  auto grOut = new TGraphErrors;
  int count = 0;
  for (int i = kmin; i < kmax; ++i) {
    grOut->SetPoint(count++, i,
                    weight1 * gr1->Eval(i) + (1. - weight1) * gr2->Eval(i));
    i += binWidth;
  }
  return grOut;
}

/// =====================================================================================
DLM_Ck* getDLMCk(TGraphErrors *gr) {
  const double *x = gr->GetX();
  DLM_Ck *out = new DLM_Ck(gr->GetN() - 1, x[0], x[gr->GetN() - 1]);
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    out->SetBinContent(i, yVal);
  }
  out->Update();
  return out;
}

/// =====================================================================================
TH2F* TransformToMeV(const TH2F *input) {
  auto histMeV = new TH2F(
      TString::Format("%s_MeV", input->GetName()).Data(), input->GetTitle(),
      input->GetNbinsX(), 1000. * input->GetXaxis()->GetBinLowEdge(1),
      1000. * input->GetXaxis()->GetBinUpEdge(input->GetNbinsX()),
      input->GetNbinsY(), 1000. * input->GetYaxis()->GetBinLowEdge(1),
      1000. * input->GetYaxis()->GetBinUpEdge(input->GetNbinsY()));
  for (int iBinX = 1; iBinX <= input->GetNbinsX(); iBinX++) {
    for (int iBinY = 1; iBinY <= input->GetNbinsY(); iBinY++) {
      histMeV->SetBinContent(iBinX, iBinY, input->GetBinContent(iBinX, iBinY));
    }
  }
  return histMeV;
}

/// =====================================================================================
double getBootstrapFromVec(const std::vector<double> &vec) {
  return vec.at(gRandom->Uniform() * vec.size());
}
double i = 0;

/// =====================================================================================
TGraphAsymmErrors* getBootstrapGraph(TGraphAsymmErrors *gr) {
  auto grOut = (TGraphAsymmErrors*) gr->Clone(
      Form("bootstrap_%s_%i", gr->GetName(),
           int(gRandom->Uniform() * 10000.f)));
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    grOut->SetPoint(i, xVal, gRandom->Gaus(yVal, gr->GetErrorY(i)));
  }
  return grOut;
}

/// =====================================================================================
TGraphErrors* getBootstrapGraph(TGraphErrors *gr) {
  auto grOut = (TGraphErrors*) gr->Clone(
      Form("bootstrap_%s_%i", gr->GetName(),
           int(gRandom->Uniform() * 10000.f)));
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    grOut->SetPoint(i, xVal, gRandom->Gaus(yVal, gr->GetErrorY(i)));
  }
  return grOut;
}

/// =====================================================================================
void getImprovedStartParamsPol3(TGraphAsymmErrors *gr, const double range,
                                std::vector<double> &params) {
  auto pol1 = new TF1("pol1", "pol1", 0, range);
  gr->Fit(pol1, "RQ");
  auto pol2 = new TF1("pol2", "pol2", 0, range);
  pol2->SetParameter(0, pol1->GetParameter(0));
  pol2->SetParameter(1, pol1->GetParameter(1));
  gr->Fit(pol2, "RQ");
  auto pol3 = new TF1("pol3", "pol3", 0, range);
  pol3->FixParameter(0, pol2->GetParameter(0));
  pol3->FixParameter(1, pol2->GetParameter(1));
  pol3->FixParameter(2, pol2->GetParameter(2));
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(2);
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(1);
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(0);
  gr->Fit(pol3, "RQ");
  params = { { pol3->GetParameter(0), pol3->GetParameter(1), pol3->GetParameter(
      2), pol3->GetParameter(3) } };
  delete pol3;
  delete pol2;
  delete pol1;
}

/// =====================================================================================
void getImprovedStartParamsGaussPol1(TGraphAsymmErrors *gr, const double range,
                                     std::vector<double> &params) {
  auto pol1 = new TF1("pol1", "pol1", 200, range);
  gr->Fit(pol1, "RQ");
  auto total = new TF1("total", "gaus(0) + pol1(3)", 0, range);
  total->SetParameter(0, 1);
  total->FixParameter(1, 0);
  total->SetParameter(2, 100);
  total->FixParameter(3, pol1->GetParameter(0));
  total->FixParameter(4, pol1->GetParameter(1));
  gr->Fit(total, "RQ");
  total->ReleaseParameter(3);
  gr->Fit(total, "RQ");
  total->ReleaseParameter(4);
  gr->Fit(total, "RQ");
  total->ReleaseParameter(1);
  total->SetParLimits(1, -25, 25);
  gr->Fit(total, "RQ");
  params = { { total->GetParameter(0), total->GetParameter(1), total
      ->GetParameter(2), total->GetParameter(3), total->GetParameter(4) } };
  delete total;
  delete pol1;
}

/// =====================================================================================
void fitDminus(TString InputDir, TString trigger, TString OutputDir) {
  const int nBoot = 1000;

  DreamPlot::SetStyle();

  TRandom3 rangen(0);

  int nArguments = 9;
  TString varList = TString::Format(
      "BootID:ppRadius:primaryContrib:flatContrib:dstarContrib:sidebandContrib:"
      "chi2Local:ndf:nSigma200").Data();
  auto ntResult = new TNtuple("fitResult", "fitResult", varList.Data());
  auto tupleSideband = new TNtuple("sideband", "sideband", "kstar:cf");
  auto tupleTotalFit = new TNtuple("totalFit", "totalFit", "kstar:cf");

  float ntBuffer[nArguments];
  bool useBaseline = true;

  TString graphfilename = TString::Format("%s/Fit_pDminus.root",
                                          OutputDir.Data());
  auto param = new TFile(graphfilename, "RECREATE");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  auto calibFile = TFile::Open(
      TString::Format("%s/dstar.root", CalibBaseDir.Data()));
  auto decayKindematicsDstar = (TH2F*) calibFile->Get("histSmearDmeson");

  auto momResFile = TFile::Open(
      TString::Format("%s/momRes_Dmesons.root", CalibBaseDir.Data()).Data());
  auto histMomentumResolutionpDplus = TransformToMeV(
      (TH2F*) momResFile->Get("pDplus"));
  auto histMomentumResolutionpDminus = TransformToMeV(
      (TH2F*) momResFile->Get("pDminus"));
  // since the resolution is the same for p-D+ and p-D- we add them
  auto momentumResolution = histMomentumResolutionpDplus;
  momentumResolution->Add(histMomentumResolutionpDminus);

  std::vector < TGraphAsymmErrors > grCFvec;

  const double normLower = 0.8;
  const double normUpper = 1.;
  const int rebin = 10;
  TString dataGrName = "Graph_from_hCk_Reweighted_0MeV";
  TString InputFileName = InputDir;
  InputFileName += "/AnalysisResults.root";
  TString systVar = "0";

  auto grCF = GetCorrelationGraph(InputFileName, trigger, systVar, dataGrName,
                                  normLower, normUpper, rebin);
  auto grCFSidebandLeft = GetCorrelationGraph(InputFileName, trigger, "7",
                                              dataGrName, normLower, normUpper,
                                              rebin);
  auto grCFSidebandRight = GetCorrelationGraph(InputFileName, trigger, "3",
                                               dataGrName, normLower, normUpper,
                                               rebin);
  std::vector<double> startParams;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  std::vector<double> sidebandFitRange = { { 700, 600, 800 } };

  nBins = 80;
  kmin = 50;
  kmax = 450;
  binWidth = (kmax - kmin) / double(nBins);

  /// Femtoscopic radius systematic variations
  const double rCoreLow = 0.74;
  const double rCoreDefault = 0.81;
  const double rCoreUp = 0.89;
  std::vector<double> sourceSize = { { rCoreDefault, rCoreLow, rCoreUp } };

  /// Lambda parameters

  /// Proton
  const double protonPurity = 0.984;
  const double protonPrimary = 0.862;
  const double protonLambda = 0.0974;
  const std::vector<double> protonSecondary = { { protonLambda
      / (1. - protonPrimary), protonLambda / (1. - protonPrimary) * 0.8,
      protonLambda / (1. - protonPrimary) * 1.2 } };

  /// D-
  const double DmesonPurity = 0.65;   // rough estimate, to be improved

  // From syst. files, evaluated at pT = 2.8 GeV/c
  // Beauty feeding: 0.0578102 +/- 0.0332809
  const std::vector<double> Bfeeddown = {
      { 0.058, 0.058 - 0.033, 0.058 + 0.033 } };

  // From syst. files, evaluated at pT = 2.8 GeV/c
  // D* feeding: 0.287128 +/- 0.0263223
  const std::vector<double> DstarFeeding = { { 0.287, 0.287 - 0.026, 0.287
      + 0.026 } };

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the model, fitter, etc.

  auto tidyCats = new TidyCats();
  CATS catsCoulombOnly, catsDstar;
  tidyCats->GetCatsProtonDminus(&catsCoulombOnly, nBins, kmin, kmax,
                                TidyCats::pCoulombOnly, TidyCats::sResonance);
  tidyCats->GetCatsProtonDstarminus(&catsDstar, nBins, kmin, kmax,
                                    TidyCats::pCoulombOnly,
                                    TidyCats::sResonance);

  catsCoulombOnly.SetAnaSource(0, rCoreDefault);
  catsCoulombOnly.KillTheCat();
  catsDstar.SetAnaSource(0, rCoreDefault);
  catsDstar.KillTheCat();

  auto DLM_Coulomb = new DLM_Ck(1, 0, catsCoulombOnly);
  auto DLM_pDstar = new DLM_Ck(1, 0, catsDstar);

  auto grSidebandLeft = new TGraphErrors(nBins);
  auto grSidebandRight = new TGraphErrors(nBins);
  int count = 0;
  for (double i = kmin; i <= kmax;) {
    grSidebandLeft->SetPoint(count, i, 0);
    grSidebandRight->SetPoint(count, i, 0);
    i += binWidth;
    count++;
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Systematic variations

  // first iteration is the default
  for (int iBoot = 0; iBoot < nBoot;) {
    std::cout << "\r Processing progress: " << iBoot << std::flush;

    const double femtoRad =
        (iBoot == 0) ? sourceSize.at(0) : getBootstrapFromVec(sourceSize);
    const double bFeed =
        (iBoot == 0) ? Bfeeddown.at(0) : getBootstrapFromVec(Bfeeddown);
    const double dStarFeed =
        (iBoot == 0) ? DstarFeeding.at(0) : getBootstrapFromVec(DstarFeeding);

    const double DmesonPrimary = 1.f - bFeed - dStarFeed;

    const Particle dmeson(DmesonPurity, DmesonPrimary,
                          { { dStarFeed, bFeed } });

    const double lambdaFeed =
        (iBoot == 0) ?
            protonSecondary.at(0) : getBootstrapFromVec(protonSecondary);
    const Particle proton(
        protonPurity,
        protonPrimary,
        { { (1. - protonPrimary) * lambdaFeed, (1. - protonPrimary)
            * (1 - lambdaFeed) } });

    const CATSLambdaParam lambdaParam(proton, dmeson);

    const double primaryContrib = lambdaParam.GetLambdaParam(
        CATSLambdaParam::Primary);
    const double pDstarContrib = lambdaParam.GetLambdaParam(
        CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0);
    double sidebandContrib = lambdaParam.GetLambdaParam(
        CATSLambdaParam::Primary, CATSLambdaParam::Fake, 0, 0);
    sidebandContrib += lambdaParam.GetLambdaParam(CATSLambdaParam::FeedDown,
                                                  CATSLambdaParam::Fake, 0, 0);
    sidebandContrib += lambdaParam.GetLambdaParam(CATSLambdaParam::FeedDown,
                                                  CATSLambdaParam::Fake, 1, 0);
    sidebandContrib += lambdaParam.GetLambdaParam(CATSLambdaParam::Fake,
                                                  CATSLambdaParam::Fake);
    const double flatContrib = 1.f - primaryContrib - sidebandContrib
        - pDstarContrib;

    if (iBoot == 0) {
      std::cout << "Lambda parameters for p-D-\n";
      std::cout << " Primary  " << primaryContrib * 100. << "\n";
      std::cout << " Sideband " << sidebandContrib * 100. << "\n";
      std::cout << " p-D*     " << pDstarContrib * 100. << "\n";
      std::cout << " Flat     " << flatContrib * 100. << "\n";
    }

    auto grCFBootstrapSidebandLeft =
        (iBoot == 0) ? &grCFSidebandLeft : getBootstrapGraph(&grCFSidebandLeft);
    auto grCFBootstrapSidebandRight =
        (iBoot == 0) ?
            &grCFSidebandRight : getBootstrapGraph(&grCFSidebandRight);

    const double sidebandRange =
        (iBoot == 0) ?
            sidebandFitRange.at(0) : getBootstrapFromVec(sidebandFitRange);

//    auto fitSidebandLeft = new TF1("fitSidebandLeft", "pol3", 0, sidebandRange);
//    getImprovedStartParamsPol3(grCFBootstrapSidebandLeft, sidebandRange,
//                               startParams);
    auto fitSidebandLeft = new TF1("fitSidebandLeft", "gaus(0) + pol1(3)", 0,
                                   sidebandRange);
    getImprovedStartParamsGaussPol1(grCFBootstrapSidebandLeft, sidebandRange,
                                    startParams);

    fitSidebandLeft->SetParameters(&startParams[0]);
    fitSidebandLeft->SetParLimits(0, -1e4, 1e4);
    fitSidebandLeft->SetParLimits(1, -25, 25);

    int workedLeft = grCFBootstrapSidebandLeft->Fit(fitSidebandLeft, "RQ");

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandLeft,
                                                          0.683);

//    auto fitSidebandRight = new TF1("fitSidebandRight", "pol3", 0,
//                                    sidebandRange);
//    getImprovedStartParamsPol3(grCFBootstrapSidebandRight, sidebandRange,
//                               startParams);
    auto fitSidebandRight = new TF1("fitSidebandRight", "gaus(0) + pol1(3)", 0,
                                    sidebandRange);
    getImprovedStartParamsGaussPol1(grCFBootstrapSidebandRight, sidebandRange,
                                    startParams);

    fitSidebandRight->SetParameters(&startParams[0]);
    fitSidebandRight->SetParLimits(0, -1e4, 1e4);
    fitSidebandRight->SetParLimits(1, -25, 25);
    int workedRight = grCFBootstrapSidebandRight->Fit(fitSidebandRight, "RQ");
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandRight,
                                                          0.683);

    // stop here when the fit fails!
    if (workedLeft != 0 || workedRight != 0) {
      std::cout << "Sideband parametrization failed - repeating\n";
      delete fitSidebandLeft;
      delete fitSidebandRight;
      if (iBoot != 0) {
        delete grCFBootstrapSidebandLeft;
        delete grCFBootstrapSidebandRight;
      }
      continue;
    }

    auto totalSideband = WeightedMean(grSidebandLeft, grSidebandRight, 0.51);  // weight extracted from background integral

    DLM_Coulomb->SetSourcePar(0, femtoRad);
    DLM_Coulomb->Update();
    DLM_pDstar->SetSourcePar(0, femtoRad);
    DLM_pDstar->Update();
    auto DLM_sideband = getDLMCk(totalSideband);

    DLM_CkDecomposition CkDec_Coulomb("pDminusTotal", 3, *DLM_Coulomb,
                                      momentumResolution);
    DLM_CkDecomposition CkDec_pDstar("pDstarminus", 0, *DLM_pDstar, nullptr);
    DLM_CkDecomposition CkDec_sideband("sideband", 0, *DLM_sideband, nullptr);

    CkDec_Coulomb.AddContribution(0, primaryContrib,
                                  DLM_CkDecomposition::cFeedDown, &CkDec_pDstar,
                                  decayKindematicsDstar);
    CkDec_Coulomb.AddContribution(1, flatContrib, DLM_CkDecomposition::cFake);
    CkDec_Coulomb.AddContribution(2, sidebandContrib,
                                  DLM_CkDecomposition::cFeedDown,
                                  &CkDec_sideband);
    CkDec_Coulomb.Update();

    TGraph *grFitTotal = new TGraph();
    TGraph *grFitTotalFine = new TGraph();

    double Chi2 = 0;
    double EffNumBins = 0;

    static double mom, dataY, dataErr, theoryX, theoryY;

    auto grCFBootstrap = (iBoot == 0) ? &grCF : getBootstrapGraph(&grCF);

    count = 0;
    for (int i = kmin; i < kmax;) {
      double mom = i;
      grFitTotalFine->SetPoint(count, mom, CkDec_Coulomb.EvalCk(mom));
      tupleTotalFit->Fill(mom, CkDec_Coulomb.EvalCk(mom));
      tupleSideband->Fill(mom, totalSideband->Eval(mom));
      i += binWidth;
      count++;
    }

    for (int iPoint = 0; iPoint <= grCF.GetN(); ++iPoint) {
      grCFBootstrap->GetPoint(iPoint, mom, dataY);
      if (mom > 1000) {
        break;
      }

      dataErr = grCFBootstrap->GetErrorY(iPoint);

      theoryY = CkDec_Coulomb.EvalCk(mom);
      grFitTotal->SetPoint(iPoint, mom, theoryY);

      if (mom < 200) {
        Chi2 += (dataY - theoryY) * (dataY - theoryY) / (dataErr * dataErr);
        ++EffNumBins;
      }
    }

    double pval = TMath::Prob(Chi2, round(EffNumBins));
    double nSigma = TMath::Sqrt(2) * TMath::ErfcInverse(pval);

    param->cd();
    param->mkdir(TString::Format("Graph_%i", iBoot));
    param->cd(TString::Format("Graph_%i", iBoot));
    grCFBootstrap->Write("dataBootstrap");

//    if (iBoot == 0) {
    grCFBootstrapSidebandLeft->Write("dataBootstrapSidebandLeft");
    fitSidebandLeft->Write("sidebandFitLeft");
    grSidebandLeft->Write("sidebandErrorLeft");

    grCFBootstrapSidebandRight->Write("dataBootstrapSidebandRight");
    fitSidebandRight->Write("sidebandFitRight");
    grSidebandRight->Write("sidebandErrorRight");

//    }

    totalSideband->Write("sidebandWeightedMean");
    grFitTotal->Write("TotalFit");
    grFitTotalFine->Write("fitFineGrain");

    param->cd();
    ntBuffer[0] = iBoot;
    ntBuffer[1] = femtoRad;
    ntBuffer[2] = primaryContrib;
    ntBuffer[3] = flatContrib;
    ntBuffer[4] = pDstarContrib;
    ntBuffer[5] = sidebandContrib;
    ntBuffer[6] = Chi2;
    ntBuffer[7] = (float) EffNumBins;
    ntBuffer[8] = nSigma;
    ntResult->Fill(ntBuffer);

    delete fitSidebandLeft;
    delete fitSidebandRight;

    if (iBoot != 0) {
      delete grCFBootstrapSidebandLeft;
      delete grCFBootstrapSidebandRight;
      delete grCFBootstrap;
    }

    delete DLM_sideband;

    delete totalSideband;
    delete grFitTotal;
    delete grFitTotalFine;

    ++iBoot;
  }

/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Exit through the gift shop...

  param->cd();
  ntResult->Write();
  tupleTotalFit->Write();
  tupleSideband->Write();

  grCF.Write("dataDefault");
  grCFSidebandLeft.Write("sidebandLeftDefault");
  grCFSidebandRight.Write("sidebandRightDefault");

  param->Close();

  delete ntResult;
  delete param;
  delete tidyCats;
  return;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  TString InputDir = argv[1];
  TString trigger = argv[2];
  TString OutputDir = argv[3];

  fitDminus(InputDir, trigger, OutputDir);
}
