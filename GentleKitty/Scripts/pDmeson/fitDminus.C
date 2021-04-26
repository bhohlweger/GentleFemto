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
#include "TROOT.h"

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
template<typename T>
T getBootstrapFromVec(const std::vector<T> &vec) {
  return vec.at(gRandom->Uniform() * vec.size());
}

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
  params = { {pol3->GetParameter(0), pol3->GetParameter(1), pol3->GetParameter(
          2), pol3->GetParameter(3)}};
  delete pol3;
  delete pol2;
  delete pol1;
}

/// =====================================================================================
TGraphErrors* EvalBootstrap(TNtuple *tuple, TList* debug, TString OutputDir,
                            TString potName) {
  auto grOut = new TGraphErrors(Form("gr_%s", tuple->GetName()));
  int count = 0;
  for (double kstar = kmin; kstar < kmax;) {
    tuple->Draw("cf >> h(10000, 0, 10)",
                Form("TMath::Abs(kstar - %.3f) < 1e-1", kstar), "N");
    auto hist = (TH1F*) gROOT->FindObject("h");
    if (hist->GetEntries() == 0) {
      kstar += binWidth;
      continue;
    }
    double error = hist->GetRMS();

    tuple->Draw("cf >> h2",
                Form("TMath::Abs(kstar - %.3f) < 1e-1 && BootID == 0", kstar));
    auto hist2 = (TH1F*) gROOT->FindObject("h2");
    if (hist2->GetEntries() == 0) {
      kstar += binWidth;
      continue;
    }
    double cfDefault = hist2->GetMean();
    grOut->SetPoint(count, kstar, cfDefault);
    grOut->SetPointError(count, 0, error);
    if (debug) {
      auto c = new TCanvas();
      DreamPlot::SetStyleHisto(hist);
      hist->SetTitle(
          Form("#it{k}* = %.1f MeV/#it{c} ;#Delta C(#it{k}*); Entries", kstar));
      hist->Draw();
      hist->GetXaxis()->SetRangeUser(hist->GetMean() - 10 * hist->GetRMS(),
                                     hist->GetMean() + 10 * hist->GetRMS());
      auto gr = new TGraphErrors();
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1);
      gr->SetMarkerColor(kRed+2);
      gr->SetLineColor(kRed+2);
      gr->SetPoint(0, cfDefault, 0.5 * hist->GetMaximum());
      gr->SetPointError(0, error, 0);
      gr->Draw("pe same");
      c->Print(
          Form("%s/Debug/Debug_%s_%s_%i.pdf", OutputDir.Data(),
               tuple->GetName(), potName.Data(), int(kstar)));
      delete gr;
      delete c;
      debug->Add(hist);
    }

    delete hist;
    delete hist2;

    count++;
    kstar += binWidth;
  }
  return grOut;
}

/// =====================================================================================
void fitDminus(TString InputDir, TString trigger, TString OutputDir,
               int errorVar, int potential) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int nBoot = 1000;

  int nSystVars = -1;
  TString errorName;
  if (errorVar == 0) {
    std::cout << "PROCESSING TOTAL ERRORS\n";
    errorName = "totErr";
    nSystVars = 20;
  } else if (errorVar == 1) {
    std::cout << "PROCESSING STATISTICAL ERRORS ONLY\n";
    errorName = "statErr";
    nSystVars = 0;
  }

  TString potName;
  if (potential == 0) {
    potName = "flat";
  } else if (potential == 1) {
    potName = "Coulomb";
  } else if (potential == 2) {
    potName = "Haidenbauer";
  } else {
    std::cout << "ERROR: Potential not defined \n";
    return;
  }

  DreamPlot::SetStyle();

  TRandom3 rangen(0);

  int nArguments = 11;
  TString varList =
      TString::Format(
          "BootID:ppRadius:primaryContrib:flatContrib:dstarContrib:sidebandContrib:chi2SidebandLeft:chi2SidebandRight:"
          "chi2Local:ndf:nSigma200").Data();
  auto ntResult = new TNtuple("fitResult", "fitResult", varList.Data());
  auto tupleSideband = new TNtuple("sideband", "sideband", "kstar:cf:BootID");
  auto tupleTotalFit = new TNtuple("totalFit", "totalFit", "kstar:cf:BootID");
  auto tupleSidebandLeft = new TNtuple("sidebandLeft", "sidebandLeft",
                                       "kstar:cf:BootID");
  auto tupleSidebandRight = new TNtuple("sidebandRight", "sidebandRight",
                                        "kstar:cf:BootID");

  float ntBuffer[nArguments];
  bool useBaseline = true;

  TString graphfilename = TString::Format("%s/Fit_pDminus_%s_%s.root",
                                          OutputDir.Data(), errorName.Data(),
                                          potName.Data());
  auto param = new TFile(graphfilename, "RECREATE");
  param->mkdir("bootVars");

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

  std::vector<TGraphAsymmErrors> grCFvec, grSBLeftvec, grSBRightvec;

  const double normLower = 1.5;
  const double normUpper = 2.;
  const int rebin = 10;
  TString dataGrName = "Graph_from_hCk_Reweighted_0MeV";
  TString InputFileName = InputDir;
  InputFileName += "/AnalysisResults.root";

  for (int i = 0; i <= nSystVars; ++i) {
    grCFvec.push_back(
        GetCorrelationGraph(InputFileName, "HM_CharmFemto_", Form("%i", i),
                            dataGrName, normLower, normUpper, rebin));
    grSBLeftvec.push_back(
        GetCorrelationGraph(InputFileName, "HM_CharmFemto_SBLeft_",
                            Form("%i", i), dataGrName, normLower, normUpper,
                            rebin));
    grSBRightvec.push_back(
        GetCorrelationGraph(InputFileName, "HM_CharmFemto_SBRight_",
                            Form("%i", i), dataGrName, normLower, normUpper,
                            rebin));
  }

  std::vector<double> startParams;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  std::vector<double> sidebandFitRange = { { 1000, 800, 1200 } };

  nBins = (potential == 0) ? 160 : 80;
  kmin = 55;
  kmax = (potential == 0) ? 905 : 455;
  binWidth = (kmax - kmin) / double(nBins);

  /// Femtoscopic radius systematic variations
  const double rLow = 0.81;
  const double rDefault = 0.89;
  const double rUp = 0.97;
  std::vector<double> sourceSize = { { rDefault, rLow, rUp } };

  /// Lambda parameters

  /// Proton
  const double protonPurity = 0.984;
  const double protonPrimary = 0.862;
  const double protonLambda = 0.0974;
  const std::vector<double> protonSecondary = { { protonLambda
      / (1. - protonPrimary), protonLambda / (1. - protonPrimary) * 0.8,
      protonLambda / (1. - protonPrimary) * 1.2 } };

  // obtain D-meson properties from Fabrizios files
  std::vector<double> DmesonPurity, DmesonPurityErr;
  std::vector<double> Bfeeddown, BfeeddownErr;
  std::vector<double> DstarFeeding, DstarFeedingErr;


  TString fileAppendix;
  for (int i = 0; i <= nSystVars; ++i) {
    if (i == 0) {
      fileAppendix = "centralcuts";
    } else if (i > 0 && i < 6) {
      fileAppendix = "loose_1";
    } else if (i > 5 && i < 11) {
      fileAppendix = "loose_2";
    } else if (i > 10 && i < 16) {
      fileAppendix = "tight_1";
    } else if (i > 15 && i < 21) {
      fileAppendix = "tight_2";
    }
    auto file = TFile::Open(Form("%s/Fractions_masswindow_2sigma_sphericity_0_1_%s.root", InputDir.Data(), fileAppendix.Data()));
    auto histFracStat = (TH1F*)file->Get("hFractions_DmPr");
    auto grFracSyst = (TGraphAsymmErrors*)file->Get("gFractionsSyst_DmPr");
    double purity = histFracStat->GetBinContent(1);
    double purityStatErr = histFracStat->GetBinError(1);
    double puritySystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(0);
    DmesonPurity.push_back(purity);
    DmesonPurityErr.push_back(std::sqrt(purityStatErr*purityStatErr + puritySystErr*puritySystErr));

    double bfeed = histFracStat->GetBinContent(2);
    double bfeedStatErr = histFracStat->GetBinError(2);
    double bfeedSystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(1);
    Bfeeddown.push_back(bfeed);
    BfeeddownErr.push_back(std::sqrt(bfeedStatErr*bfeedStatErr + bfeedSystErr*bfeedSystErr));

    double dstarfeed = histFracStat->GetBinContent(3);
    double dstarfeedStatErr = histFracStat->GetBinError(3);
    double dstarfeedSystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(2);
    DstarFeeding.push_back(dstarfeed);
    DstarFeedingErr.push_back(std::sqrt(dstarfeedStatErr*dstarfeedStatErr + dstarfeedSystErr*dstarfeedSystErr));
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the model, fitter, etc.

  auto tidyCats = new TidyCats();
  CATS cats, catsDstar;
  if (potential == 1) {
    tidyCats->GetCatsProtonDminus(&cats, nBins, kmin, kmax,
                                  TidyCats::pDCoulombOnly,
                                  TidyCats::sGaussian);
    cats.SetAnaSource(0, rDefault);
    cats.KillTheCat();
    DLM_Coulomb = new DLM_Ck(1, 0, cats);
  } else if (potential == 2) {
    tidyCats->GetCatsProtonDminus(&cats, nBins, kmin, kmax,
                                  TidyCats::pDminusHaidenbauer,
                                  TidyCats::sGaussian);
    cats.SetAnaSource(0, rDefault);
    cats.KillTheCat();
  }
  tidyCats->GetCatsProtonDstarminus(&catsDstar, nBins, kmin, kmax,
                                    TidyCats::pDCoulombOnly,
                                    TidyCats::sGaussian);
  catsDstar.SetAnaSource(0, rDefault);
  catsDstar.KillTheCat();

  DLM_Ck *DLM_Coulomb;
  if (potential == 0) {
    DLM_Coulomb = new DLM_Ck(1, 0, nBins, kmin, kmax, Flat_Residual);
    DLM_Coulomb->Update();
  } else {
    DLM_Coulomb = new DLM_Ck(1, 0, cats);
  }
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
    if (iBoot % 10 == 0) {
      std::cout << "\r Processing progress: "
                << double(iBoot) / double(nBoot) * 100.f << "%" << std::flush;
    }

    const double femtoRad =
        (iBoot == 0 || errorVar == 1) ?
            sourceSize.at(0) : getBootstrapFromVec(sourceSize);

    // since purity etc. depends on the cut variation, this needs to be done consistently
    int nSystVar = (iBoot == 0) ? 0 : gRandom->Uniform() * grCFvec.size();
    auto grVarCF = grCFvec.at(nSystVar);
    auto grVarCFSidebandLeft = grSBLeftvec.at(nSystVar);
    auto grVarCFSidebandRight = grSBRightvec.at(nSystVar);
    
    double dMesonPur = DmesonPurity.at(nSystVar);
    double bFeed = Bfeeddown.at(nSystVar);
    double dStarFeed = DstarFeeding.at(nSystVar);

    // now let's sample the uncertainties
    if (iBoot != 0) {
      dMesonPur += DmesonPurityErr.at(nSystVar) * std::round(gRandom->Uniform(-1.4999, 1.4999));
      bFeed += BfeeddownErr.at(nSystVar) * std::round(gRandom->Uniform(-1.4999, 1.4999));
      dStarFeed += DstarFeedingErr.at(nSystVar) * std::round(gRandom->Uniform(-1.4999, 1.4999));
    }
    
    const double DmesonPrimary = 1.f - bFeed - dStarFeed;

    const Particle dmeson(dMesonPur, DmesonPrimary, { { dStarFeed, bFeed } });

    const double lambdaFeed =
        (iBoot == 0 || errorVar == 1) ?
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
        (iBoot == 0) ?
            &grSBLeftvec.at(0) : getBootstrapGraph(&grVarCFSidebandLeft);
    auto grCFBootstrapSidebandRight =
        (iBoot == 0) ?
            &grSBRightvec.at(0) : getBootstrapGraph(&grVarCFSidebandRight);

    const double sidebandRange =
        (iBoot == 0) ?
            sidebandFitRange.at(0) : getBootstrapFromVec(sidebandFitRange);

    auto fitSidebandLeft = new TF1("fitSidebandLeft", "pol3", 0, sidebandRange);
    getImprovedStartParamsPol3(grCFBootstrapSidebandLeft, sidebandRange,
                               startParams);
    fitSidebandLeft->SetParameters(&startParams[0]);

    int workedLeft = grCFBootstrapSidebandLeft->Fit(fitSidebandLeft, "RQ");
    const double chiSqSidebandLeft = fitSidebandLeft->GetChisquare() / double(fitSidebandLeft->GetNDF());

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandLeft,
                                                          0.683);

    auto fitSidebandRight = new TF1("fitSidebandRight", "pol3", 0,
                                    sidebandRange);
    getImprovedStartParamsPol3(grCFBootstrapSidebandRight, sidebandRange,
                               startParams);
    fitSidebandRight->SetParameters(&startParams[0]);

    int workedRight = grCFBootstrapSidebandRight->Fit(fitSidebandRight, "RQ");
    const double chiSqSidebandRight = fitSidebandRight->GetChisquare() / double(fitSidebandRight->GetNDF());

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandRight,
                                                          0.683);

    // stop here when the fit fails!
    if (workedLeft != 0 || workedRight != 0 || chiSqSidebandLeft > 50 || chiSqSidebandRight > 50) {
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

    DLM_CkDecomposition CkDec_Coulomb(
        "pDminusTotal", 3, *DLM_Coulomb,
        (potential == 0) ? nullptr : momentumResolution);  // mom res not necessary for flat
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
    TGraph *grCFRaw = new TGraph();

    double Chi2 = 0;
    double EffNumBins = 0;

    static double mom, dataY, dataErr, theoryX, theoryY;

    auto grCFBootstrap =
        (iBoot == 0) ? &grCFvec.at(0) : getBootstrapGraph(&grVarCF);

    count = 0;
    for (int i = kmin; i < kmax;) {
      double mom = i;
      grCFRaw->SetPoint(count, mom, DLM_Coulomb->Eval(mom));
      grFitTotalFine->SetPoint(count, mom, CkDec_Coulomb.EvalCk(mom));
      tupleTotalFit->Fill(mom, CkDec_Coulomb.EvalCk(mom), iBoot);
      tupleSideband->Fill(mom, totalSideband->Eval(mom), iBoot);
      tupleSidebandLeft->Fill(mom, fitSidebandLeft->Eval(mom), iBoot);
      tupleSidebandRight->Fill(mom, fitSidebandRight->Eval(mom), iBoot);
      i += binWidth;
      count++;
    }

    for (int iPoint = 0; iPoint <= grVarCF.GetN(); ++iPoint) {
      grVarCF.GetPoint(iPoint, mom, dataY);
      if (mom > 1000) {
        break;
      }

      dataErr = grVarCF.GetErrorY(iPoint);

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
    param->mkdir(TString::Format("bootVars/Graph_%i", iBoot));
    param->cd(TString::Format("bootVars/Graph_%i", iBoot));
    grCFBootstrap->Write("dataBootstrap");

    grCFBootstrapSidebandLeft->Write("dataBootstrapSidebandLeft");
    fitSidebandLeft->Write("sidebandFitLeft");
    grSidebandLeft->Write("sidebandErrorLeft");

    grCFBootstrapSidebandRight->Write("dataBootstrapSidebandRight");
    fitSidebandRight->Write("sidebandFitRight");
    grSidebandRight->Write("sidebandErrorRight");

    totalSideband->SetName("totalSideband");
    totalSideband->Write("sidebandWeightedMean");
    grFitTotal->SetName("TotalFit");
    grFitTotal->Write("TotalFit");
    grFitTotalFine->SetName("TotalFitFine");
    grFitTotalFine->Write("fitFineGrain");

    grCFRaw->Write("raw");

    param->cd();
    ntBuffer[0] = iBoot;
    ntBuffer[1] = femtoRad;
    ntBuffer[2] = primaryContrib;
    ntBuffer[3] = flatContrib;
    ntBuffer[4] = pDstarContrib;
    ntBuffer[5] = sidebandContrib;
    ntBuffer[6] = chiSqSidebandLeft;
    ntBuffer[7] = chiSqSidebandRight;
    ntBuffer[8] = Chi2;
    ntBuffer[9] = (float) EffNumBins;
    ntBuffer[10] = nSigma;
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

  grCFvec.at(0).Write("dataDefault");
  grSBLeftvec.at(0).Write("sidebandLeftDefault");
  grSBRightvec.at(0).Write("sidebandRightDefault");

  auto list = new TList();
  auto fitFull = EvalBootstrap(tupleTotalFit, list, OutputDir, potName);
  auto sidebandFull = EvalBootstrap(tupleSideband, nullptr, OutputDir, potName);
  auto sidebandLeft = EvalBootstrap(tupleSidebandLeft, nullptr, OutputDir,
                                    potName);
  auto sidebandRight = EvalBootstrap(tupleSidebandRight, nullptr, OutputDir,
                                     potName);
  fitFull->Write("fitFull");
  sidebandFull->Write("sidebandFull");
  sidebandLeft->Write("sidebandLeft");
  sidebandRight->Write("sidebandRight");

  param->Close();

  delete ntResult;
  delete param;
  delete tidyCats;
  return;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  TString InputDir = argv[1];
  TString trigger = argv[2];
  TString OutputDir = argv[3];
  int errorVar = atoi(argv[4]);
  int potential = atoi(argv[5]);

  fitDminus(InputDir, trigger, OutputDir, errorVar, potential);
}
