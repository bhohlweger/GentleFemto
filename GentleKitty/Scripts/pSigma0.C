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

/// Number of parameters for the sideband fit
const int nSidebandPars = 6;

/// =====================================================================================
/// Fit for the sidebands
auto sidebandFit =
    [ ] (double *x, double *p) {
      return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] *x[0] + std::exp(p[4] + p[5] * x[0]);
    };

/// =====================================================================================
/// Function to cast the nice lambda from above to CATS...
double sidebandFitCATS(const double &Momentum, const double *SourcePar,
                       const double *PotPar) {
  double *x = new double[1];
  x[0] = Momentum;
  double *p = const_cast<double*>(PotPar);
  return sidebandFit(x, p);
}

/// =====================================================================================
void PrintVars(const std::vector<double> &vec) {
  for (const auto &it : vec) {
    std::cout << " " << it;
  }
  std::cout << "\n";
}

/// =====================================================================================
void FitSigma0(const unsigned& NumIter, TString InputDir, TString SystInputDir,
               TString trigger, TString suffix, TString OutputDir,
               const int potential, std::vector<double> params) {
  bool batchmode = true;
  bool debugPlots = false;
  double d0, REf0inv, IMf0inv, deltap0, deltap1, deltap2, etap0, etap1, etap2;

  DreamPlot::SetStyle();
  bool fastPlot = (NumIter == 0) ? true : false;
  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  int nArguments = 34;
  TString varList =
      TString::Format(
          "IterID:femtoFitRange:BLSlope:ppRadius:bl_a:bl_a_err:bl_b:bl_b_err:"
          "sb_p0:sb_p0_err:sb_p1:sb_p1_err:sb_p2:sb_p2_err:sb_p3:sb_p3_err:sb_p4:sb_p4_err:sb_p5:sb_p5_err:"
          "primaryContrib:fakeContrib:SBnormDown:SBnormUp:"
          "chi2NDFGlobal:pvalGlobal:chi2Local:ndf:chi2NDF:pval:nSigma250:nSigma200:nSigma150:CFneg")
          .Data();
  if (potential == 0) {
    varList += ":d0:REf0inv:IMf0inv";
    nArguments += 3;
  } else if (potential == 1) {
    varList += ":deltap0:deltap1:deltap2:etap0:etap1:etap2";
    nArguments += 6;
  }

  TNtuple* ntResult = new TNtuple("fitResult", "fitResult", varList.Data());

  Float_t ntBuffer[nArguments];
  int iterID = 0;
  bool useBaseline = true;

  TString graphfilename;
  if (potential == 0) {
    if (params.size() != 3) {
      std::cout << "ERROR: Wrong number of scattering parameters\n";
      return;
    }
    d0 = params[0];
    REf0inv = params[1];
    IMf0inv = params[2];

    graphfilename = TString::Format("%s/Param_pSigma0_%i_%.3f_%.3f_%.3f.root",
                                    OutputDir.Data(), potential, d0, REf0inv,
                                    IMf0inv);
  } else if (potential == 1) {
    if (params.size() != 6) {
      std::cout << "ERROR: Wrong number of parameters for delta/eta\n";
      return;
    }
    deltap0 = params[0];
    deltap1 = params[1];
    deltap2 = params[2];
    etap0 = params[3];
    etap1 = params[4];
    etap2 = params[5];

    graphfilename = TString::Format(
        "%s/Param_pSigma0_%i_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.root",
        OutputDir.Data(), potential, deltap0, deltap1, deltap2, etap0, etap1,
        etap2);

  } else {
    graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                    OutputDir.Data(), potential);
  }
  auto param = new TFile(graphfilename, "RECREATE");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  std::vector<TH1F*> histSysVar;
  std::vector<float> puritySigma0;

  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                       suffix.Data());
  CATSinput->CountPairs(InputDir.Data(), trigger.Data(), suffix.Data());
  CATSinput->ObtainCFs(10, 250, 400);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHist = CATSinput->GetCF("pSigma0", dataHistName.Data());
  if (!dataHist) {
    std::cerr << "ERROR pSigma0 fitter: p-Sigma0 histogram missing\n";
    return;
  }
  histSysVar.push_back(dataHist);
  puritySigma0.push_back(CATSinput->GetSigma0Purity());

  auto sidebandHistUp = CATSinput->GetCF("pSigmaSBUp",
                                         "hCk_ReweightedpSigmaSBUpMeV_0");
  auto sidebandHistLow = CATSinput->GetCF("pSigmaSBLow",
                                          "hCk_ReweightedpSigmaSBLowMeV_0");
  if (!sidebandHistLow || !sidebandHistUp) {
    std::cerr << "ERROR pSigma0 fitter: p-pSigmaSB histogram missing\n";
    return;
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input for the systematic variations

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  for (int i = 1; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadSigma0CorrelationFile(SystInputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->CountPairs(SystInputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    puritySigma0.push_back(CATSinputVar->GetSigma0Purity());
    CATSinputVar->ObtainCFs(10, 250, 400);
    auto dataHistVar = CATSinputVar->GetCF("pSigma0", dataHistName.Data());
    histSysVar.push_back(dataHistVar);
    delete CATSinputVar;
  }

  if (puritySigma0.size() != histSysVar.size()) {
    std::cout << "ERROR: No matching hist/purity found \n";
    return;
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.
  const int binwidth = dataHist->GetBinWidth(1);
  int NumMomBins_pSigma = int(600 / binwidth);
  double kMin_pSigma = dataHist->GetBinCenter(1) - binwidth / 2.;
  double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;
  int NumMomBins_pSigma_draw = 25;
  double kMin_pSigma_draw = -9.99;
  double kMax_pSigma_draw = 490.01;
  const float drawBinWidth = float(kMax_pSigma_draw - kMin_pSigma_draw)
      / float(NumMomBins_pSigma_draw);

  std::cout << "kMin: " << kMin_pSigma << std::endl;
  std::cout << "kMax: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins: " << NumMomBins_pSigma << std::endl;
  std::cout << "DRAWING\n";
  std::cout << "kMin: " << kMin_pSigma_draw << std::endl;
  std::cout << "kMax: " << kMax_pSigma_draw << std::endl;
  std::cout << "NumMomBins: " << NumMomBins_pSigma_draw << std::endl;
  std::cout << "Binwidth : " << drawBinWidth << std::endl;

  // pp radius systematic variations
  const double resonancesRadius = 1.154;
  const double radiusLower = 1.113;
  const double radiusUpper = 1.196;
  const std::vector<double> sourceSize = { { resonancesRadius, radiusLower,
      radiusUpper } };

  // femto fit region systematic variations
  const std::vector<double> femtoFitRegionUp = { { 550, 500, 600 } };

  /// Lambda parameters
  const double protonPurity = 0.9943;
  const double protonPrimary = 0.877;
  const double protonLambda = 0.089;

  // proton secondary contribution systematic variations
  const std::vector<double> protonSecondary = { { protonLambda
      / (1. - protonPrimary), protonLambda / (1. - protonPrimary) * 0.8,
      protonLambda / (1. - protonPrimary) * 1.2 } };
  std::vector<CATSLambdaParam> lambdaParams;

  const double sigmaPrimary = 1.;
  const Particle sigma0default(puritySigma0[0], sigmaPrimary, { { 0 } });
  const Particle protondefault(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary[0], (1. - protonPrimary)
          * (1 - protonSecondary[0]) } });

  const CATSLambdaParam lambdaParamDefault(sigma0default, protondefault);

  const float lambdaParamDefaultPrimary = lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Primary);
  float lambdaParamDefaultSideband = lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::Primary, 0, 0);
  lambdaParamDefaultSideband += lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 0);
  lambdaParamDefaultSideband += lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 1);
  lambdaParamDefaultSideband += lambdaParamDefault.GetLambdaParam(
      CATSLambdaParam::Fake, CATSLambdaParam::Fake);

  std::cout << "Lambda parameters for the default case\n";
  std::cout << " Primary  " << lambdaParamDefaultPrimary << "\n";
  std::cout << " Sideband " << lambdaParamDefaultSideband << "\n";
  std::cout << " Flat     "
            << 1 - lambdaParamDefaultPrimary - lambdaParamDefaultSideband
            << "\n";

  // sideband fit normalization range systematic variation
  const std::vector<double> sidebandNormDown = { { 300, 250, 350 } };
  const std::vector<double> sidebandNormUp = { { 500, 450, 550 } };

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Fit the sideband
  auto side = new SidebandSigma();
  side->SetRebin(10);
  side->SetSideBandFile(InputDir.Data(), trigger.Data(), suffix.Data());

  /// Prefit for the baseline
  auto PrefitDefault = (TH1F*) dataHist->Clone(
      Form("%s_prefit", dataHist->GetName()));
  auto funct_0_Default = new TF1("myPol0", "pol0", 250, 750);
  PrefitDefault->Fit(funct_0_Default, "FSNRMQ");

  TF1* funct_1_Default = new TF1("myPol1", "pol1", 250, 750);
  PrefitDefault->Fit(funct_1_Default, "FSNRMQ");
  gMinuit->SetErrorDef(1);  // 1 corresponds to 1 sigma contour, 4 to 2 sigma
  auto grPrefitContour_Default = (TGraph*) gMinuit->Contour(40, 0, 1);

  std::vector<double> prefit_a_Default;
  prefit_a_Default.emplace_back(funct_0_Default->GetParameter(0));
  prefit_a_Default.emplace_back(funct_1_Default->GetParameter(0));
  prefit_a_Default.emplace_back(
      TMath::MinElement(grPrefitContour_Default->GetN(),
                        grPrefitContour_Default->GetX()));
  prefit_a_Default.emplace_back(
      TMath::MaxElement(grPrefitContour_Default->GetN(),
                        grPrefitContour_Default->GetX()));

  std::vector<double> prefit_b_Default;
  prefit_b_Default.emplace_back(0);
  prefit_b_Default.emplace_back(funct_1_Default->GetParameter(1));
  prefit_b_Default.emplace_back(
      grPrefitContour_Default->Eval(prefit_a_Default[2]));
  prefit_b_Default.emplace_back(
      grPrefitContour_Default->Eval(prefit_a_Default[3]));

  std::cout << "Result of the prefit to constrain the baseline\n";
  for (size_t i = 0; i < prefit_a_Default.size(); ++i) {
    std::cout << i << " a: " << prefit_a_Default[i] << " b: "
              << prefit_b_Default[i] << "\n";
  }

  if (debugPlots) {
    auto c = new TCanvas();
    DreamPlot::SetStyleGraph(grPrefitContour_Default);
    grPrefitContour_Default->SetTitle("; #it{a}; #it{b}");
    grPrefitContour_Default->SetLineStyle(2);
    auto grDots = new TGraph();
    DreamPlot::SetStyleGraph(grDots, 20, kRed + 2);
    for (size_t i = 0; i < prefit_a_Default.size(); ++i) {
      grDots->SetPoint(i, prefit_a_Default[i], prefit_b_Default[i]);
    }
    grPrefitContour_Default->Draw("AL");
    grDots->Draw("PEsame");
    c->Print(Form("%s/Prefit_%i.pdf", OutputDir.Data(), potential));
    delete grDots;
    delete c;
  }

  delete funct_0_Default;
  delete funct_1_Default;
  delete PrefitDefault;

  if (NumIter != 0) {
    std::cout << "\n\nStarting the systematic variations\n";
    std::cout << "Number of systematic variations of the data: "
              << histSysVar.size() << "\n";
    std::cout << "Number of variations of the fit region: "
              << femtoFitRegionUp.size() << "\n";
    PrintVars(femtoFitRegionUp);
    std::cout << "Number of variations of the baseline:   "
              << prefit_a_Default.size() << "\n";
    PrintVars(prefit_a_Default);
    PrintVars(prefit_b_Default);
    std::cout << "Number of variations of the source size : "
              << sourceSize.size() << "\n";
    PrintVars(sourceSize);
    std::cout << "Number of variations of the sideband normalization: "
              << sidebandNormDown.size() << "\n";
    PrintVars(sidebandNormDown);
    PrintVars(sidebandNormUp);
    std::cout << "Number of variations of the lambda param: "
              << protonSecondary.size() << "\n";
    PrintVars(protonSecondary);
    std::cout << "\n";
  }

  // Set up the model, fitter, etc.
  DLM_Ck* Ck_pSigma0;
  DLM_Ck* Ck_pSigma0_draw;
  CATS AB_pSigma0;
  CATS AB_pSigma0_draw;
  if (potential == 0) {  //  Effective Lednicky with scattering parameters
    std::cout << "Running with scattering parameters - d0 = " << d0
              << " fm - Re(f0^-1) = " << REf0inv << " fm^-1 - Im(f0^-1) = "
              << IMf0inv << "\n";
    Ck_pSigma0 = new DLM_Ck(1, 3, NumMomBins_pSigma, kMin_pSigma, kMax_pSigma,
                            ComplexLednicky_Singlet_InvScatLen);
    Ck_pSigma0_draw = new DLM_Ck(1, 3, NumMomBins_pSigma_draw, kMin_pSigma_draw,
                                 kMax_pSigma_draw,
                                 ComplexLednicky_Singlet_InvScatLen);
  } else if (potential == 1) {  // Effective Lednicky with delta/eta
    std::cout << "Running with delta/eta parametrization \n";
    std::cout << "Delta: " << deltap0 << " " << deltap1 << " " << deltap2
              << "\n";
    std::cout << "Eta  : " << etap0 << " " << etap1 << " " << etap2 << "\n";
    Ck_pSigma0 = new DLM_Ck(1, 6, NumMomBins_pSigma, kMin_pSigma, kMax_pSigma,
                            LednickySingletScatAmplitude);
    Ck_pSigma0_draw = new DLM_Ck(1, 6, NumMomBins_pSigma_draw, kMin_pSigma_draw,
                                 kMax_pSigma_draw,
                                 LednickySingletScatAmplitude);
  } else if (potential == 2) {  // Lednicky coupled channel model fss2
    std::cout << "Running with coupled Lednicky \n";
    Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins_pSigma, kMin_pSigma, kMax_pSigma,
                            Lednicky_gauss_Sigma0);
    Ck_pSigma0_draw = new DLM_Ck(1, 0, NumMomBins_pSigma_draw, kMin_pSigma_draw,
                                 kMax_pSigma_draw, Lednicky_gauss_Sigma0);
  } else if (potential == 3) {  // Haidenbauer WF
    std::cout << "Running with the Haidenbauer chiEFT potential \n";
    // Haidenbauer is valid up to 350 MeV, therefore we have to adopt
    double kMax_pSigma_Haidenbauer = kMin_pSigma;
    int NumMomBins_pSigma_Haidenbauer = 0;
    while (kMax_pSigma_Haidenbauer < 348 - binwidth) {
      kMax_pSigma_Haidenbauer += binwidth;
      ++NumMomBins_pSigma_Haidenbauer;
    }
    tidy->GetCatsProtonSigma0(&AB_pSigma0, NumMomBins_pSigma_Haidenbauer,
                              kMin_pSigma, kMax_pSigma_Haidenbauer,
                              TidyCats::sGaussian,
                              TidyCats::pSigma0Haidenbauer);

    tidy->GetCatsProtonSigma0(&AB_pSigma0_draw, NumMomBins_pSigma_draw,
                              kMin_pSigma_draw, kMax_pSigma_draw,
                              TidyCats::sGaussian,
                              TidyCats::pSigma0Haidenbauer);
    AB_pSigma0.KillTheCat();
    AB_pSigma0_draw.KillTheCat();
    Ck_pSigma0 = new DLM_Ck(1, 0, AB_pSigma0);
    Ck_pSigma0_draw = new DLM_Ck(1, 0, AB_pSigma0_draw);
  } else if (potential == 4) {  // ESC16 WF
    std::cout << "Running with the ESC16 potential \n";
    // ESC16 is valid up to 400 MeV, therefore we have to adopt
    double kMax_pSigma_ESC16 = kMin_pSigma;
    int NumMomBins_pSigma_ESC16 = 0;
    while (kMax_pSigma_ESC16 < 398 - binwidth) {
      kMax_pSigma_ESC16 += binwidth;
      ++NumMomBins_pSigma_ESC16;
    }
    tidy->GetCatsProtonSigma0(&AB_pSigma0, NumMomBins_pSigma_ESC16, kMin_pSigma,
                              kMax_pSigma_ESC16, TidyCats::sGaussian,
                              TidyCats::pSigma0ESC16);

    tidy->GetCatsProtonSigma0(&AB_pSigma0_draw, NumMomBins_pSigma_draw,
                              kMin_pSigma_draw, kMax_pSigma_draw,
                              TidyCats::sGaussian, TidyCats::pSigma0ESC16);
    AB_pSigma0.KillTheCat();
    AB_pSigma0_draw.KillTheCat();
    Ck_pSigma0 = new DLM_Ck(1, 0, AB_pSigma0);
    Ck_pSigma0_draw = new DLM_Ck(1, 0, AB_pSigma0_draw);
  }

  float counter = 0;
  float total = histSysVar.size() * femtoFitRegionUp.size()
      * prefit_a_Default.size() * sidebandNormDown.size();
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Systematic variations

  // 1. Systematic variations of the data
  for (size_t systDataIter = 0; systDataIter < histSysVar.size();
      ++systDataIter) {
    auto currentHist = histSysVar[systDataIter];

    /// Prefit for the baseline
    auto Prefit = (TH1F*) currentHist->Clone(
        Form("%i_%s_prefit", int(systDataIter), currentHist->GetName()));
    auto funct_0 = new TF1("myPol0", "pol0", 250, 750);
    Prefit->Fit(funct_0, "FSNRMQ");

    TF1* funct_1 = new TF1("myPol1", "pol1", 250, 750);
    Prefit->Fit(funct_1, "FSNRMQ");
    gMinuit->SetErrorDef(1);  // 1 corresponds to 1 sigma contour, 4 to 2 sigma
    auto grPrefitContour = (TGraph*) gMinuit->Contour(40, 0, 1);

    static std::vector<double> prefit_a;
    prefit_a.clear();
    prefit_a.emplace_back(funct_0->GetParameter(0));
    prefit_a.emplace_back(funct_1->GetParameter(0));
    prefit_a.emplace_back(
        TMath::MinElement(grPrefitContour->GetN(), grPrefitContour->GetX()));
    prefit_a.emplace_back(
        TMath::MaxElement(grPrefitContour->GetN(), grPrefitContour->GetX()));

    static std::vector<double> prefit_b;
    prefit_b.clear();
    prefit_b.emplace_back(0);
    prefit_b.emplace_back(funct_1->GetParameter(1));
    prefit_b.emplace_back(grPrefitContour->Eval(prefit_a[2]));
    prefit_b.emplace_back(grPrefitContour->Eval(prefit_a[3]));

    delete funct_0;
    delete funct_1;
    delete Prefit;

    /// Lambda parameters
    std::vector<CATSLambdaParam> lambdaParams;

    const double sigmaPurity = puritySigma0[systDataIter];
    const Particle sigma0(sigmaPurity, sigmaPrimary, { { 0 } });

    for (size_t lambdaIter = 0; lambdaIter < protonSecondary.size();
        ++lambdaIter) {
      const Particle proton(
          protonPurity,
          protonPrimary,
          { { (1. - protonPrimary) * protonSecondary[lambdaIter], (1.
              - protonPrimary) * (1 - protonSecondary[lambdaIter]) } });

      lambdaParams.push_back( { sigma0, proton });
    }

    // 2. Femto fit range
    for (size_t femtoFitIter = 0; femtoFitIter < femtoFitRegionUp.size();
        ++femtoFitIter) {

      // 3. Baseline treatment
      for (size_t blIter = 0; blIter < prefit_a.size(); ++blIter) {

        if (blIter == 0) {
          useBaseline = false;  //use baseline
        } else {
          useBaseline = true;  // no baseline
        }

        // 4. Sideband normalization
        for (size_t sbNormIter = 0; sbNormIter < sidebandNormDown.size();
            ++sbNormIter) {

          side->SetNormalizationRange(sidebandNormDown[sbNormIter],
                                      sidebandNormUp[sbNormIter]);

          side->SideBandCFs();
          auto SBmerge = side->GetSideBands(5);
          auto sideband = new TF1(Form("sideband_%i", iterID), sidebandFit, 0,
                                  650, nSidebandPars);
          sideband->SetParameter(0, -1.5);
          sideband->SetParameter(1, 0.);
          sideband->SetParameter(2, 0);
          sideband->SetParameter(3, 0.1);
          sideband->SetParameter(4, 1);
          sideband->SetParameter(5, 0);
          SBmerge->Fit(sideband, "FSNRMQ");

          DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins_pSigma,
                                           kMin_pSigma, kMax_pSigma,
                                           sidebandFitCATS);

          DLM_Ck* Ck_SideBand_draw = new DLM_Ck(0, nSidebandPars,
                                                NumMomBins_pSigma_draw,
                                                kMin_pSigma_draw,
                                                kMax_pSigma_draw,
                                                sidebandFitCATS);

          for (unsigned i = 0; i < sideband->GetNumberFreeParameters(); ++i) {
            Ck_SideBand->SetPotPar(i, sideband->GetParameter(i));
            Ck_SideBand_draw->SetPotPar(i, sideband->GetParameter(i));
          }
          Ck_SideBand->Update();
          Ck_SideBand_draw->Update();

          std::cout << "\r Processing progress: "
              << TString::Format("%.1f %%", counter++ / total * 100.f).Data()
              << std::flush;

          // 5. Source size
          for (size_t sizeIter = 0; sizeIter < sourceSize.size(); ++sizeIter) {

            // 6. Lambda parameters
            for (size_t lambdaIter = 0; lambdaIter < lambdaParams.size();
                ++lambdaIter) {

              /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              /// Correlation function
              Ck_pSigma0->SetSourcePar(0, sourceSize[sizeIter]);
              Ck_pSigma0->Update();
              Ck_pSigma0_draw->SetSourcePar(0, sourceSize[sizeIter]);
              Ck_pSigma0_draw->Update();

              DLM_CkDecomposition CkDec_pSigma0("pSigma0", 2, *Ck_pSigma0,
                                                CATSinput->GetSigmaFile(1));

              DLM_CkDecomposition CkDec_pSigma0_draw(
                  "pSigma0Draw", 2, *Ck_pSigma0_draw,
                  CATSinput->GetSigmaFile(1));

              /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

              DLM_CkDecomposition CkDec_SideBand("pSigma0SideBand", 0,
                                                 *Ck_SideBand, nullptr);

              DLM_CkDecomposition CkDec_SideBand_draw("pSigma0SideBandDraw", 0,
                                                      *Ck_SideBand_draw,
                                                      nullptr);
              float sidebandContr = lambdaParams[lambdaIter].GetLambdaParam(
                  CATSLambdaParam::Fake, CATSLambdaParam::Primary, 0, 0);
              sidebandContr += lambdaParams[lambdaIter].GetLambdaParam(
                  CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 0);
              sidebandContr += lambdaParams[lambdaIter].GetLambdaParam(
                  CATSLambdaParam::Fake, CATSLambdaParam::FeedDown, 0, 1);
              sidebandContr += lambdaParams[lambdaIter].GetLambdaParam(
                  CATSLambdaParam::Fake, CATSLambdaParam::Fake);

              const float primaryContr =
                  lambdaParams[lambdaIter].GetLambdaParam(
                      CATSLambdaParam::Primary);
              CkDec_pSigma0.AddContribution(0, sidebandContr,
                                            DLM_CkDecomposition::cFake,
                                            &CkDec_SideBand);
              CkDec_pSigma0.AddContribution(1,
                                            1.f - sidebandContr - primaryContr,
                                            DLM_CkDecomposition::cFeedDown);
              CkDec_pSigma0.Update();

              CkDec_pSigma0_draw.AddContribution(0, sidebandContr,
                                                 DLM_CkDecomposition::cFake,
                                                 &CkDec_SideBand_draw);
              CkDec_pSigma0_draw.AddContribution(
                  1, 1.f - sidebandContr - primaryContr,
                  DLM_CkDecomposition::cFeedDown);
              CkDec_pSigma0_draw.Update();

              /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              /// Fitter
              DLM_Fitter1* fitter = new DLM_Fitter1(1);
              fitter->SetSystem(0, *currentHist, 1, CkDec_pSigma0, kMin_pSigma,
                                femtoFitRegionUp[femtoFitIter], 900, 900);
              fitter->SetSeparateBL(0, false);              //Simultaneous BL
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_a,
                                   prefit_a[blIter]);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_b,
                                   prefit_b[blIter]);
              fitter->AddSameSource("pSigma0SideBand", "pSigma0", 1);

              //Fit BL & Normalization
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_c, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_Cl, -1.);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0,
                                   sourceSize[sizeIter]);

              // Suppress warnings from ROOT
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot0, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot1, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot2, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot3, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot4, 0);
              fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot5, 0);
              if (potential == 0) {
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot0, REf0inv);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot1, IMf0inv);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot2, d0);
              } else if (potential == 1) {
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot0, deltap0);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot1, deltap1);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot2, deltap2);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot3, etap0);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot4, etap1);
                fitter->FixParameter("pSigma0", DLM_Fitter1::p_pot5, etap2);
              }

              fitter->GoBabyGo();

              /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              /// Get the parameters from the fit
              const double bl_a = fitter->GetParameter("pSigma0",
                                                       DLM_Fitter1::p_a);
              const double bl_a_err = fitter->GetParError("pSigma0",
                                                          DLM_Fitter1::p_a);
              const double bl_b = fitter->GetParameter("pSigma0",
                                                       DLM_Fitter1::p_b);
              const double bl_b_err = fitter->GetParError("pSigma0",
                                                          DLM_Fitter1::p_b);
              const double Cl = fitter->GetParameter("pSigma0",
                                                     DLM_Fitter1::p_c);
              const double chi2 = fitter->GetChi2Ndf();
              const double pval = fitter->GetPval();
              const bool isCFneg = fitter->CheckNegativeCk();

              TGraph FitResult_pSigma0;
              FitResult_pSigma0.SetName(
                  TString::Format("pSigma0Graph_%i", iterID));
              FitResult_pSigma0.SetTitle(
                  TString::Format("pSigma0Graph_%i", iterID));
              fitter->GetFitGraph(0, FitResult_pSigma0);

              double Chi2_pSigma0_250 = 0;
              double EffNumBins_pSigma0_250 = 0;
              double Chi2_pSigma0_200 = 0;
              double EffNumBins_pSigma0_200 = 0;
              double Chi2_pSigma0_150 = 0;
              double EffNumBins_pSigma0_150 = 0;
              int maxkStarBin = currentHist->FindBin(250);
              for (unsigned uBin = 1; uBin <= maxkStarBin; uBin++) {

                double mom = currentHist->GetBinCenter(uBin);
                //double dataX;
                double dataY;
                double dataErr;
                double theoryX;
                double theoryY;

                if (mom > femtoFitRegionUp[femtoFitIter]) {
                  continue;
                }

                FitResult_pSigma0.GetPoint(uBin - 1, theoryX, theoryY);
                if (mom != theoryX) {
                  std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX
                            << std::endl;
                  return;
                }
                dataY = currentHist->GetBinContent(uBin);
                dataErr = currentHist->GetBinError(uBin);
                if (mom < 250) {
                  Chi2_pSigma0_250 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  ++EffNumBins_pSigma0_250;
                  if (mom < 200) {
                    Chi2_pSigma0_200 += (dataY - theoryY) * (dataY - theoryY)
                        / (dataErr * dataErr);
                    ++EffNumBins_pSigma0_200;
                    if (mom < 150) {
                      Chi2_pSigma0_150 += (dataY - theoryY) * (dataY - theoryY)
                          / (dataErr * dataErr);
                      ++EffNumBins_pSigma0_150;

                    }
                  }
                }
              }
              double pvalpSigma0_250 = TMath::Prob(
                  Chi2_pSigma0_250, round(EffNumBins_pSigma0_250));
              double nSigmapSigma0_250 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalpSigma0_250);

              double pvalpSigma0_200 = TMath::Prob(
                  Chi2_pSigma0_200, round(EffNumBins_pSigma0_200));
              double nSigmapSigma0_200 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalpSigma0_200);

              double pvalpSigma0_150 = TMath::Prob(
                  Chi2_pSigma0_150, round(EffNumBins_pSigma0_150));
              double nSigmapSigma0_150 = TMath::Sqrt(2)
                  * TMath::ErfcInverse(pvalpSigma0_150);

              /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              /// Write out all the stuff

              TGraph grCFSigmaExtrapolate;
              grCFSigmaExtrapolate.SetName(Form("S0Extrapolate_%i", iterID));
              TGraph grCFSigmaRaw;
              grCFSigmaRaw.SetName(Form("Sigma0Raw_%i", iterID));
              TGraph grCFSigmaMain;
              grCFSigmaMain.SetName(Form("Sigma0Main_%i", iterID));
              TGraph grCFSigmaFeed;
              grCFSigmaFeed.SetName(Form("Sigma0Feed_%i", iterID));
              TGraph grCFSigmaSideband;
              grCFSigmaSideband.SetName(Form("Sigma0Sideband_%i", iterID));
              TGraph grCFSigmaSidebandExtrapolate;
              grCFSigmaSidebandExtrapolate.SetName(
                  Form("SBExtrapolate_%i", iterID));

              for (int i = 0; i < NumMomBins_pSigma_draw; ++i) {
                const float mom = Ck_pSigma0_draw->GetBinCenter(0, i);
                const float baseline = bl_a + bl_b * mom;
                float cf = CkDec_pSigma0_draw.EvalCk(mom);
                grCFSigmaExtrapolate.SetPoint(
                    i, mom, CkDec_pSigma0_draw.EvalCk(mom) * baseline);
              }

              for (int i = 0; i < Ck_SideBand_draw->GetNbins(); ++i) {
                const float mom = Ck_SideBand_draw->GetBinCenter(0, i);
                const float baseline = bl_a + bl_b * mom;
                grCFSigmaSidebandExtrapolate.SetPoint(
                    i,
                    mom,
                    (((Ck_SideBand_draw->Eval(mom) - 1.) * sidebandContr) + 1)
                        * baseline);
              }

              double mom, ck;
              for (unsigned int i = 0; i < FitResult_pSigma0.GetN(); ++i) {
                FitResult_pSigma0.GetPoint(i, mom, ck);
                const float baseline = bl_a + bl_b * mom;
                grCFSigmaRaw.SetPoint(i, mom,
                                      CkDec_pSigma0.EvalCk(mom) * baseline);
                grCFSigmaMain.SetPoint(
                    i,
                    mom,
                    (((CkDec_pSigma0.EvalMain(mom) - 1.) * primaryContr) + 1)
                        * baseline);
                grCFSigmaFeed.SetPoint(
                    i, mom, CkDec_pSigma0.EvalMainFeed(mom) * baseline);
                grCFSigmaSideband.SetPoint(
                    i,
                    mom,
                    (((Ck_SideBand->Eval(mom) - 1.) * sidebandContr) + 1)
                        * baseline);
              }

              param->cd();
              param->mkdir(TString::Format("Graph_%i", iterID));
              param->cd(TString::Format("Graph_%i", iterID));

              /// beautification
              DreamPlot::SetStyleHisto(currentHist, 24, kBlack);
              currentHist->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
              DreamPlot::SetStyleHisto(SBmerge, 20, kRed + 2);
              SBmerge->SetTitle(";#it{k}* (MeV/#it{c}); C_{sideband}(#it{k}*)");
              DreamPlot::SetStyleHisto(sidebandHistLow, 26, kGreen + 2);
              sidebandHistLow->SetTitle(
                  ";#it{k}* (MeV/#it{c}); C_{sideband, low}(#it{k}*)");
              DreamPlot::SetStyleHisto(sidebandHistUp, 26, kCyan + 2);
              sidebandHistUp->SetTitle(
                  ";#it{k}* (MeV/#it{c}); C_{sideband, up}(#it{k}*)");
              sideband->SetLineWidth(2);
              sideband->SetLineStyle(2);
              sideband->SetLineColor(kGray + 1);
              FitResult_pSigma0.SetLineWidth(2);
              FitResult_pSigma0.SetLineColor(kRed + 2);
              grCFSigmaSideband.SetLineWidth(2);
              grCFSigmaSideband.SetLineStyle(2);
              grCFSigmaSideband.SetLineColor(kGray + 1);

              if (fastPlot || iterID == 0) {
                dataHist->Write();
              }

              grCFSigmaExtrapolate.Write();
              grCFSigmaSidebandExtrapolate.Write();
              grCFSigmaSideband.Write();
              FitResult_pSigma0.Write(Form("Fit_%i", iterID));
              currentHist->SetName(
                  TString::Format("HistCF_Var_%i", int(systDataIter)));
              currentHist->Write();

              if (fastPlot || iterID == 0) {
                grCFSigmaRaw.Write();
                grCFSigmaMain.Write();
                grCFSigmaFeed.Write();
                grPrefitContour->Write("fitContour");
                sideband->Write(Form("SidebandFitNotScaled_%i", iterID));
                SBmerge->Write(Form("SidebandMerged_%i", iterID));
                sidebandHistLow->Write("SidebandLow");
                sidebandHistUp->Write("SidebandUp");

                auto c = new TCanvas("DefaultFit", "DefaultFit");
                currentHist->GetXaxis()->SetRangeUser(0., 600);
                currentHist->GetYaxis()->SetRangeUser(0.8, 1.6);
                currentHist->Draw();
                FitResult_pSigma0.Draw("l3same");
                grCFSigmaSideband.Draw("l3same");

                auto info = new TPaveText(0.5, useBaseline ? 0.505 : 0.58, 0.88,
                                          0.85, "blNDC");
                info->SetBorderSize(0);
                info->SetTextSize(0.04);
                info->SetFillColor(kWhite);
                info->SetTextFont(42);
                TString SOURCE_NAME = "Gauss";
                double Yoffset = 1.2;
                info->AddText(
                    TString::Format(
                        "#it{r}_{%s} = %.3f #pm %.3f fm", SOURCE_NAME.Data(),
                        fitter->GetParameter("pSigma0", DLM_Fitter1::p_sor0),
                        fitter->GetParError("pSigma0", DLM_Fitter1::p_sor0)));
                info->AddText(
                    TString::Format("#it{a} = %.3f #pm %.3f", bl_a, bl_a_err));

                if (useBaseline) {
                  info->AddText(
                      TString::Format(
                          "#it{b} = (%.3f #pm %.3f ) #times 10^{-4}",
                          bl_b * 1e4, bl_b_err * 1e4));
                }
                info->AddText(
                    TString::Format(
                        "#chi_{loc}^{2}/ndf=%.1f/%.0f = %.3f", Chi2_pSigma0_250,
                        EffNumBins_pSigma0_250,
                        Chi2_pSigma0_250 / double(EffNumBins_pSigma0_250)));
                info->AddText(
                    TString::Format("#it{p}_{val}=%.3f, n_{#sigma}=%.3f",
                                    pvalpSigma0_250, nSigmapSigma0_250));

                info->Draw("same");
                c->Write("CFplot");
                if (debugPlots) {
                  if (potential == 0) {
                    c->Print(
                        Form("%s/CF_pSigma0_%.3f_%.3f_%.3f.pdf",
                             OutputDir.Data(), d0, REf0inv, IMf0inv));
                  } else if (potential == 1) {
                    c->Print(
                        Form("%s/CF_pSigma0_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.pdf",
                             OutputDir.Data(), deltap0, deltap1, deltap2, etap0,
                             etap1, etap2));
                  } else {
                    c->Print(Form("%s/CF_pSigma0.pdf", OutputDir.Data()));
                  }
                }

                auto d = new TCanvas("SidebandFit", "SidebandFit");
                SBmerge->GetXaxis()->SetRangeUser(0., 600);
                SBmerge->GetYaxis()->SetRangeUser(0.8, 1.6);
                SBmerge->Draw();
                sideband->Draw("l3same");
                d->Write("CFsideband");
                if (debugPlots) {
                  d->Print(Form("%s/CF_pSideband.pdf", OutputDir.Data()));
                }

                auto e = new TCanvas("Sidebands", "Sidebands");
                SBmerge->Draw();
                sideband->Draw("l3same");
                sidebandHistLow->Draw("same");
                sidebandHistUp->Draw("same");
                auto leg = new TLegend(0.5, 0.6, 0.88, 0.85);
                leg->SetTextFont(42);
                leg->SetTextSize(0.05);
                leg->AddEntry(sidebandHistLow, "Sideband low", "pe");
                leg->AddEntry(sidebandHistUp, "Sideband up", "pe");
                leg->AddEntry(SBmerge, "Sideband merged", "pe");
                leg->AddEntry(sideband, "Fit", "l");
                leg->Draw("same");
                e->Write("CFsideband");
                if (debugPlots) {
                  e->Print(Form("%s/CF_pSideband_all.pdf", OutputDir.Data()));
                }

                delete c;
                delete info;
                delete d;
                delete e;
              }

              param->cd();
              ntBuffer[0] = iterID;
              ntBuffer[1] = femtoFitRegionUp[femtoFitIter];
              ntBuffer[2] = (float) blIter;
              ntBuffer[3] = sourceSize[sizeIter];
              ntBuffer[4] = bl_a;
              ntBuffer[5] = bl_a_err;
              ntBuffer[6] = bl_b;
              ntBuffer[7] = bl_b_err;
              ntBuffer[8] = sideband->GetParameter(0);
              ntBuffer[9] = sideband->GetParError(0);
              ntBuffer[10] = sideband->GetParameter(1);
              ntBuffer[11] = sideband->GetParError(1);
              ntBuffer[12] = sideband->GetParameter(2);
              ntBuffer[13] = sideband->GetParError(2);
              ntBuffer[14] = sideband->GetParameter(3);
              ntBuffer[15] = sideband->GetParError(3);
              ntBuffer[16] = sideband->GetParameter(4);
              ntBuffer[17] = sideband->GetParError(4);
              ntBuffer[18] = sideband->GetParameter(5);
              ntBuffer[19] = sideband->GetParError(5);
              ntBuffer[20] = primaryContr;
              ntBuffer[21] = sidebandContr;
              ntBuffer[22] = sidebandNormDown[sbNormIter];
              ntBuffer[23] = sidebandNormUp[sbNormIter];
              ntBuffer[24] = chi2;
              ntBuffer[25] = pval;
              ntBuffer[26] = Chi2_pSigma0_250;
              ntBuffer[27] = (float) EffNumBins_pSigma0_250;
              ntBuffer[28] = Chi2_pSigma0_250 / double(EffNumBins_pSigma0_250);
              ntBuffer[29] = pvalpSigma0_250;
              ntBuffer[30] = nSigmapSigma0_250;
              ntBuffer[31] = nSigmapSigma0_200;
              ntBuffer[32] = nSigmapSigma0_150;
              ntBuffer[33] = (float) isCFneg;
              if (potential == 0) {
                ntBuffer[34] = d0;
                ntBuffer[35] = REf0inv;
                ntBuffer[36] = IMf0inv;
              } else if (potential == 1) {
                ntBuffer[34] = deltap0;
                ntBuffer[35] = deltap1;
                ntBuffer[36] = deltap2;
                ntBuffer[37] = etap0;
                ntBuffer[38] = etap1;
                ntBuffer[39] = etap2;
              }

              ntResult->Fill(ntBuffer);
              ++iterID;

              delete fitter;

              if (NumIter == 0) {
                std::cout << "Skipping all systematic variations \n";
                goto exitThroughTheGiftShop;
              }
            }
          }
          delete sideband;
          delete Ck_SideBand_draw;
          delete Ck_SideBand;
        }
      }
    }
  }
  exitThroughTheGiftShop:

  ntResult->Write();
  param->Close();

  delete side;
  delete CATSinput;
  delete ntResult;
  delete param;
  delete tidy;
  std::cout << "\n";
  return;
}

/// =====================================================================================
void FitSigma0(char *argv[]) {
  const unsigned& NumIter = atoi(argv[1]);
  TString InputDir = argv[2];
  TString SystDir = argv[3];
  TString trigger = argv[4];
  TString suffix = argv[5];
  TString OutputDir = argv[6];
  const int potential = atoi(argv[7]);
  std::vector<double> params;
  if (potential == 0) {
    if (!argv[8] || !argv[9] || !argv[10]) {
      std::cout << "ERROR: Missing the scattering parameters\n";
      return;
    }
    params.push_back(atof(argv[8]));  // d0
    params.push_back(atof(argv[9]));  // REf0inv
    params.push_back(atof(argv[10]));  // IMf0inv
  } else if (potential == 1) {
    if (!argv[8] || !argv[9] || !argv[10] || !argv[11] || !argv[12]
        || !argv[13]) {
      std::cout << "ERROR: Missing the parameters for delta/eta\n";
      return;
    }
    params.push_back(atof(argv[8]));   // deltap0
    params.push_back(atof(argv[9]));   // deltap1
    params.push_back(atof(argv[10]));   // deltap2
    params.push_back(atof(argv[11]));  // etap0
    params.push_back(atof(argv[12]));  // etap1
    params.push_back(atof(argv[13]));  // etap2
  }
  FitSigma0(NumIter, InputDir, SystDir, trigger, suffix, OutputDir, potential,
            params);
}
