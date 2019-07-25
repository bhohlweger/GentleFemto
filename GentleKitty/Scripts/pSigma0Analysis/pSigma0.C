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
const int nSidebandPars = 4;

/// =====================================================================================
/// Fit for the sidebands
auto sidebandFit =
    [ ] (double *x, double *p) {
      return p[0] + p[1] * std::exp(-0.5*std::pow(((x[0]-p[2])/p[3]), 2));
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
void FitSigma0(TString InputDir, TString SystInputDir, TString trigger,
               TString OutputDir, const int potential) {
  bool batchmode = false;
  bool debugPlots = true;
  double d0, REf0inv, IMf0inv, deltap0, deltap1, deltap2, etap0, etap1, etap2;

  DreamPlot::SetStyle();
  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  int nArguments = 24;
  TString varList =
      TString::Format(
          "IterID:femtoFitRange:ppRadius:bl_a:bl_b:purity:primaryContrib:fakeContrib:SBfitVal:"
          "sb_p0:sb_p0_err:sb_p1:sb_p1_err:sb_p2:sb_p2_err:sb_p3:sb_p3_err:"
          "chi2Local:ndf:chi2NDF:pval:nSigma250:nSigma200:nSigma150")
          .Data();
  TNtuple* ntResult = new TNtuple("fitResult", "fitResult", varList.Data());

  Float_t ntBuffer[nArguments];
  int iterID = 0;
  bool useBaseline = true;

  TString graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                          OutputDir.Data(), potential);
  auto param = new TFile(graphfilename, "RECREATE");

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  std::vector<TH1F> histSysVar;
  std::vector<TH1F> histSidebandSysVar;
  std::vector<TGraphAsymmErrors> grSysVar;
  std::vector<TGraphAsymmErrors> grSidebandSysVar;
  std::vector<double> puritySigma0;

  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  TString dataGrName = "Graph_from_hCk_ReweightedpSigma0_0MeV";
  const int nSidebandHist = 5;
  const int rebinSideband = 10;

  DreamSystematics protonsigma(DreamSystematics::pSigma0);
  for (int i = 0; i <= protonsigma.GetNumberOfVars(); ++i) {
    auto CATSinputVar = new CATSInputSigma0();
    auto appendixVar = TString::Format("%i", i);
    CATSinputVar->ReadSigma0CorrelationFile(SystInputDir.Data(), trigger.Data(),
                                            appendixVar.Data());
    CATSinputVar->CountPairs(SystInputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    puritySigma0.push_back(CATSinputVar->GetSigma0PurityPt());
    CATSinputVar->ObtainCFs(10, 250, 400);
    auto dataHistVar = CATSinputVar->GetCF("pSigma0", dataHistName.Data());
    auto dataGrVar = CATSinputVar->GetCFGr("pSigma0", dataGrName.Data());
    histSysVar.push_back(*dataHistVar);
    grSysVar.push_back(*dataGrVar);
    delete CATSinputVar;
    auto sideVar = new SidebandSigma();
    sideVar->SetRebin(rebinSideband);
    sideVar->SetSideBandFile(InputDir.Data(), trigger.Data(),
                             appendixVar.Data());
    sideVar->SetNormalizationRange(250, 400);
    sideVar->SideBandCFs();
    histSidebandSysVar.push_back(*(sideVar->GetSideBands(nSidebandHist)));
    grSidebandSysVar.push_back(*(sideVar->GetSideBandGraph(nSidebandHist)));
    delete sideVar;
  }

  TH1F *dataHist = &histSysVar[0];

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.
  const int binwidth = dataHist->GetBinWidth(1);
  int NumMomBins_pSigma = int(600 / binwidth);
  double kMin_pSigma = dataHist->GetBinCenter(1) - binwidth / 2.;
  double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;
  int NumMomBins_pSigma_draw = 50;
  double kMin_pSigma_draw = -4.99;
  double kMax_pSigma_draw = 495.01;
  const float drawBinWidth = float(kMax_pSigma_draw - kMin_pSigma_draw)
      / float(NumMomBins_pSigma_draw);

  std::cout << "kMin: " << kMin_pSigma << std::endl;
  std::cout << "kMax: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins: " << NumMomBins_pSigma << std::endl;

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

  // sideband fit range systematic variation
  const std::vector<double> sidebandFitRange = { { 650, 600, 700 } };

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the model, fitter, etc.
  DLM_Ck* Ck_pSigma0;
  DLM_Ck* Ck_pSigma0_draw;
  CATS AB_pSigma0;
  CATS AB_pSigma0_draw;
  if (potential == 2) {  // Lednicky coupled channel model fss2
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
  } else if (potential == 5) {  // Haidenbauer WF
    std::cout << "Running with the NSC97f potential \n";
    // NSC97f is valid up to 350 MeV, therefore we have to adopt
    double kMax_pSigma_NSC97f = kMin_pSigma;
    int NumMomBins_pSigma_NSC97f = 0;
    while (kMax_pSigma_NSC97f < 348 - binwidth) {
      kMax_pSigma_NSC97f += binwidth;
      ++NumMomBins_pSigma_NSC97f;
    }
    tidy->GetCatsProtonSigma0(&AB_pSigma0, NumMomBins_pSigma_NSC97f,
                              kMin_pSigma, kMax_pSigma_NSC97f,
                              TidyCats::sGaussian, TidyCats::pSigma0NSC97f);

    tidy->GetCatsProtonSigma0(&AB_pSigma0_draw, NumMomBins_pSigma_draw,
                              kMin_pSigma_draw, kMax_pSigma_draw,
                              TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
    AB_pSigma0.KillTheCat();
    AB_pSigma0_draw.KillTheCat();
    Ck_pSigma0 = new DLM_Ck(1, 0, AB_pSigma0);
    Ck_pSigma0_draw = new DLM_Ck(1, 0, AB_pSigma0_draw);
  } else if (potential == 6) {  // flat - sideband only
    std::cout << "Running with a flat correlation function = sideband only \n";
    Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins_pSigma, kMin_pSigma, kMax_pSigma,
                            Flat_Residual);
    Ck_pSigma0_draw = new DLM_Ck(1, 0, NumMomBins_pSigma_draw, kMin_pSigma_draw,
                                 kMax_pSigma_draw, Flat_Residual);

  }

  float counter = 0;
  float total = histSysVar.size() * femtoFitRegionUp.size() * 4
      * sidebandFitRange.size();

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Systematic variations

  // 1. Systematic variations of the data
  for (size_t systDataIter = 0; systDataIter < histSysVar.size();
      ++systDataIter) {
    TH1F* currentHist = &histSysVar[systDataIter];
    TGraphAsymmErrors *currentGr = &grSysVar[systDataIter];
    TH1F* currentSidebandHist = &histSidebandSysVar[systDataIter];
    TGraphAsymmErrors *currentSidebandGr = &grSidebandSysVar[systDataIter];
    kMin_pSigma = currentHist->GetBinCenter(1) - binwidth / 2.;

    /// Prefit for the baseline
    auto Prefit = (TH1F*) currentSidebandHist->Clone(
        Form("%i_%s_prefit", int(systDataIter),
             currentSidebandHist->GetName()));
    auto funct_0 = new TF1("myPol0", "pol0", 250, 600);
    Prefit->Fit(funct_0, "FSNRMQ");

    TF1* funct_1 = new TF1("myPol1", "pol1", 250, 600);
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
    static std::vector<CATSLambdaParam> lambdaParams;
    lambdaParams.resize(0);
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
        for (size_t sbNormIter = 0; sbNormIter < sidebandFitRange.size();
            ++sbNormIter) {

          auto sideband = new TF1(Form("sideband_%i", iterID), sidebandFit, 0,
                                  sidebandFitRange[sbNormIter], nSidebandPars);
          sideband->SetParameter(0, 1.);
          sideband->SetParameter(1, 0.4);
          sideband->SetParLimits(1, 0.4, 0.9);
          sideband->SetParameter(2, -70);
          sideband->SetParameter(3, 135);
          currentSidebandGr->Fit(sideband, "NRMQEX0");

          DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins_pSigma,
                                           kMin_pSigma, kMax_pSigma,
                                           sidebandFitCATS);

          DLM_Ck* Ck_SideBand_draw = new DLM_Ck(0, nSidebandPars,
                                                NumMomBins_pSigma_draw,
                                                kMin_pSigma_draw,
                                                kMax_pSigma_draw,
                                                sidebandFitCATS);

          for (unsigned i = 0; i < nSidebandPars; ++i) {
            Ck_SideBand->SetPotPar(i, sideband->GetParameter(i));
            Ck_SideBand_draw->SetPotPar(i, sideband->GetParameter(i));
          }
          Ck_SideBand->Update();
          Ck_SideBand_draw->Update();

          std::cout
              << "\r Processing progress: "
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
              double mom, dataY, dataErr, theoryX, theoryY;

              for (int iPoint = 0; iPoint <= currentGr->GetN(); ++iPoint) {
                currentGr->GetPoint(iPoint, mom, dataY);
                dataErr = currentGr->GetErrorY(iPoint);

                if (mom > femtoFitRegionUp[femtoFitIter]) {
                  continue;
                }

                theoryY = FitResult_pSigma0.Eval(mom);
                if (mom < 250) {
                  Chi2_pSigma0_250 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  ++EffNumBins_pSigma0_250;
                }
                if (mom < 200) {
                  Chi2_pSigma0_200 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  ++EffNumBins_pSigma0_200;
                }
                if (mom < 150) {
                  Chi2_pSigma0_150 += (dataY - theoryY) * (dataY - theoryY)
                      / (dataErr * dataErr);
                  ++EffNumBins_pSigma0_150;
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
              TGraph grCFSigmaGenuineSidebandExtrapolate;
              grCFSigmaGenuineSidebandExtrapolate.SetName(
                  Form("GenuineSideBand_%i", iterID));

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
                grCFSigmaGenuineSidebandExtrapolate.SetPoint(
                    i, mom, Ck_SideBand_draw->Eval(mom));
              }

              double ck;
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

              currentGr->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
              currentSidebandGr->SetTitle(
                  ";#it{k}* (MeV/#it{c}); C_{p#minus(#Lambda#gamma)}(#it{k}*)");

              if (iterID == 0) {
                /// beautification
                DreamPlot::SetStyleGraph(currentGr, 24, kBlue + 3);
                DreamPlot::SetStyleGraph(currentSidebandGr, 20, kRed + 2);
                currentGr->GetXaxis()->SetTitleSize(28);
                currentGr->GetYaxis()->SetTitleSize(28);
                currentSidebandGr->GetXaxis()->SetTitleSize(28);
                currentSidebandGr->GetYaxis()->SetTitleSize(28);
                sideband->SetLineWidth(2);
                sideband->SetLineStyle(2);
                sideband->SetLineColor(kGray + 1);
                FitResult_pSigma0.SetLineWidth(2);
                FitResult_pSigma0.SetLineColor(kRed + 2);
                grCFSigmaSideband.SetLineWidth(2);
                grCFSigmaSideband.SetLineStyle(2);
                grCFSigmaSideband.SetLineColor(kGray + 1);

                if (debugPlots) {
                  auto f = new TCanvas("baseline");
                  currentSidebandGr->Draw("APEZ");
                  currentSidebandGr->GetXaxis()->SetRangeUser(0, 600);
                  currentSidebandGr->GetYaxis()->SetTitle("C(#it{k}*)");
                  currentGr->Draw("pezsame");
                  currentGr->GetXaxis()->SetRangeUser(0, 600);
                  const int nBins = prefit_a.size();
                  TF1* baselines[nBins];
                  for (int i = 0; i < nBins; ++i) {
                    baselines[i] = new TF1(Form("baseline_%i", i), "pol1", 0,
                                           900);
                    baselines[i]->SetParameter(0, prefit_a[i]);
                    baselines[i]->SetParameter(1, prefit_b[i]);
                    baselines[i]->SetLineColor(kGreen + 2);
                    baselines[i]->Draw("same");
                  }
                  auto leg2 = new TLegend(0.4, 0.68, 0.6, 0.85);
                  leg2->SetTextFont(42);
                  leg2->SetTextSize(0.05);
                  leg2->AddEntry(currentGr, "p#minus#Sigma^{0}", "pe");
                  leg2->AddEntry(currentSidebandGr,
                                 "p#minus(#Lambda#gamma) sideband (unscaled)",
                                 "pe");
                  leg2->AddEntry(baselines[0], "Baseline fits", "l");
                  leg2->Draw("same");
                  f->Print(Form("%s/CF_baseline.pdf", OutputDir.Data()));
                  delete f;

                    auto g = new TCanvas();
                  DreamPlot::SetStyleGraph (grPrefitContour);
                  grPrefitContour->SetTitle("; #it{a}; #it{b}");
                  grPrefitContour->SetLineStyle(2);
                  auto grDots = new TGraph();
                  DreamPlot::SetStyleGraph(grDots, 20, kRed + 2);
                  for (size_t i = 0; i < prefit_a.size(); ++i) {
                    grDots->SetPoint(i, prefit_a[i],
                                     prefit_b[i]);
                  }
                  grPrefitContour->Draw("AL");
                  grDots->Draw("PEsame");
                  g->Print(
                      Form("%s/Prefit_%i.pdf", OutputDir.Data(), potential));
                  delete grDots;
                  delete g;
                }

                currentGr->Write();
                grCFSigmaRaw.Write();
                grCFSigmaMain.Write();
                grCFSigmaFeed.Write();
                grPrefitContour->Write("fitContour");
                sideband->Write(Form("SidebandFitNotScaled_%i", iterID));
                currentSidebandGr->SetName(Form("SidebandMerged_%i", iterID));
                currentSidebandGr->Write();
              }

              grCFSigmaExtrapolate.Write();
              grCFSigmaSidebandExtrapolate.Write();
              grCFSigmaGenuineSidebandExtrapolate.Write();
              grCFSigmaSideband.Write();
              FitResult_pSigma0.Write(Form("Fit_%i", iterID));
              currentGr->SetName(
                  TString::Format("HistCF_Var_%i", int(systDataIter)));
              currentGr->Write();

              param->cd();
              ntBuffer[0] = iterID;
              ntBuffer[1] = femtoFitRegionUp[femtoFitIter];
              ntBuffer[2] = sourceSize[sizeIter];
              ntBuffer[3] = bl_a;
              ntBuffer[4] = bl_b;
              ntBuffer[5] = sigmaPurity;
              ntBuffer[6] = primaryContr;
              ntBuffer[7] = sidebandContr;
              ntBuffer[8] = sidebandFitRange[sbNormIter];
              ntBuffer[9] = sideband->GetParameter(0);
              ntBuffer[10] = sideband->GetParError(0);
              ntBuffer[11] = sideband->GetParameter(1);
              ntBuffer[12] = sideband->GetParError(1);
              ntBuffer[13] = sideband->GetParameter(2);
              ntBuffer[14] = sideband->GetParError(2);
              ntBuffer[15] = sideband->GetParameter(3);
              ntBuffer[16] = sideband->GetParError(3);
              ntBuffer[17] = Chi2_pSigma0_250;
              ntBuffer[18] = (float) EffNumBins_pSigma0_250;
              ntBuffer[19] = Chi2_pSigma0_250 / double(EffNumBins_pSigma0_250);
              ntBuffer[20] = pvalpSigma0_250;
              ntBuffer[21] = nSigmapSigma0_250;
              ntBuffer[22] = nSigmapSigma0_200;
              ntBuffer[23] = nSigmapSigma0_150;
              ntResult->Fill(ntBuffer);
              ++iterID;

              delete fitter;
            }
          }
          delete sideband;
          delete Ck_SideBand_draw;
          delete Ck_SideBand;
        }
      }
    }
  }

  ntResult->Write();
  param->Close();

  delete CATSinput;
  delete ntResult;
  delete param;
  delete tidy;
  std::cout << "\n";
  return;
}

/// =====================================================================================
void FitSigma0(char *argv[]) {
  TString InputDir = argv[1];
  TString SystDir = argv[2];
  TString trigger = argv[3];
  TString OutputDir = argv[4];
  const int potential = atoi(argv[5]);
  FitSigma0(InputDir, SystDir, trigger, OutputDir, potential);
}
