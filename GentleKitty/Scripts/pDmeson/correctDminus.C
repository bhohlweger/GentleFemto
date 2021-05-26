#include "DmesonTools.h"

/// =====================================================================================
void correctDminus(TString InputDir, TString trigger, int errorVar) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const int nBoot = 2500;

  TString OutputDir = InputDir;
  OutputDir += "/fit/";

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

  DreamPlot::SetStyle();

  int nArguments = 25;
  TString varList =
      TString::Format(
          "BootID:systID:ppRadius:primaryContrib:flatContrib:dstarContrib:beautyContrib:sidebandContrib:chi2SidebandLeft:chi2SidebandRight:"
          "ndfBkg:chi2bkg:nSigmaBkf:ndfSig:chi2Coulomb:nSigmaCoulomb:chi2Haidenbauer:nSigmaHaidenbauer:chi2Model1:nSigmaModel1:chi2Model3:nSigmaModel3:chi2Model4:nSigmaModel4")
          .Data();
  auto ntResult = new TNtuple("fitResult", "fitResult", varList.Data());
  auto tupleSideband = new TNtuple("sideband", "sideband", "kstar:cf:BootID");
  auto tupleTotalFit = new TNtuple("totalFit", "totalFit", "kstar:cf:BootID");
  auto tupleCorrected = new TNtuple("correctedCF", "correctedCF",
                                    "kstar:cf:BootID");
  auto tupleSidebandLeft = new TNtuple("sidebandLeft", "sidebandLeft",
                                       "kstar:cf:BootID");
  auto tupleSidebandRight = new TNtuple("sidebandRight", "sidebandRight",
                                        "kstar:cf:BootID");
  auto tupleCoulomb = new TNtuple("coulomb", "coulomb", "kstar:cf:BootID");
  auto tupleHaidenbauer = new TNtuple("haidenbauer", "haidenbauer",
                                      "kstar:cf:BootID");
  auto tupleModel1 = new TNtuple("model1", "model1", "kstar:cf:BootID");
  auto tupleModel3 = new TNtuple("model3", "model3", "kstar:cf:BootID");
  auto tupleModel4 = new TNtuple("model4", "model4", "kstar:cf:BootID");

  float ntBuffer[nArguments];
  bool useBaseline = true;

  TString graphfilename = TString::Format("%s/Fit_pDminus_corrected_%s.root",
                                          OutputDir.Data(), errorName.Data());
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

  auto meDistForFeedDown = GetMEDist(InputFileName, "HM_CharmFemto_", "0",
                                     dataGrName, normLower, normUpper, 1);

  std::vector<double> startParams;

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  std::vector<double> sidebandFitRange = { { 1500, 1200, 1800 } };

  const int nBins = 95;
  const double kmin = 0;
  const double kmax = 950;
  double binWidth = (kmax - kmin) / double(nBins);

  const int nBinsModel = 80;
  const double kminModel = 10;
  const double kmaxModel = 410;
  double binWidthModel = (kmaxModel - kminModel) / double(nBinsModel);

  /// Femtoscopic radius systematic variations
  const double rErr = 0.08;
  const double rDefault = 0.89;
  const double rScale = 0.84;
  std::vector<double> sourceSize;
  if (errorVar == 1) {
    sourceSize = { {rDefault}};
  } else {
    sourceSize = { {rDefault, rScale * rDefault - rErr, rDefault + rErr}};
  }

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
    auto file = TFile::Open(
        Form("%s/Fractions_masswindow_2sigma_sphericity_0_1_%s.root",
             InputDir.Data(), fileAppendix.Data()));
    auto histFracStat = (TH1F*) file->Get("hFractions_DmPr");
    auto grFracSyst = (TGraphAsymmErrors*) file->Get("gFractionsSyst_DmPr");
    double purity = histFracStat->GetBinContent(1);
    double purityStatErr = histFracStat->GetBinError(1);
    double puritySystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(0);
    DmesonPurity.push_back(purity);
    DmesonPurityErr.push_back(
        std::sqrt(
            purityStatErr * purityStatErr + puritySystErr * puritySystErr));

    double bfeed = histFracStat->GetBinContent(2);
    double bfeedStatErr = histFracStat->GetBinError(2);
    double bfeedSystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(1);
    Bfeeddown.push_back(bfeed);
    BfeeddownErr.push_back(
        std::sqrt(bfeedStatErr * bfeedStatErr + bfeedSystErr * bfeedSystErr));

    double dstarfeed = histFracStat->GetBinContent(3);
    double dstarfeedStatErr = histFracStat->GetBinError(3);
    double dstarfeedSystErr = (errorVar == 1) ? 0 : grFracSyst->GetErrorY(2);
    DstarFeeding.push_back(dstarfeed);
    DstarFeedingErr.push_back(
        std::sqrt(
            dstarfeedStatErr * dstarfeedStatErr
                + dstarfeedSystErr * dstarfeedSystErr));
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the model, fitter, etc.

  DLM_Ck *DLM_CFflat = new DLM_Ck(1, 0, nBins, kmin, kmax, Flat_Residual);
  DLM_CFflat->Update();

  auto tidyCats = new TidyCats();
  CATS catsDstar, catsDminusCoulomb, catsDminusHaidenbauer;

  // Coulomb
  tidyCats->GetCatsProtonDminus(&catsDminusCoulomb, nBinsModel, kminModel,
                                kmaxModel, TidyCats::pDCoulombOnly,
                                TidyCats::sGaussian);
  catsDminusCoulomb.SetAnaSource(0, rDefault);
  catsDminusCoulomb.KillTheCat();
  DLM_Ck* DLM_Coulomb = new DLM_Ck(1, 0, catsDminusCoulomb);

  // Haidenbauer
  tidyCats->GetCatsProtonDminus(&catsDminusHaidenbauer, nBinsModel, kminModel,
                                kmaxModel, TidyCats::pDminusHaidenbauer,
                                TidyCats::sGaussian);
  catsDminusHaidenbauer.SetAnaSource(0, rDefault);
  catsDminusHaidenbauer.KillTheCat();
  //DLM_Ck* DLM_Haidenbauer = new DLM_Ck(1, 0, catsDminusHaidenbauer);

  // now store all of them for further computation
  std::vector<DLM_CkDecomposition*> coulombModels, haidenbauerModels,
      model1Models, model3Models;
  std::vector<TGraphErrors*> model4Models;
  TString HomeDir = gSystem->GetHomeDirectory().c_str();

  for (const auto &sourceRad : sourceSize) {
    DLM_Coulomb->SetSourcePar(0, sourceRad);
    DLM_Coulomb->Update();
    DLM_CkDecomposition *CkDec_CFCoulomb = new DLM_CkDecomposition(
        Form("pDminusCoulomb_%.3f", sourceRad), 0, *DLM_Coulomb,
        momentumResolution);
    CkDec_CFCoulomb->Update();
    coulombModels.push_back(CkDec_CFCoulomb);

    auto grHaidenbauer = new TGraph(
        TString::Format(
            "%s/CERNHome/D-mesons/Analysis/Models/Haidenbauer_%.2f_fm.dat",
            HomeDir.Data(), sourceRad));
    auto DLM_Haidenbauer = getDLMCk(grHaidenbauer);

    DLM_CkDecomposition *CkDec_CFHaidenbauer = new DLM_CkDecomposition(
        Form("pDminusHaidenbauer_%f", sourceRad), 0, *DLM_Haidenbauer,
        momentumResolution);
    CkDec_CFHaidenbauer->Update();
    haidenbauerModels.push_back(CkDec_CFHaidenbauer);
  }

  // TODO For now we have the Yuki CF only for 0.9 fm! To be updated
  const double sourceRadYuki = 0.9;
  auto grYukiModel1 = getCkFromYuki(3, sourceRadYuki);  // model 1 = potential 3
  auto DLM_YukiModel1 = getDLMCk(grYukiModel1);
  DLM_CkDecomposition *CkDec_Model1 = new DLM_CkDecomposition(
      Form("pDminusModel1_%f", sourceRadYuki), 0, *DLM_YukiModel1,
      momentumResolution);
  model1Models.push_back(CkDec_Model1);

  auto grYukiModel3 = getCkFromYuki(4, sourceRadYuki);  // model 3 = potential 4
  auto DLM_YukiModel3 = getDLMCk(grYukiModel3);
  DLM_CkDecomposition *CkDec_Model3 = new DLM_CkDecomposition(
      Form("pDminusModel2_%f", sourceRadYuki), 0, *DLM_YukiModel3,
      momentumResolution);
  model3Models.push_back(CkDec_Model3);

  auto grYukiModel4 = getCkFromYuki(5, sourceRadYuki);  // model 4 = potential 5 - here things need to be done a bit differently as we have uncertainties on the scat. params that need to be bootstrapped
  model4Models.push_back(grYukiModel4);

  tidyCats->GetCatsProtonDstarminus(&catsDstar, nBins, kmin, kmax,
                                    TidyCats::pDCoulombOnly,
                                    TidyCats::sGaussian);
  catsDstar.SetAnaSource(0, rDefault);
  catsDstar.KillTheCat();

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

    const int femtoIt =
        (iBoot == 0 || errorVar == 1) ?
            0 : gRandom->Uniform() * sourceSize.size();
    const double femtoRad = sourceSize.at(femtoIt);

    // since purity etc. depends on the cut variation, this needs to be done consistently
    int nSystVar =
        (iBoot == 0 || errorVar == 1) ? 0 : gRandom->Uniform() * grCFvec.size();
    auto grVarCF = grCFvec.at(nSystVar);
    auto grVarCFSidebandLeft = grSBLeftvec.at(nSystVar);
    auto grVarCFSidebandRight = grSBRightvec.at(nSystVar);

    double dMesonPur = DmesonPurity.at(nSystVar);
    double bFeed = Bfeeddown.at(nSystVar);
    double dStarFeed = DstarFeeding.at(nSystVar);

    // now let's sample the uncertainties
    if (iBoot != 0) {
      dMesonPur += DmesonPurityErr.at(nSystVar)
          * std::round(gRandom->Uniform(-1.4999, 1.4999));
      bFeed += BfeeddownErr.at(nSystVar)
          * std::round(gRandom->Uniform(-1.4999, 1.4999));
      dStarFeed += DstarFeedingErr.at(nSystVar)
          * std::round(gRandom->Uniform(-1.4999, 1.4999));
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
        (iBoot == 0 || errorVar == 1) ?
            sidebandFitRange.at(0) : getBootstrapFromVec(sidebandFitRange);

    auto fitSidebandLeft = new TF1("fitSidebandLeft", "pol3", 0, sidebandRange);
    getImprovedStartParamsPol3(grCFBootstrapSidebandLeft, sidebandRange,
                               startParams);
    fitSidebandLeft->SetParameters(&startParams[0]);

    int workedLeft = grCFBootstrapSidebandLeft->Fit(fitSidebandLeft, "RQ", "",
                                                    200, sidebandRange);
    const double chiSqSidebandLeft = fitSidebandLeft->GetChisquare()
        / double(fitSidebandLeft->GetNDF());

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandLeft,
                                                          0.683);

    auto fitSidebandRight = new TF1("fitSidebandRight", "pol3", 0,
                                    sidebandRange);
    getImprovedStartParamsPol3(grCFBootstrapSidebandRight, sidebandRange,
                               startParams);
    fitSidebandRight->SetParameters(&startParams[0]);

    int workedRight = grCFBootstrapSidebandRight->Fit(fitSidebandRight, "RQ",
                                                      "", 200, sidebandRange);
    const double chiSqSidebandRight = fitSidebandRight->GetChisquare()
        / double(fitSidebandRight->GetNDF());

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grSidebandRight,
                                                          0.683);

    // stop here when the fit fails!
    if (workedLeft != 0 || workedRight != 0 || chiSqSidebandLeft > 50
        || chiSqSidebandRight > 50) {
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

    // Set the radii

    DLM_pDstar->SetSourcePar(0, femtoRad);
    DLM_pDstar->Update();
    auto DLM_sideband = getDLMCk(totalSideband);

    DLM_CkDecomposition CkDec_CFflat("pDminusTotal", 3, *DLM_CFflat, nullptr);
    CkDec_CFflat.AddPhaseSpace(0, &meDistForFeedDown);
    DLM_CkDecomposition CkDec_pDstar("pDstarminus", 0, *DLM_pDstar,
                                     momentumResolution);
    CkDec_pDstar.AddPhaseSpace(&meDistForFeedDown);
    DLM_CkDecomposition CkDec_sideband("sideband", 0, *DLM_sideband, nullptr);

    CkDec_CFflat.AddContribution(0, pDstarContrib,
                                 DLM_CkDecomposition::cFeedDown, &CkDec_pDstar,
                                 decayKindematicsDstar);
    CkDec_CFflat.AddContribution(1, flatContrib, DLM_CkDecomposition::cFake);
    CkDec_CFflat.AddContribution(2, sidebandContrib,
                                 DLM_CkDecomposition::cFeedDown,
                                 &CkDec_sideband);
    CkDec_CFflat.Update();

    TGraph *grFitTotalFine = new TGraph();
    TGraph *grCFRaw = new TGraph();

    auto grCFBootstrap =
        (iBoot == 0) ? &grCFvec.at(0) : getBootstrapGraph(&grVarCF);

    count = 0;
    for (int i = kmin; i < kmax;) {
      double mom = i;
      grCFRaw->SetPoint(count, mom, DLM_CFflat->Eval(mom));
      grFitTotalFine->SetPoint(count, mom, CkDec_CFflat.EvalCk(mom));
      tupleTotalFit->Fill(mom, CkDec_CFflat.EvalCk(mom), iBoot);
      tupleSideband->Fill(mom, totalSideband->Eval(mom), iBoot);
      tupleSidebandLeft->Fill(mom, fitSidebandLeft->Eval(mom), iBoot);
      tupleSidebandRight->Fill(mom, fitSidebandRight->Eval(mom), iBoot);
      i += binWidth;
      count++;
    }

    // Correct the correlation function and do the model comparison
    auto coulombModelIt = coulombModels.at(femtoIt);
    auto haidenbauerModelIt = haidenbauerModels.at(femtoIt);
    auto model1ModelIt = model1Models.at(0);  // TODO - FIX THE RADIUS
    auto model3ModelIt = model3Models.at(0);  // TODO - FIX THE RADIUS
    auto model4GraphIt = model4Models.at(0);  // TODO - FIX THE RADIUS
    auto DLM_YukiModel4 = getDLMCk(getBootstrapGraph(model4GraphIt));
    auto model4ModelIt = new DLM_CkDecomposition("pDminusModel4", 0,
                                                 *DLM_YukiModel4,
                                                 momentumResolution);  // TODO eventually remove the momRes

    for (int i = kminModel; i < kmaxModel; i += binWidthModel) {
      double mom = i;
      tupleCoulomb->Fill(mom, coulombModelIt->EvalCk(mom), iBoot);
      tupleHaidenbauer->Fill(mom, haidenbauerModelIt->EvalCk(mom), iBoot);
      tupleModel1->Fill(mom, model1ModelIt->EvalCk(mom), iBoot);
      tupleModel3->Fill(mom, model3ModelIt->EvalCk(mom), iBoot);
      tupleModel4->Fill(mom, model4ModelIt->EvalCk(mom), iBoot);
    }

    TGraph *grFitTotal = new TGraph();

    double EffNumBinsFlat = 0;
    double Chi2Flat = 0;
    double EffNumBins = 0;
    double Chi2Coulomb = 0;
    double Chi2Haidenbauer = 0;
    double Chi2Model1 = 0;
    double Chi2Model3 = 0;
    double Chi2Model4 = 0;

    static double mom, dataY, dataBootstrapY, dataErr, corrBootstrapY,
        corrDataY, corrErr, flatY, coulombY, haidenbauerY, model1Y, model3Y,
        model4Y, sbmom, sbY, sbErr;

    for (int iPoint = 0; iPoint <= grVarCF.GetN(); ++iPoint) {
      grCFBootstrap->GetPoint(iPoint, mom, dataBootstrapY);
      grVarCF.GetPoint(iPoint, mom, dataY);
      dataErr = grVarCF.GetErrorY(iPoint);
      if (mom > 900) {
        break;
      }

      for (int isb = 0; isb <= totalSideband->GetN(); ++isb) {
        totalSideband->GetPoint(isb, sbmom, sbY);
        if (std::abs(sbmom - mom) < 2) {
          sbErr = totalSideband->GetErrorY(isb);
          break;
        }
      }

      corrBootstrapY = 1.f
          + 1.f / primaryContrib * (dataBootstrapY - CkDec_CFflat.EvalCk(mom));

      tupleCorrected->Fill(mom, corrBootstrapY, iBoot);

      flatY = CkDec_CFflat.EvalCk(mom);
      grFitTotal->SetPoint(iPoint, mom, flatY);

      if (mom < 200) {  // evaluate the chi2 of the models
        ++EffNumBins;

        coulombY = coulombModelIt->EvalCk(mom);
        haidenbauerY = haidenbauerModelIt->EvalCk(mom);
        model1Y = model1ModelIt->EvalCk(mom);
        model3Y = model3ModelIt->EvalCk(mom);
        model4Y = model4ModelIt->EvalCk(mom);

        corrDataY = 1.f
            + 1.f / primaryContrib * (dataY - CkDec_CFflat.EvalCk(mom));
        corrErr = std::sqrt(
            dataErr * dataErr / (primaryContrib * primaryContrib)
                + sbErr * sbErr * sidebandContrib * sidebandContrib
                    / (primaryContrib * primaryContrib));
        Chi2Coulomb += (corrDataY - coulombY) * (corrDataY - coulombY)
            / (corrErr * corrErr);

        Chi2Haidenbauer += (corrDataY - haidenbauerY)
            * (corrDataY - haidenbauerY) / (corrErr * corrErr);

        Chi2Model1 += (corrDataY - model1Y) * (corrDataY - model1Y)
            / (corrErr * corrErr);

        Chi2Model3 += (corrDataY - model3Y) * (corrDataY - model3Y)
            / (corrErr * corrErr);

        Chi2Model4 += (corrDataY - model4Y) * (corrDataY - model4Y)
            / (corrErr * corrErr);

      } else {  // evaluate the chi2 of the background description
        Chi2Flat += (dataY - flatY) * (dataY - flatY) / (dataErr * dataErr);
        ++EffNumBinsFlat;

      }
    }

    double pvalFlat = TMath::Prob(Chi2Flat, round(EffNumBinsFlat));
    double nSigmaFlat = TMath::Sqrt(2) * TMath::ErfcInverse(pvalFlat);
    double pvalCoulomb = TMath::Prob(Chi2Coulomb, round(EffNumBins));
    double nSigmaCoulomb = TMath::Sqrt(2) * TMath::ErfcInverse(pvalCoulomb);
    double pvalHaidenbauer = TMath::Prob(Chi2Haidenbauer, round(EffNumBins));
    double nSigmaHaidenbauer = TMath::Sqrt(2)
        * TMath::ErfcInverse(pvalHaidenbauer);
    double pvalModel1 = TMath::Prob(Chi2Model1, round(EffNumBins));
    double nSigmaModel1 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalModel1);
    double pvalModel3 = TMath::Prob(Chi2Model3, round(EffNumBins));
    double nSigmaModel3 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalModel3);
    double pvalModel4 = TMath::Prob(Chi2Model4, round(EffNumBins));
    double nSigmaModel4 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalModel4);

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
    ntBuffer[1] = nSystVar;
    ntBuffer[2] = femtoRad;
    ntBuffer[3] = primaryContrib;
    ntBuffer[4] = flatContrib;
    ntBuffer[5] = pDstarContrib;
    ntBuffer[6] = bFeed;
    ntBuffer[7] = sidebandContrib;
    ntBuffer[8] = chiSqSidebandLeft;
    ntBuffer[9] = chiSqSidebandRight;
    ntBuffer[10] = (float) EffNumBinsFlat;
    ntBuffer[11] = Chi2Flat;
    ntBuffer[12] = nSigmaFlat;
    ntBuffer[13] = (float) EffNumBins;
    ntBuffer[14] = Chi2Coulomb;
    ntBuffer[15] = nSigmaCoulomb;
    ntBuffer[16] = Chi2Haidenbauer;
    ntBuffer[17] = nSigmaHaidenbauer;
    ntBuffer[18] = Chi2Model1;
    ntBuffer[19] = nSigmaModel1;
    ntBuffer[20] = Chi2Model3;
    ntBuffer[21] = nSigmaModel3;
    ntBuffer[22] = Chi2Model4;
    ntBuffer[23] = nSigmaModel4;
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

    delete DLM_YukiModel4;
    delete model4ModelIt;

    ++iBoot;
  }

/// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/// Exit through the gift shop...

  param->cd();
  ntResult->Write();
  tupleTotalFit->Write();
  tupleSideband->Write();
  tupleCorrected->Write();
  tupleCoulomb->Write();
  tupleHaidenbauer->Write();
  tupleModel1->Write();
  tupleModel3->Write();
  tupleModel4->Write();

  grCFvec.at(0).Write("dataDefault");
  grSBLeftvec.at(0).Write("sidebandLeftDefault");
  grSBRightvec.at(0).Write("sidebandRightDefault");
  meDistForFeedDown.Write("meDistDefault");

  auto list = new TList();
  auto fitFull = EvalBootstrap(tupleTotalFit, list, OutputDir, "background",
                               kmin, kmax, binWidth);
  auto dataCorrected = EvalBootstrap(tupleCorrected, &grCFvec.at(0), list,
                                     OutputDir, "corrected");
  auto sidebandFull = EvalBootstrap(tupleSideband, nullptr, OutputDir,
                                    "sidebandFull", kmin, kmax, binWidth);
  auto sidebandLeft = EvalBootstrap(tupleSidebandLeft, nullptr, OutputDir,
                                    "sidebandLeft", kmin, kmax, binWidth);
  auto sidebandRight = EvalBootstrap(tupleSidebandRight, nullptr, OutputDir,
                                     "sidebandRight", kmin, kmax, binWidth);

  auto coulomb = EvalBootstrap(tupleCoulomb, list, OutputDir, "coulomb",
                               kminModel, kmaxModel, binWidthModel, false);  // use full width as we have only the three radii to vary
  auto haidenbauer = EvalBootstrap(tupleHaidenbauer, list, OutputDir,
                                   "haidenbauer", kminModel, kmaxModel,
                                   binWidthModel, false);  // use full width as we have only the three radii to vary
  auto model1 = EvalBootstrap(tupleModel1, list, OutputDir, "model1", kminModel,
                              kmaxModel, binWidthModel, false);  // use full width as we have only the three radii to vary
  auto model3 = EvalBootstrap(tupleModel3, list, OutputDir, "model3", kminModel,
                              kmaxModel, binWidthModel, false);  // use full width as we have only the three radii to vary
  auto model4 = EvalBootstrap(tupleModel4, list, OutputDir, "model4", kminModel,
                              kmaxModel, binWidthModel);  // here also the scat. params are varied so we use the RMS

  fitFull->Write("fitFull");
  sidebandFull->Write("sidebandFull");
  sidebandLeft->Write("sidebandLeft");
  sidebandRight->Write("sidebandRight");

  dataCorrected->Write("correctedData");

  coulomb->Write("Ck_coulomb");
  haidenbauer->Write("Ck_haidenbauer");
  model1->Write("Ck_model1");
  model3->Write("Ck_model3");
  model4->Write("Ck_model4");

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
  int errorVar = atoi(argv[3]);

  correctDminus(InputDir, trigger, errorVar);
}
