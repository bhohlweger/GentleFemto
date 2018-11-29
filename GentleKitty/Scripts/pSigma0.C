#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
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
#include "SideBandFit.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"
#include "TNtuple.h"

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
void FitSigma0(const unsigned& NumIter, TString InputDir, TString appendix,
               TString ppFile, TString OutputDir) {
  DreamPlot::SetStyle();
  bool fastPlot = (NumIter == 0) ? true : false;
  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  TNtuple* ntResult =
      new TNtuple(
          "fitResult",
          "fitResult",
          "IterID:FitReg_pSigma:BLSlope:ppRadius:"
          "bl_a:bl_a_err:bl_b:bl_b_err:chi2NDFGlobal:pvalGlobal:chi2Local:ndf:chi2NDF:pval:nSigma");

  Float_t ntBuffer[13];
  int ppRadius;
  int iterID = 0;
  bool useBaselineSlope = true;

  TString graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                          OutputDir.Data(), NumIter);
  auto param = new TFile(graphfilename, "RECREATE");

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
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  auto dataHist = CATSinput->GetCF("pSigma0", dataHistName.Data());
  if (!dataHist) {
    std::cerr << "ERROR pSigma0 fitter: p-Sigma0 histogram missing\n";
    return;
  }
  auto sidebandHistUp = CATSinput->GetCF("pSigmaSBUp",
                                         "hCk_ReweightedpSigmaSBUpMeV_0");

  auto sidebandHistLow = CATSinput->GetCF("pSigmaSBLow",
                                          "hCk_ReweightedpSigmaSBLowMeV_0");
  if (!sidebandHistLow || !sidebandHistUp) {
    std::cerr << "ERROR pSigma0 fitter: p-pSigmaSB histogram missing\n";
    return;
  }

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.
  const int binwidth = 40;
  const unsigned NumMomBins_pSigma = 20;
  double kMin_pSigma = dataHist->GetBinCenter(1)
      - dataHist->GetBinWidth(1) / 2.f;
  const double kMax_pSigma = kMin_pSigma + binwidth * NumMomBins_pSigma;

  std::cout << "kMin_pSigma: " << kMin_pSigma << std::endl;
  std::cout << "kMax_pSigma: " << kMax_pSigma << std::endl;
  std::cout << "Binwidth: " << binwidth << std::endl;
  std::cout << "NumMomBins_pSigma: " << NumMomBins_pSigma << std::endl;

  // TODO TO BE DONE - Get the radius from the pp file
  double sourceSize[3];
  sourceSize[0] = 1.36;
  sourceSize[1] = 1.36 * 0.8;
  sourceSize[2] = 1.36 * 1.2;

  double FemtoRegion_pSigma[3][2];
  FemtoRegion_pSigma[0][0] = kMin_pSigma;
  FemtoRegion_pSigma[0][1] = 500;
  FemtoRegion_pSigma[1][0] = kMin_pSigma;
  FemtoRegion_pSigma[1][1] = 550;
  FemtoRegion_pSigma[2][0] = kMin_pSigma;
  FemtoRegion_pSigma[2][1] = 600;

  /// Lambda parameters
  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;

  const double sigmaPurity = 0.199;
  const double sigmaPrimary = 1.;

  const Particle proton(
      protonPurity,
      protonPrimary,
      { { protonLambda, (1. - protonPrimary)
          * (1. - protonLambda / (1. - protonPrimary)) } });

  const Particle sigma0(0.2, sigmaPrimary, { { 0 } });

  const CATSLambdaParam lambdaParams(proton, sigma0);
  lambdaParams.PrintLambdaParams();

  const float primFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Primary);
  const float fakeFrac = lambdaParams.GetLambdaParam(CATSLambdaParam::Fake);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Fit the sideband
  SideBandFit* side = new SideBandFit();
  auto SBmerge = side->AddCF(sidebandHistUp, sidebandHistLow, "SBmerge");
  auto sideband = new TF1("sideband", sidebandFit, kMin_pSigma, 650,
                          nSidebandPars);
  SBmerge->Fit(sideband, "FSNRMQ");

  /// Prefit for the baseline
  float p_a_prefit[4];
  float p_b_prefit[4];
  auto Prefit = (TH1F*) dataHist->Clone(Form("%s_prefit", dataHist->GetName()));
  auto funct_0 = new TF1("myPol0", "pol0", 250, 750);
  Prefit->Fit(funct_0, "FSNRMQ");
  p_a_prefit[0] = funct_0->GetParameter(0);
  p_b_prefit[0] = 0;

  TF1* funct_1 = new TF1("myPol1", "pol1", 250, 750);
  Prefit->Fit(funct_1, "FSNRMQ");
  gMinuit->SetErrorDef(4);  //note 4 and not 2!
  p_a_prefit[1] = funct_1->GetParameter(0);
  p_b_prefit[1] = funct_1->GetParameter(1);

  auto grPrefitContour = (TGraph*) gMinuit->Contour(40, 0, 1);

  p_a_prefit[2] = TMath::MinElement(grPrefitContour->GetN(),
                                    grPrefitContour->GetX());
  p_b_prefit[2] = grPrefitContour->Eval(p_a_prefit[2]);
  p_a_prefit[3] = TMath::MaxElement(grPrefitContour->GetN(),
                                    grPrefitContour->GetX());
  p_b_prefit[3] = grPrefitContour->Eval(p_a_prefit[3]);

  std::cout << "Result of the prefit to constrain the baseline\n";
  for (int i = 0; i < 4; ++i) {
    std::cout << i << " a: " << p_a_prefit[i] << " b: " << p_b_prefit[i]
              << "\n";
  }
  delete Prefit;

  std::cout << "\n\nStarting the systematic variations\n\n";

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// "Systematic" variations
  for (unsigned int FitReg_pSigma = 0; FitReg_pSigma < 3; ++FitReg_pSigma) {
    for (unsigned int BaselineSlope = 0; BaselineSlope < 4; ++BaselineSlope) {
      for (unsigned int SourceSizeIter = 0; SourceSizeIter < 3;
          ++SourceSizeIter) {

        if (BaselineSlope == 0) {
          useBaselineSlope = false;  //use baseline
        } else {
          useBaselineSlope = true;  // no baseline
        }

        std::cout << "\n\n_________________________\n";
        std::cout << "Processing iteration " << iterID << "\n";
        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Correlation function
        DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins_pSigma, kMin_pSigma,
                                        kMax_pSigma, Lednicky_gauss_Sigma0);
        Ck_pSigma0->SetSourcePar(0, sourceSize[SourceSizeIter]);
        Ck_pSigma0->Update();

        // TODO TO BE FIXED - MOMENTUM RESOLUTION MATRIX!
        DLM_CkDecomposition CkDec_pSigma0("pSigma0", 1, *Ck_pSigma0,
                                          CATSinput->GetSigmaFile(0));

        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Sidebands
        DLM_Ck* Ck_SideBand = new DLM_Ck(0, nSidebandPars, NumMomBins_pSigma,
                                         kMin_pSigma, kMax_pSigma,
                                         sidebandFitCATS);
        std::cout << "Result of the sideband fit\n";
        for (unsigned i = 0; i < sideband->GetNumberFreeParameters(); ++i) {
          std::cout << i << " " << sideband->GetParameter(i) << std::endl;
          Ck_SideBand->SetPotPar(i, sideband->GetParameter(i));
        }
        Ck_SideBand->Update();

        DLM_CkDecomposition CkDec_SideBand("pSigma0SideBand", 0, *Ck_SideBand,
                                           nullptr);
        CkDec_pSigma0.AddContribution(
            0, lambdaParams.GetLambdaParam(CATSLambdaParam::Fake),
            DLM_CkDecomposition::cFake, &CkDec_SideBand);
        CkDec_pSigma0.Update();

        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Fitter
        DLM_Fitter1* fitter = new DLM_Fitter1(1);
        fitter->SetSystem(0, *dataHist, 1, CkDec_pSigma0,
                          FemtoRegion_pSigma[FitReg_pSigma][0],
                          FemtoRegion_pSigma[FitReg_pSigma][1], 1000, 1000);
        fitter->SetSeparateBL(0, false);  //Simultaneous BL
        if (useBaselineSlope) {
          fitter->FixParameter("pSigma0", DLM_Fitter1::p_a,
                               p_a_prefit[BaselineSlope]);
          fitter->FixParameter("pSigma0", DLM_Fitter1::p_b,
                               p_b_prefit[BaselineSlope]);
          std::cout << "Use slope for baseline\n";
        } else {
          fitter->FixParameter("pSigma0", DLM_Fitter1::p_a, p_a_prefit[0]);
          fitter->FixParameter("pSigma0", DLM_Fitter1::p_b, p_b_prefit[0]);
        }
        fitter->AddSameSource("pSigma0SideBand", "pSigma0", 1);

        //Fit BL & Normalization
        fitter->FixParameter("pSigma0", DLM_Fitter1::p_c, 0);
        fitter->FixParameter("pSigma0", DLM_Fitter1::p_Cl, -1.);
        fitter->FixParameter("pSigma0", DLM_Fitter1::p_sor0,
                             sourceSize[SourceSizeIter]);

        fitter->GoBabyGo();

        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Get the parameters from the fit
        const double bl_a = fitter->GetParameter("pSigma0", DLM_Fitter1::p_a);
        const double bl_a_err = fitter->GetParError("pSigma0",
                                                    DLM_Fitter1::p_a);
        const double bl_b = fitter->GetParameter("pSigma0", DLM_Fitter1::p_b);
        const double bl_b_err = fitter->GetParError("pSigma0",
                                                    DLM_Fitter1::p_b);
        const double Cl = fitter->GetParameter("pSigma0", DLM_Fitter1::p_c);
        const double chi2 = fitter->GetChi2Ndf();
        const double pval = fitter->GetPval();

        TGraph FitResult_pSigma0;
        FitResult_pSigma0.SetName(TString::Format("pSigma0Graph"));
        fitter->GetFitGraph(0, FitResult_pSigma0);

        double Chi2_pSigma0 = 0;
        double EffNumBins_pSigma0 = 0;
        int maxkStarBin = dataHist->FindBin(250);
        for (unsigned uBin = 1; uBin <= maxkStarBin; uBin++) {

          double mom = dataHist->GetBinCenter(uBin);
          //double dataX;
          double dataY;
          double dataErr;
          double theoryX;
          double theoryY;

          if (mom > FemtoRegion_pSigma[FitReg_pSigma][1]) {
            continue;
          }

          FitResult_pSigma0.GetPoint(uBin - 1, theoryX, theoryY);
          if (mom != theoryX) {
            std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX
                      << std::endl;
          }
          dataY = dataHist->GetBinContent(uBin);
          dataErr = dataHist->GetBinError(uBin);
          Chi2_pSigma0 += (dataY - theoryY) * (dataY - theoryY)
              / (dataErr * dataErr);
          ++EffNumBins_pSigma0;
        }
        double pvalpSigma0 = TMath::Prob(Chi2_pSigma0,
                                         round(EffNumBins_pSigma0));
        double nSigmapSigma0 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalpSigma0);

        std::cout << "=============\n";
        std::cout << "Fitter output\n";
        std::cout << "BL a  " << bl_a << " " << bl_a_err << "\n";
        std::cout << "BL b  " << bl_b << " " << bl_b_err << "\n";
        std::cout << "Cl    " << Cl << "\n";
        std::cout << "Chi2\n";
        std::cout << " glob " << chi2 << "\n";
        std::cout << " loc  " << Chi2_pSigma0 / round(EffNumBins_pSigma0)
                  << "\n";
        std::cout << "p-val\n";
        std::cout << " glob " << pval << "\n";
        std::cout << " loc  " << pvalpSigma0 << "\n";
        std::cout << "=============\n";

        /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        /// Write out all the stuff

        TGraph grCFSigmaRaw;
        grCFSigmaRaw.SetName("Sigma0Raw");
        TGraph grCFSigmaMain;
        grCFSigmaMain.SetName("Sigma0Main");
        TGraph grCFSigmaFeed;
        grCFSigmaFeed.SetName("Sigma0Feed");
        TGraph grCFSigmaSideband;
        grCFSigmaSideband.SetName("Sigma0Sideband");
        for (unsigned int i = 0; i < Ck_pSigma0->GetNbins(); ++i) {
          grCFSigmaRaw.SetPoint(
              i, Ck_pSigma0->GetBinCenter(i),
              CkDec_pSigma0.EvalCk(Ck_pSigma0->GetBinCenter(i)));
          grCFSigmaMain.SetPoint(
              i,
              Ck_pSigma0->GetBinCenter(i),
              ((CkDec_pSigma0.EvalMain(Ck_pSigma0->GetBinCenter(i)) - 1.)
                  * lambdaParams.GetLambdaParam(CATSLambdaParam::Primary)) + 1);
          grCFSigmaFeed.SetPoint(
              i, Ck_pSigma0->GetBinCenter(i),
              CkDec_pSigma0.EvalMainFeed(Ck_pSigma0->GetBinCenter(i)));
          grCFSigmaSideband.SetPoint(
              i,
              Ck_pSigma0->GetBinCenter(i),
              ((Ck_SideBand->Eval(Ck_pSigma0->GetBinCenter(i)) - 1.)
                  * lambdaParams.GetLambdaParam(CATSLambdaParam::Fake)) + 1);
        }

        TString outfilename = TString::Format("%s/Graph_pSigma0_%i_%i.root",
                                              OutputDir.Data(), NumIter,
                                              iterID);
        auto out = new TFile(outfilename, "RECREATE");

        /// beautification
        DreamPlot::SetStyleHisto(dataHist, 24, kBlack);
        dataHist->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
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

        grPrefitContour->Write("fitContour");
        grCFSigmaRaw.Write();
        grCFSigmaMain.Write();
        grCFSigmaFeed.Write();
        grCFSigmaSideband.Write();
        sidebandHistLow->Write("SidebandLow");
        sidebandHistUp->Write("SidebandUp");
        sideband->Write("SidebandFitNotScaled");
        SBmerge->Write("SidebandMerged");
        dataHist->Write("CF");
        FitResult_pSigma0.Write("Fit");

        if (fastPlot || iterID == 0) {
          auto c = new TCanvas();
          dataHist->GetXaxis()->SetRangeUser(0., 600);
          dataHist->GetYaxis()->SetRangeUser(0.8, 1.6);
          dataHist->Draw();
          FitResult_pSigma0.Draw("l3same");
          grCFSigmaSideband.Draw("l3same");

          auto info = new TPaveText(0.5, 0.5, 0.88, 0.85, "blNDC");
          info->SetBorderSize(0);
          info->SetTextSize(0.04);
          info->SetFillColor(kWhite);
          info->SetTextFont(42);
          TString SOURCE_NAME = "Gauss";
          double Yoffset = 1.2;
          info->AddText(
              TString::Format(
                  "r_{%s} = %.3f #pm %.3f fm", SOURCE_NAME.Data(),
                  fitter->GetParameter("pSigma0", DLM_Fitter1::p_sor0),
                  fitter->GetParError("pSigma0", DLM_Fitter1::p_sor0)));
          info->AddText(TString::Format("a = %.3f #pm %.3f", bl_a, bl_a_err));

          if (useBaselineSlope) {
            info->AddText(
                TString::Format("b = (%.3f #pm %.3f ) #times 10^{-4}",
                                bl_b * 1e4, bl_b_err * 1e4));
          }
          info->AddText(
              TString::Format("Local #chi^{2}/ndf=%.1f/%.0f = %.2f", Chi2_pSigma0,
                              EffNumBins_pSigma0,
                              Chi2_pSigma0 / double(EffNumBins_pSigma0)));
          info->AddText(
              TString::Format("p_{val}=%.3f, n_{#sigma}=%.3f", pvalpSigma0,
                              nSigmapSigma0));

          info->Draw("same");
          c->Write("CFplot");
          c->Print(Form("%s/CF_pSigma0.pdf", OutputDir.Data()));

          auto d = new TCanvas();
          SBmerge->GetXaxis()->SetRangeUser(0., 600);
          SBmerge->GetYaxis()->SetRangeUser(0.8, 1.6);
          SBmerge->Draw();
          sideband->Draw("l3same");
          d->Write("CFsideband");
          d->Print(Form("%s/CF_pSideband.pdf", OutputDir.Data()));

          auto e = new TCanvas();
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
          e->Print(Form("%s/CF_pSideband_all.pdf", OutputDir.Data()));
        }
        out->Close();

        param->cd();
        ntBuffer[0] = iterID;
        ntBuffer[1] = FemtoRegion_pSigma[FitReg_pSigma][1];
        ntBuffer[2] = (float) BaselineSlope;
        ntBuffer[3] = sourceSize[SourceSizeIter];
        ntBuffer[4] = bl_a;
        ntBuffer[5] = bl_a_err;
        ntBuffer[6] = bl_b;
        ntBuffer[7] = bl_b_err;
        ntBuffer[8] = chi2;
        ntBuffer[9] = pval;
        ntBuffer[10] = Chi2_pSigma0;
        ntBuffer[11] = (float) EffNumBins_pSigma0;
        ntBuffer[12] = Chi2_pSigma0 / double(EffNumBins_pSigma0);
        ntBuffer[13] = pvalpSigma0;
        ntBuffer[14] = nSigmapSigma0;

//  ntBuffer[2] = vFrac_pXim_pXi1530;
//  ntBuffer[3] = tOut;
//  ntBuffer[5] = varSideNorm;

        ntResult->Fill(ntBuffer);
        ++iterID;

        delete Ck_SideBand;
        delete Ck_pSigma0;
        delete fitter;
        delete out;

        if (NumIter == 0) {
          goto exitThroughTheGiftShop;
        }
      }
    }
  }
  exitThroughTheGiftShop: ntResult->Write();
  param->Close();
  return;
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  FitSigma0(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5]);
  return 0;
}
