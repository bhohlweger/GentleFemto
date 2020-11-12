#include "TApplication.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

#include "TidyCats.h"
#include "DreamPlot.h"
#include "CATSLambdaParam.h"
#include "CATSInput.h"
#include "TStyle.h"

#include "DLM_CkDecomposition.h"

/// =====================================================================================
void FillSourceGraph(CATS &kitty, TGraph *gr) {
  for (double i = 0; i < 100; ++i) {
    gr->SetPoint(i, i * 0.2, kitty.EvaluateTheSource(0, i * 0.2, 0));
  }
}

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}
/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, DLM_CkDecomposition &ckGraph, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ckGraph.EvalCk(mom));
  }
}

/// =====================================================================================
void FillCkGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  TApplication *app = new TApplication("app", 0, 0);
  DreamPlot::SetStyle();

  /// -----------------------------------------------------------------------------------
  /// Parameters

  const double kmin = 0;
  const double kmax = 500;
  const int nBins = 50;

  /// Lambda parameters
  const double protonPurity = 0.994;
  const double protonPrimary = 0.823;
  const double protonLambda = 0.125;

  const double DmesonPurity = 0.35;   // rough estimate, to be improved
  const double Bfeeddown = 0.9;  // rough estimate, to be improved TODO - charged fractions?
  const double DstarFeeding = 1. / 3.;  // rough estimate, to be improved
  const double DmesonPrimary = Bfeeddown * 1.f / (1.f + DstarFeeding);
  const double DstarContribution = Bfeeddown * DstarFeeding
      / (1.f + DstarFeeding);

  const Particle dmeson(DmesonPurity, DmesonPrimary,
                        { { DstarContribution, 1.f - Bfeeddown } });
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonLambda, (1. - protonPrimary)
          * (1 - protonLambda) } });

  const CATSLambdaParam lambdaParam(proton, dmeson);
  lambdaParam.PrintLambdaParams();

  const double primary = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary);
  const double pDstar = lambdaParam.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   0);

  /// Femto radius
  const double rCoreLow = 0.86;
  const double rCoreDefault = 0.91;
  const double rCoreUp = 0.96;

  /// -----------------------------------------------------------------------------------
  /// Calibration
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  auto calibFile = TFile::Open(
      TString::Format("%s/dstar.root", CalibBaseDir.Data()));
  auto histDecayKindematicsDstar = (TH2F*) calibFile->Get("histSmearDmeson");
  DreamPlot::SetStyleHisto(histDecayKindematicsDstar);
  auto asdasd = new TCanvas("Decay momentum D*+", "Decay momentum D*+");
  asdasd->SetRightMargin(0.16);
  histDecayKindematicsDstar->Draw("colz");

  // TODO Momentum resolution
  TH2F* histMomentumResolution = nullptr;
//  DreamPlot::SetStyleHisto(histMomentumResolution);

  /// -----------------------------------------------------------------------------------

  auto gaussFit =
      new TF1(
          "gaus",
          [&](double *x, double *p) {return
            4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
            std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
          },
          0, 10, 1);
  gaussFit->SetParameter(0, 1.15);
  gaussFit->SetNpx(1000);
  gaussFit->SetLineColorAlpha(kBlue + 3, 0.8);
  gaussFit->SetLineWidth(2);
  gaussFit->SetLineStyle(2);

  /// -----------------------------------------------------------------------------------
  // For the Dplus

  auto grSourceDplus = new TGraph();
  DreamPlot::SetStyleGraph(grSourceDplus, 20, kBlue + 3, 0.8);
  grSourceDplus->SetLineWidth(2);
  grSourceDplus->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");

  auto grDplusCoulomb = new TGraph();
  DreamPlot::SetStyleGraph(grDplusCoulomb, 20, kBlue + 3, 0.8);
  grDplusCoulomb->SetLineWidth(2);
  grDplusCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDplusCoulombLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDplusCoulombLambda, 20, kBlue + 3, 0.8);
  grDplusCoulombLambda->SetLineWidth(2);
  grDplusCoulombLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  auto grDstarplusGenuine = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusGenuine, 20, kBlue + 3, 0.8);
  grDstarplusGenuine->SetLineWidth(2);
  grDstarplusGenuine->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarplusSmeared = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusSmeared, 20, kGreen + 3, 0.8);
  grDstarplusSmeared->SetLineWidth(2);
  grDstarplusSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarplusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusLambda, 20, kRed + 2, 0.8);
  grDstarplusLambda->SetLineWidth(2);
  grDstarplusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  CATS catsDplusCoulombOnly, catsDstarplus;
  TidyCats *tidyCats = new TidyCats();
  tidyCats->GetCatsProtonDplus(&catsDplusCoulombOnly, nBins, kmin, kmax,
                               TidyCats::pCoulombOnly, TidyCats::sResonance);
  tidyCats->GetCatsProtonDstarplus(&catsDstarplus, nBins, kmin, kmax,
                                   TidyCats::pCoulombOnly,
                                   TidyCats::sResonance);

  catsDplusCoulombOnly.SetAnaSource(0, rCoreDefault);
  catsDplusCoulombOnly.KillTheCat();
  catsDstarplus.SetAnaSource(0, rCoreDefault);
  catsDstarplus.KillTheCat();

  FillSourceGraph(catsDplusCoulombOnly, grSourceDplus);
  FillCkGraph(catsDplusCoulombOnly, grDplusCoulomb);

  auto cPlus = new TCanvas("p-D+ source", "p-D+ source");
  grSourceDplus->Draw("AL3");
  grSourceDplus->Fit(gaussFit, "Q", "R", 0, 12);
  std::cout << "p-D+ radius " << gaussFit->GetParameter(0) << " fm\n";

  auto dPlus = new TCanvas("p-D+ Coulomb", "p-D+ Coulomb");
  grDplusCoulomb->Draw("AL3");

  auto DLM_pDplusCoulomb = new DLM_Ck(1, 0, catsDplusCoulombOnly);
  auto DLM_pDstarplus = new DLM_Ck(1, 0, catsDstarplus);

  DLM_CkDecomposition CkDec_pDplusCoulomb("pDplusCoulomb", 2,
                                          *DLM_pDplusCoulomb,
                                          histMomentumResolution);
  DLM_CkDecomposition CkDec_pDstarplus("pDstarplus", 0, *DLM_pDstarplus,
                                       nullptr);
  CkDec_pDplusCoulomb.AddContribution(0, pDstar, DLM_CkDecomposition::cFeedDown,
                                      &CkDec_pDstarplus,
                                      histDecayKindematicsDstar);
  CkDec_pDplusCoulomb.AddContribution(1, 1.f - primary - pDstar,
                                      DLM_CkDecomposition::cFake);
  CkDec_pDplusCoulomb.Update();
  FillCkGraph(DLM_pDplusCoulomb, CkDec_pDplusCoulomb, grDplusCoulombLambda);

  DLM_CkDecomposition CkDec_pDstarplus_Smeared("pDstarplus_smear", 0,
                                               *DLM_pDstarplus,
                                               histDecayKindematicsDstar);
  DLM_CkDecomposition CkDec_pDstarplus_Lambda("pDstarplus_lambda", 1,
                                              *DLM_pDstarplus,
                                              histDecayKindematicsDstar);
  CkDec_pDstarplus_Lambda.AddContribution(0, 1.f - pDstar,
                                          DLM_CkDecomposition::cFake);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus_Smeared, grDstarplusSmeared);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus, grDstarplusGenuine);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus_Lambda, grDstarplusLambda);

  auto dPlusLambda = new TCanvas("p-D+ Experimental", "p-D+ Experimental");
  grDplusCoulombLambda->Draw("AL3");

  auto dDstarPlus = new TCanvas("p-D*+ Transformation", "p-D*+ Transformation");
  grDstarplusSmeared->Draw("AL3");
  grDstarplusGenuine->Draw("L3");
  grDstarplusLambda->Draw("L3");
  auto legDstarplus = new TLegend(0.45, 0.25, 0.45 + 0.3, 0.55);
  legDstarplus->SetBorderSize(0);
  legDstarplus->SetTextFont(42);
  legDstarplus->SetHeader("p#minus#kern[-0.95]{ }D*^{+} Coulomb");
  legDstarplus->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDstarplus->AddEntry(grDstarplusGenuine, "Genuine", "l");
  legDstarplus->AddEntry(
      grDstarplusSmeared,
      "p#minus#kern[-0.95]{ }D*^{+} #rightarrow p#minus#kern[-0.95]{ }D^{+}",
      "l");
  legDstarplus->AddEntry(
      grDstarplusLambda,
      TString::Format(
          "p#minus#kern[-0.95]{ }D*^{+} #rightarrow p#minus#kern[-0.95]{ }D^{+} (#lambda = %.2f %%)",
          pDstar * 100.),
      "l");
  legDstarplus->Draw("same");

  /// -----------------------------------------------------------------------------------
  /// For the Dminus

  auto grSourceDminus = new TGraph();
  DreamPlot::SetStyleGraph(grSourceDminus, 20, kBlue + 3, 0.8);
  grSourceDminus->SetLineWidth(2);
  grSourceDminus->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");

  auto grDminusCoulomb = new TGraph();
  DreamPlot::SetStyleGraph(grDminusCoulomb, 20, kBlue + 3, 0.8);
  grDminusCoulomb->SetLineWidth(2);
  grDminusCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDminusCoulombLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDminusCoulombLambda, 20, kBlue + 3, 0.8);
  grDminusCoulombLambda->SetLineWidth(2);
  grDminusCoulombLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  auto grDstarminusGenuine = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusGenuine, 20, kBlue + 3, 0.8);
  grDstarminusGenuine->SetLineWidth(2);
  grDstarminusGenuine->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarminusSmeared = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusSmeared, 20, kGreen + 3, 0.8);
  grDstarminusSmeared->SetLineWidth(2);
  grDstarminusSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarminusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusLambda, 20, kRed + 2, 0.8);
  grDstarminusLambda->SetLineWidth(2);
  grDstarminusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  CATS catsDminusCoulombOnly, catsDstarminus;
  tidyCats->GetCatsProtonDminus(&catsDminusCoulombOnly, nBins, kmin, kmax,
                               TidyCats::pCoulombOnly, TidyCats::sResonance);
  tidyCats->GetCatsProtonDstarminus(&catsDstarminus, nBins, kmin, kmax,
                                   TidyCats::pCoulombOnly,
                                   TidyCats::sResonance);

  catsDminusCoulombOnly.SetAnaSource(0, rCoreDefault);
  catsDminusCoulombOnly.KillTheCat();
  catsDstarminus.SetAnaSource(0, rCoreDefault);
  catsDstarminus.KillTheCat();

  FillSourceGraph(catsDminusCoulombOnly, grSourceDminus);
  FillCkGraph(catsDminusCoulombOnly, grDminusCoulomb);

  auto cminus = new TCanvas("p-D- source", "p-D- source");
  grSourceDminus->Draw("AL3");
  grSourceDminus->Fit(gaussFit, "Q", "R", 0, 12);
  std::cout << "p-D- radius " << gaussFit->GetParameter(0) << " fm\n";

  auto dminus = new TCanvas("p-D- Coulomb", "p-D- Coulomb");
  grDminusCoulomb->Draw("AL3");

  auto DLM_pDminusCoulomb = new DLM_Ck(1, 0, catsDminusCoulombOnly);
  auto DLM_pDstarminus = new DLM_Ck(1, 0, catsDstarminus);

  DLM_CkDecomposition CkDec_pDminusCoulomb("pDminusCoulomb", 2,
                                          *DLM_pDminusCoulomb,
                                          histMomentumResolution);
  DLM_CkDecomposition CkDec_pDstarminus("pDstarminus", 0, *DLM_pDstarminus,
                                       nullptr);
  CkDec_pDminusCoulomb.AddContribution(0, pDstar, DLM_CkDecomposition::cFeedDown,
                                      &CkDec_pDstarminus,
                                      histDecayKindematicsDstar);
  CkDec_pDminusCoulomb.AddContribution(1, 1.f - primary - pDstar,
                                      DLM_CkDecomposition::cFake);
  CkDec_pDminusCoulomb.Update();
  FillCkGraph(DLM_pDminusCoulomb, CkDec_pDminusCoulomb, grDminusCoulombLambda);

  DLM_CkDecomposition CkDec_pDstarminus_Smeared("pDstarminus_smear", 0,
                                               *DLM_pDstarminus,
                                               histDecayKindematicsDstar);
  DLM_CkDecomposition CkDec_pDstarminus_Lambda("pDstarminus_lambda", 1,
                                              *DLM_pDstarminus,
                                              histDecayKindematicsDstar);
  CkDec_pDstarminus_Lambda.AddContribution(0, 1.f - pDstar,
                                          DLM_CkDecomposition::cFake);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus_Smeared, grDstarminusSmeared);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus, grDstarminusGenuine);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus_Lambda, grDstarminusLambda);

  auto dminusLambda = new TCanvas("p-D- Experimental", "p-D- Experimental");
  grDminusCoulombLambda->Draw("AL3");

  auto dDstar = new TCanvas("p-D*- Transformation", "p-D*- Transformation");
  grDstarminusSmeared->Draw("AL3");
  grDstarminusGenuine->Draw("L3");
  grDstarminusLambda->Draw("L3");
  auto legDstarminus = new TLegend(0.45, 0.25, 0.45 + 0.3, 0.55);
  legDstarminus->SetBorderSize(0);
  legDstarminus->SetTextFont(42);
  legDstarminus->SetHeader("p#minus#kern[-0.95]{ }D*^{#minus} Coulomb");
  legDstarminus->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDstarminus->AddEntry(grDstarminusGenuine, "Genuine", "l");
  legDstarminus->AddEntry(
      grDstarminusSmeared,
      "p#minus#kern[-0.95]{ }D*^{#minus} #rightarrow p#minus#kern[-0.95]{ }D^{#minus}",
      "l");
  legDstarminus->AddEntry(
      grDstarminusLambda,
      TString::Format(
          "p#minus#kern[-0.95]{ }D*^{#minus} #rightarrow p#minus#kern[-0.95]{ }D^{#minus} (#lambda = %.2f %%)",
          pDstar * 100.),
      "l");
  legDstarminus->Draw("same");
  app->Run();
}
