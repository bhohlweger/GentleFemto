#include "TApplication.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"

#include "TidyCats.h"
#include "DreamPlot.h"
#include "CATSLambdaParam.h"
#include "CATSInput.h"
#include "TStyle.h"
#include "TLatex.h"

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
TH2F* TransformToMeV(const TH2F* input) {
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
int main(int argc, char *argv[]) {
  TApplication *app = new TApplication("app", 0, 0);
  DreamPlot::SetStyle();

  /// -----------------------------------------------------------------------------------
  /// Parameters

  const double kmin = 0;
  const double kmax = 300;
  const int nBins = 30;

  /// Lambda parameters
  const double protonPurityDplus = 0.980983;
  const double protonPrimaryDplus = 0.863268;
  const double protonLambdaDplus = 0.0966698;

  const double protonPurityDminus = 0.983773;
  const double protonPrimaryDminus = 0.862303;
  const double protonLambdaDminus = 0.0973517;

  const double DmesonPurity = 0.61669272;
  const double Bfeeddown = 0.076733507;
  const double DstarFeeding = 0.276;
  const double DmesonPrimary = 1.f - Bfeeddown - DstarFeeding;
  const Particle dmeson(DmesonPurity, DmesonPrimary,
                        { { DstarFeeding,  Bfeeddown} });

  const Particle protonDplus(
      protonPurityDplus,
      protonPrimaryDplus,
      { { (1. - protonPrimaryDplus) * protonLambdaDplus, (1. - protonPrimaryDplus)
          * (1 - protonLambdaDplus) } });
  const Particle protonDminus(
      protonPurityDminus,
      protonPrimaryDminus,
      { { (1. - protonPrimaryDminus) * protonLambdaDminus, (1. - protonPrimaryDminus)
          * (1 - protonLambdaDminus) } });

  const CATSLambdaParam lambdaParampDplus(protonDplus, dmeson);
  lambdaParampDplus.PrintLambdaParams();

  const double primaryDplus = lambdaParampDplus.GetLambdaParam(CATSLambdaParam::Primary);
  const double pDstarDplus = lambdaParampDplus.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   0);
  const double flatDplus = 1.f - primaryDplus - pDstarDplus;
  
  std::cout << "Genuine p-D+ " << primaryDplus << "\n";
  std::cout << "p-D* -> p-D+ " << pDstarDplus << "\n";
  std::cout << "Flat in p-D+ " << flatDplus << "\n";

  const CATSLambdaParam lambdaParampDminus(protonDminus, dmeson);
  lambdaParampDminus.PrintLambdaParams();

  const double primaryDminus = lambdaParampDminus.GetLambdaParam(CATSLambdaParam::Primary);
  const double pDstarDminus = lambdaParampDminus.GetLambdaParam(CATSLambdaParam::Primary,
                                                   CATSLambdaParam::FeedDown, 0,
                                                   0);
  const double flatDminus = 1.f - primaryDminus - pDstarDminus;

  std::cout << "Genuine p-D- " << primaryDminus << "\n";
  std::cout << "p-D* -> p-D- " << pDstarDminus << "\n";
  std::cout << "Flat in p-D- " << flatDminus << "\n";
  
  /// Femto radius
  const double rCorepDplusLow = 0.72;
  const double rCorepDplusDefault = 0.8;
  const double rCorepDplusUp = 0.88;

  const double rCorepDminusLow = 0.72;
  const double rCorepDminusDefault = 0.80;
  const double rCorepDminusUp = 0.88;

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
  asdasd->Print("DecayTransformationDstar.pdf");
  
  auto calibFileBeauty = TFile::Open(
      TString::Format("%s/beauty.root", CalibBaseDir.Data()));
  auto histDecayKindematicsBeauty = (TH2F*) calibFileBeauty->Get("histSmearBmeson");
  auto beautyCan = new TCanvas();
  beautyCan->SetRightMargin(0.16);
  histDecayKindematicsBeauty->Draw("colz");
  beautyCan->Print("DecayTransformationB.pdf");
  
  // Momentum resolution
  auto momResFile = TFile::Open(TString::Format("%s/momRes_Dmesons.root", CalibBaseDir.Data()).Data());
  auto histMomentumResolutionpDplus = TransformToMeV((TH2F*)momResFile->Get("pDplus"));
  auto histMomentumResolutionpDminus = TransformToMeV((TH2F*)momResFile->Get("pDminus"));
  DreamPlot::SetStyleHisto(histMomentumResolutionpDplus);
  DreamPlot::SetStyleHisto(histMomentumResolutionpDminus);

  // since the resolution is the same for p-D+ and p-D- we add them (see control histogram!)
  auto histMomentumResolution = histMomentumResolutionpDplus;
  histMomentumResolution->Add(histMomentumResolutionpDminus);

  auto cMomRes1 = new TCanvas("Momentum resolution", "Momentum resolution");
  histMomentumResolution->Draw("colz");

  auto cMomRes2 = new TCanvas("Mom. res. - ratio p-Dplus p-Dminus", "Mom. res. - ratio p-Dplus p-Dminus");
  auto ratio = (TH2F*)histMomentumResolutionpDminus->Clone("ratp");
  ratio->Divide(histMomentumResolutionpDplus);
  ratio->Draw("colz");

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
  auto grDplusTotal = new TGraph();
  DreamPlot::SetStyleGraph(grDplusTotal, 20, kBlue + 3, 0.8);
  grDplusTotal->SetLineWidth(2);
  grDplusTotal->SetLineStyle(2);
  grDplusTotal->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDplusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDplusLambda, 20, kBlue + 3, 0.8);
  grDplusLambda->SetLineWidth(2);
  grDplusLambda->SetLineStyle(2);
  grDplusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  
  auto grDstarplusGenuine = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusGenuine, 20, kGreen + 3, 0.8);
  grDstarplusGenuine->SetLineWidth(2);
  grDstarplusGenuine->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarplusSmeared = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusSmeared, 20, kGreen + 3, 0.8);
  grDstarplusSmeared->SetLineWidth(2);
  grDstarplusSmeared->SetLineStyle(3);
  grDstarplusSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarplusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDstarplusLambda, 20, kGreen + 3, 0.8);
  grDstarplusLambda->SetLineWidth(2);
  grDstarplusLambda->SetLineStyle(2);
  grDstarplusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  CATS catsDplusCoulombOnly, catsDstarplus;
  TidyCats *tidyCats = new TidyCats();
  tidyCats->GetCatsProtonDplus(&catsDplusCoulombOnly, nBins, kmin, kmax,
                               TidyCats::pDCoulombOnly, TidyCats::sResonance);
  tidyCats->GetCatsProtonDstarplus(&catsDstarplus, nBins, kmin, kmax,
                                   TidyCats::pDCoulombOnly,
                                   TidyCats::sResonance);

  catsDplusCoulombOnly.SetAnaSource(0, rCorepDplusDefault);
  catsDplusCoulombOnly.KillTheCat();
  catsDstarplus.SetAnaSource(0, rCorepDplusDefault);
  catsDstarplus.KillTheCat();

  FillSourceGraph(catsDplusCoulombOnly, grSourceDplus);
  FillCkGraph(catsDplusCoulombOnly, grDplusCoulomb);

  auto cPlus = new TCanvas("p-D+ source", "p-D+ source");
  grSourceDplus->Draw("AL3");
  grSourceDplus->Fit(gaussFit, "Q", "R", 0, 12);
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(23);
  text.SetTextFont(43);
  const double dplusRad = gaussFit->GetParameter(0);
  text.DrawLatex(0.4, 0.8, Form("#it{r}_{eff} = %.2f fm", dplusRad));
  cPlus->Print("SourceDplus.pdf");
  std::cout << "p-D+ radius " << dplusRad << " fm\n";

  auto DLM_pDplusCoulomb = new DLM_Ck(1, 0, catsDplusCoulombOnly);
  auto DLM_pDstarplus = new DLM_Ck(1, 0, catsDstarplus);

  DLM_CkDecomposition CkDec_pDplusCoulomb("pDplusCoulomb", 2,
                                          *DLM_pDplusCoulomb,
                                          histMomentumResolution);
  DLM_CkDecomposition CkDec_pDstarplus("pDstarplus", 0, *DLM_pDstarplus,
                                       nullptr);
  CkDec_pDplusCoulomb.AddContribution(0, pDstarDplus, DLM_CkDecomposition::cFeedDown,
                                      &CkDec_pDstarplus,
                                      histDecayKindematicsDstar);
  CkDec_pDplusCoulomb.AddContribution(1, flatDplus,
                                      DLM_CkDecomposition::cFake);
  CkDec_pDplusCoulomb.Update();
  FillCkGraph(DLM_pDplusCoulomb, CkDec_pDplusCoulomb, grDplusTotal);

  DLM_CkDecomposition CkDec_pDplusLambda("pDplusCoulomb", 1,
                                         *DLM_pDplusCoulomb,
                                         histMomentumResolution);
  CkDec_pDplusLambda.AddContribution(0, primaryDplus,
                                     DLM_CkDecomposition::cFake);
  CkDec_pDplusLambda.Update();
  FillCkGraph(DLM_pDplusCoulomb, CkDec_pDplusLambda, grDplusLambda);

  DLM_CkDecomposition CkDec_pDstarplus_Smeared("pDstarplus_smear", 0,
                                               *DLM_pDstarplus,
                                               histDecayKindematicsDstar);
  DLM_CkDecomposition CkDec_pDstarplus_Lambda("pDstarplus_lambda", 1,
                                              *DLM_pDstarplus,
                                              histDecayKindematicsDstar);
  CkDec_pDstarplus_Lambda.AddContribution(0, 1.f - pDstarDplus,
                                          DLM_CkDecomposition::cFake);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus_Smeared, grDstarplusSmeared);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus, grDstarplusGenuine);
  FillCkGraph(DLM_pDstarplus, CkDec_pDstarplus_Lambda, grDstarplusLambda);


  auto dPlusLambda = new TCanvas("p-D+ Experimental", "p-D+ Experimental");
  grDplusTotal->Draw("AL3");
  auto legDplus2 = new TLegend(0.45, 0.25, 0.45 + 0.3, 0.25+0.125);
  legDplus2->SetBorderSize(0);
  legDplus2->SetTextFont(42);
  legDplus2->SetHeader(
      TString::Format("p#minus#kern[-0.95]{ }D^{+}, #it{r}_{eff} = %.2f fm",
                      dplusRad));
  legDplus2->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDplus2->AddEntry(grDplusTotal,
                      "p#minus#kern[-0.95]{ }D^{+} (Coulomb only)", "l");
  legDplus2->Draw("same");
  dPlusLambda->Print("dplusExp.pdf");
  
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
          pDstarDplus * 100.),
      "l");
  legDstarplus->Draw("same");
  dDstarPlus->Print("Transformation_Dstar_Dplus.pdf");

  auto cDplusTotal = new TCanvas("p-D+ total", "p-D+ total");
  grDplusTotal->Draw("AL3");
  grDplusTotal->GetXaxis()->SetNdivisions(504);
  grDplusTotal->GetYaxis()->SetNdivisions(504);
  grDplusTotal->SetLineStyle(0);
  grDplusTotal->GetXaxis()->SetRangeUser(0, 300);
  grDplusTotal->GetYaxis()->SetRangeUser(0.7, 1.1);
  grDplusLambda->Draw("L3");
  grDstarplusLambda->Draw("L3");
  auto legDplus = new TLegend(0.45, 0.25, 0.45 + 0.3, 0.55);
  legDplus->SetBorderSize(0);
  legDplus->SetTextFont(42);
  legDplus->SetHeader(TString::Format("p#minus#kern[-0.95]{ }D^{+}, #it{r}_{eff} = %.2f fm", dplusRad));
  legDplus->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDplus->AddEntry(grDplusTotal, "Total", "l");
  legDplus->AddEntry(
      grDplusLambda,
      TString::Format("Genuine p#minus#kern[-0.95]{ }D^{+} (#lambda = %.2f %%)", primaryDplus * 100.f),
      "l");
  legDplus->AddEntry(
      grDstarplusLambda,
      TString::Format(
          "p#minus#kern[-0.95]{ }D*^{+} #rightarrow p#minus#kern[-0.95]{ }D^{+} (#lambda = %.2f %%)",
          pDstarDplus * 100.),
      "l");
  legDplus->Draw("same");
  cDplusTotal->Print("p-Dplus_total.pdf");

  auto dPlus = new TCanvas("p-D+ Coulomb", "p-D+ Coulomb");
  grDplusCoulomb->Draw("AL3");
  grDplusLambda->Draw("L3");
  dPlus->Print("dplusCoul.pdf");
  
  /// -----------------------------------------------------------------------------------
  /// For the Dminus

  auto grSourceDminus = new TGraph();
  DreamPlot::SetStyleGraph(grSourceDminus, 20, kBlue + 3, 0.8);
  grSourceDminus->SetLineWidth(2);
  grSourceDminus->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");

  auto grSourceDminusUp = new TGraph();
  auto grSourceDminusLow = new TGraph();
  
  auto grDminusCoulomb = new TGraph();
  DreamPlot::SetStyleGraph(grDminusCoulomb, 20, kBlue + 3, 0.8);
  grDminusCoulomb->SetLineWidth(2);
  grDminusCoulomb->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDminusHaidenbauer = new TGraph();
  DreamPlot::SetStyleGraph(grDminusHaidenbauer, 20, kRed + 3, 0.8);
  grDminusHaidenbauer->SetLineWidth(2);
  grDminusHaidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDminusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDminusLambda, 20, kBlue + 3, 0.8);
  grDminusLambda->SetLineWidth(2);
  grDminusLambda->SetLineStyle(2);
  grDminusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDminusTotal = new TGraph();
  DreamPlot::SetStyleGraph(grDminusTotal, 20, kBlue + 3, 0.8);
  grDminusTotal->SetLineWidth(2);
  grDminusTotal->SetLineStyle(2);
  grDminusTotal->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  auto grDstarminusGenuine = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusGenuine, 20, kGreen + 3, 0.8);
  grDstarminusGenuine->SetLineWidth(2);
  grDstarminusGenuine->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarminusSmeared = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusSmeared, 20, kGreen + 3, 0.8);
  grDstarminusSmeared->SetLineWidth(2);
  grDstarminusSmeared->SetLineStyle(3);
  grDstarminusSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grDstarminusLambda = new TGraph();
  DreamPlot::SetStyleGraph(grDstarminusLambda, 20, kGreen + 3, 0.8);
  grDstarminusLambda->SetLineWidth(2);
  grDstarminusLambda->SetLineStyle(2);
  grDstarminusLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");

  CATS catsDminusCoulombOnly, catsDstarminus, catsDminusHaidenbauer;
  tidyCats->GetCatsProtonDminus(&catsDminusCoulombOnly, nBins, kmin, kmax,
                               TidyCats::pDCoulombOnly, TidyCats::sResonance);
  tidyCats->GetCatsProtonDstarminus(&catsDstarminus, nBins, kmin, kmax,
                                   TidyCats::pDCoulombOnly,
                                   TidyCats::sResonance);
  tidyCats->GetCatsProtonDminus(&catsDminusHaidenbauer, nBins, kmin, kmax,
                               TidyCats::pDminusHaidenbauer, TidyCats::sResonance);

  catsDminusCoulombOnly.SetAnaSource(0, rCorepDminusLow);
  catsDminusCoulombOnly.KillTheCat();
  FillSourceGraph(catsDminusCoulombOnly, grSourceDminusLow);
  grSourceDminusLow->Fit(gaussFit, "Q", "R", 0, 12);
  const double dminusRadLow = gaussFit->GetParameter(0);
  catsDminusCoulombOnly.SetAnaSource(0, rCorepDminusUp);
  catsDminusCoulombOnly.KillTheCat();
  FillSourceGraph(catsDminusCoulombOnly, grSourceDminusUp);
  grSourceDminusUp->Fit(gaussFit, "Q", "R", 0, 12);
  const double dminusRadUp = gaussFit->GetParameter(0);
  
  catsDminusCoulombOnly.SetAnaSource(0, rCorepDminusDefault);
  catsDminusCoulombOnly.KillTheCat();
  catsDstarminus.SetAnaSource(0, rCorepDminusDefault);
  catsDstarminus.KillTheCat();
  catsDminusHaidenbauer.SetAnaSource(0, rCorepDminusDefault);
  catsDminusHaidenbauer.KillTheCat();

  FillSourceGraph(catsDminusCoulombOnly, grSourceDminus);
  FillCkGraph(catsDminusCoulombOnly, grDminusCoulomb);
  FillCkGraph(catsDminusHaidenbauer, grDminusHaidenbauer);
  
  auto cminus = new TCanvas("p-D- source", "p-D- source");
  grSourceDminus->Draw("AL3");
  grSourceDminus->GetXaxis()->SetRangeUser(0, 10);
  grSourceDminus->Fit(gaussFit, "Q", "R", 0, 12);
  const double dminusRad = gaussFit->GetParameter(0);
  double err = 0.5 * (dminusRadUp-dminusRadLow);
  auto leg = new TLegend(0.5, 0.72, 0.8, 0.87);
  leg->SetTextFont(42);
  leg->SetTextSize(0.9 * gStyle->GetTextSize());
  leg->AddEntry(grSourceDminus, Form("#it{r}_{core} = (%.2f #pm 0.08) fm", rCorepDminusDefault), "l");
  leg->AddEntry(gaussFit, Form("#it{r}_{eff} = (%.2f #pm %.2f) fm", dminusRad, err));
  leg->Draw("same");
  cminus->Print("Source_Dminus.pdf");
  std::cout << "p-D- radius " << dminusRad << "" << dminusRadLow << " " << dminusRadUp << " fm\n";

  auto cCompareDminus = new TCanvas("Compare", "Compare");
  grDminusCoulomb->Draw("AL3");
  grDminusCoulomb->GetXaxis()->SetRangeUser(0, 300);
  grDminusCoulomb->GetYaxis()->SetRangeUser(0.6, 1.5);
  grDminusHaidenbauer->Draw("L3");
  auto legComp = new TLegend(0.45, 0.735, 0.45 + 0.3, 0.9);
  legComp->SetBorderSize(0);
  legComp->SetTextFont(42);
  legComp->SetHeader(
      TString::Format(
          "p#minus#kern[-0.95]{ }D^{#minus}, #it{r}_{eff} = %.2f fm",
          dminusRad));
  legComp->SetTextSize(gStyle->GetTextSize() * 0.9);
  legComp->AddEntry(grDminusCoulomb, "Coulomb-only", "l");
  legComp->AddEntry(grDminusHaidenbauer, "Coulomb+strong", "l");
  legComp->Draw("same");
  cCompareDminus->Print("Dminus_model.pdf");

  auto DLM_pDminusCoulomb = new DLM_Ck(1, 0, catsDminusCoulombOnly);
  auto DLM_pDstarminus = new DLM_Ck(1, 0, catsDstarminus);

  DLM_CkDecomposition CkDec_pDminusCoulomb("pDminusTotal", 2,
                                          *DLM_pDminusCoulomb,
                                          histMomentumResolution);
  DLM_CkDecomposition CkDec_pDstarminus("pDstarminus", 0, *DLM_pDstarminus,
                                       nullptr);
  CkDec_pDminusCoulomb.AddContribution(0, pDstarDminus, DLM_CkDecomposition::cFeedDown,
                                      &CkDec_pDstarminus,
                                      histDecayKindematicsDstar);
  CkDec_pDminusCoulomb.AddContribution(1, flatDminus,
                                      DLM_CkDecomposition::cFake);
  CkDec_pDminusCoulomb.Update();
  FillCkGraph(DLM_pDminusCoulomb, CkDec_pDminusCoulomb, grDminusTotal);

  DLM_CkDecomposition CkDec_pDminusLambda("pDminusCoulomb", 1,
                                          *DLM_pDminusCoulomb,
                                          histMomentumResolution);
  CkDec_pDminusLambda.AddContribution(0, 1.f - primaryDminus,
                                      DLM_CkDecomposition::cFake);
  CkDec_pDminusLambda.Update();
  FillCkGraph(DLM_pDminusCoulomb, CkDec_pDminusLambda, grDminusLambda);
  
  DLM_CkDecomposition CkDec_pDstarminus_Smeared("pDstarminus_smear", 0,
                                               *DLM_pDstarminus,
                                               histDecayKindematicsDstar);
  DLM_CkDecomposition CkDec_pDstarminus_Lambda("pDstarminus_lambda", 1,
                                              *DLM_pDstarminus,
                                              histDecayKindematicsDstar);
  CkDec_pDstarminus_Lambda.AddContribution(0, 1.f - pDstarDminus,
                                          DLM_CkDecomposition::cFake);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus_Smeared, grDstarminusSmeared);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus, grDstarminusGenuine);
  FillCkGraph(DLM_pDstarminus, CkDec_pDstarminus_Lambda, grDstarminusLambda);

  auto dminusLambda = new TCanvas("p-D- Experimental", "p-D- Experimental");
  grDminusTotal->Draw("AL3");
  auto legDminus2 = new TLegend(0.45, 0.9-0.125, 0.45 + 0.3, 0.9);
  legDminus2->SetBorderSize(0);
  legDminus2->SetTextFont(42);
  legDminus2->SetHeader(
      TString::Format("p#minus#kern[-0.95]{ }D^{#minus}, #it{r}_{eff} = %.2f fm",
                      dminusRad));
  legDminus2->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDminus2->AddEntry(grDminusTotal,
                      "p#minus#kern[-0.95]{ }D^{#minus} (Coulomb only)", "l");
  legDminus2->Draw("same");
  dminusLambda->Print("dminusExp.pdf");

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
          pDstarDminus * 100.),
      "l");
  legDstarminus->Draw("same");
  dDstar->Print("Transformation_Dstar_Dminus.pdf");

  auto cDminusTotal = new TCanvas("p-D- total", "p-D- total");
  grDminusTotal->Draw("AL3");
  grDminusTotal->GetXaxis()->SetNdivisions(504);
  grDminusTotal->GetYaxis()->SetNdivisions(504);
  grDminusTotal->SetLineStyle(0);
  grDminusTotal->GetXaxis()->SetRangeUser(0, 300);
  grDminusTotal->GetYaxis()->SetRangeUser(0.9, 2.5);
  grDminusLambda->Draw("L3");
  grDstarminusLambda->Draw("L3");
  auto legDminus = new TLegend(0.45, 0.6, 0.45 + 0.3, 0.9);
  legDminus->SetBorderSize(0);
  legDminus->SetTextFont(42);
  legDminus->SetHeader(TString::Format("p#minus#kern[-0.95]{ }D^{#minus}, #it{r}_{eff} = %.2f fm", dminusRad));
  legDminus->SetTextSize(gStyle->GetTextSize() * 0.9);
  legDminus->AddEntry(grDminusTotal, "Total", "l");
  legDminus->AddEntry(
      grDminusLambda,
      TString::Format("Genuine p#minus#kern[-0.95]{ }D^{#minus} (#lambda = %.2f %%)", primaryDminus * 100.f),
      "l");
  legDminus->AddEntry(
      grDstarminusLambda,
      TString::Format(
          "p#minus#kern[-0.95]{ }D*^{#minus} #rightarrow p#minus#kern[-0.95]{ }D^{#minus} (#lambda = %.2f %%)",
          pDstarDminus * 100.),
      "l");
  legDminus->Draw("same");
  cDminusTotal->Print("p-Dminus_total.pdf");

  auto dminus = new TCanvas("p-D- Coulomb", "p-D- Coulomb");
  grDminusCoulomb->Draw("AL3");
  grDminusLambda->Draw("L3");
  dminus->Print("dminusCoul.pdf");
  
  /// -----------------------------------------------------------------------------------
  /// Feeding from beauty
  
  auto grBpGenuine = new TGraph();
  DreamPlot::SetStyleGraph(grBpGenuine, 20, kRed + 3, 0.8);
  grBpGenuine->SetLineWidth(2);
  grBpGenuine->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grBpSmeared = new TGraph();
  DreamPlot::SetStyleGraph(grBpSmeared, 20, kRed + 3, 0.8);
  grBpSmeared->SetLineWidth(2);
  grBpSmeared->SetLineStyle(3);
  grBpSmeared->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grBpLambda = new TGraph();
  DreamPlot::SetStyleGraph(grBpLambda, 20, kRed + 2, 0.8);
  grBpLambda->SetLineWidth(2);
  grBpLambda->SetLineStyle(2);
  grBpLambda->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  
  const double protonMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double bminusMass =
      TDatabasePDG::Instance()->GetParticle(-521)->Mass() * 1000;

  CATS beautyCats;
  CATSparameters* cPars = nullptr;
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, 1.);
  beautyCats.SetAnaSource(GaussSource, *cPars);
  beautyCats.SetUseAnalyticSource(true);
  beautyCats.SetThetaDependentSource(false);
  beautyCats.SetExcludeFailedBins(false);
  beautyCats.SetMomBins(150, 0, 1500);
  beautyCats.SetNumChannels(1);
  beautyCats.SetNumPW(0, 1);
  beautyCats.SetSpin(0, 0);
  beautyCats.SetChannelWeight(0, 1.);
  beautyCats.SetQ1Q2(-1);
  beautyCats.SetPdgId(2212, 521);
  beautyCats.SetRedMass((protonMass * bminusMass) / (protonMass + bminusMass));
  beautyCats.KillTheCat();

  auto DLM_pB = new DLM_Ck(1, 0, beautyCats);
  DLM_CkDecomposition CkDec_pBGenuine("pBgen", 0, *DLM_pB, nullptr);
  CkDec_pBGenuine.Update();
  
  DLM_CkDecomposition CkDec_pB("pB", 1, *DLM_pB,
                                histDecayKindematicsBeauty);
  CkDec_pB.AddContribution(0, 1.f - Bfeeddown, DLM_CkDecomposition::cFake);
  CkDec_pB.Update();
  FillCkGraph(DLM_pB, CkDec_pB, grBpLambda);
  
  DLM_CkDecomposition CkDec_pB_Smeared("pB_smear", 0, *DLM_pB, histDecayKindematicsBeauty);
  CkDec_pB_Smeared.Update();
  FillCkGraph(DLM_pB, CkDec_pB_Smeared, grBpSmeared);
  FillCkGraph(DLM_pB, CkDec_pBGenuine, grBpGenuine);

  auto dB = new TCanvas("p-B- Transformation", "p-B- Transformation");
  grBpGenuine->Draw("AL3");
  grBpGenuine->GetXaxis()->SetRangeUser(0, 1000);
  grBpGenuine->GetYaxis()->SetRangeUser(0.8, 2);
  grBpLambda->Draw("L3");
  grBpSmeared->Draw("L3");
  auto legB = new TLegend(0.45, 0.55, 0.45 + 0.3, 0.85);
  legB->SetBorderSize(0);
  legB->SetTextFont(42);
  legB->SetHeader("p#minus#kern[-0.95]{ }B^{-} Coulomb");
  legB->SetTextSize(gStyle->GetTextSize() * 0.9);
  legB->AddEntry(grBpGenuine, "Genuine", "l");
  legB->AddEntry(grBpSmeared,
      "p#minus#kern[-0.95]{ }B^{-} #rightarrow p#minus#kern[-0.95]{ }D^{-}",
      "l");
  legB->AddEntry(grBpLambda,
      TString::Format(
          "p#minus#kern[-0.95]{ }B^{-} #rightarrow p#minus#kern[-0.95]{ }D^{-} (#lambda = %.2f %%)",
          (Bfeeddown) * 100.),
      "l");
  legB->Draw("same");
  dB->Print("Transformation_Bminus_Dminus.pdf");

  app->Run();
}
