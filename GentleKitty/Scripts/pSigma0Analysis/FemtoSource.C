#include "TApplication.h"
#include "DreamPlot.h"
#include "TDatabasePDG.h"
#include "TidyCats.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

/// =====================================================================================
void FillSourceGraph(CATS &kitty, TGraph *gr) {
  for (double i = 0; i < 100; ++i) {
    gr->SetPoint(i, i * 0.2, kitty.EvaluateTheSource(0, i * 0.2, 0));
  }
}

/// =====================================================================================
int main(int argc, char* argv[]) {
  TApplication *app = new TApplication("app", 0, 0);

  DreamPlot::SetStyle();

  const double rCore = 0.87;
  const double rCorepp = 1.17;

  CATS pSigmaResonantSource;
  TidyCats *tidySigma = new TidyCats();
  tidySigma->GetCatsProtonSigma0(&pSigmaResonantSource, 100, 0, 250, TidyCats::sResonance, TidyCats::pSigma0Haidenbauer);

  auto grSigmaSource = new TGraph();
  DreamPlot::SetStyleGraph(grSigmaSource, 20, kBlue + 3, 0.8);
  grSigmaSource->SetLineWidth(2);
  grSigmaSource->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");
  pSigmaResonantSource.SetAnaSource(0, rCore);  // this is the radius for this mT
  pSigmaResonantSource.KillTheCat();
  FillSourceGraph(pSigmaResonantSource, grSigmaSource);

  auto grSigmaSourceUp = new TGraph();
  pSigmaResonantSource.SetAnaSource(0, rCore+0.06);  // uncertainties
  pSigmaResonantSource.KillTheCat();
  FillSourceGraph(pSigmaResonantSource, grSigmaSourceUp);

  auto grSigmaSourceLow = new TGraph();
  pSigmaResonantSource.SetAnaSource(0, rCore-0.06);  // uncertainties
  pSigmaResonantSource.KillTheCat();
  FillSourceGraph(pSigmaResonantSource, grSigmaSourceLow);

  delete tidySigma;

  CATS ppResonantSource;
  TidyCats *tidyProton= new TidyCats();
  tidyProton->GetCatsProtonProton(&ppResonantSource, 100, 0, 250, TidyCats::sResonance);
  ppResonantSource.SetAnaSource(0, rCorepp);  // this is the radius for this mT
  ppResonantSource.KillTheCat();

  auto grPPSource = new TGraph();
  DreamPlot::SetStyleGraph(grPPSource, 20, kTeal + 2, 0.8);
  grPPSource->SetLineWidth(2);
  grPPSource->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");
  FillSourceGraph(ppResonantSource, grPPSource);
  delete tidyProton;


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Plotting

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

  grSigmaSource->Fit(gaussFit, "Q", "RQ", 0, 12);
  const double sigmaRadDefault = gaussFit->GetParameter(0);
  grSigmaSourceUp->Fit(gaussFit, "Q", "RQ", 0, 12);
  const double sigmaRadUp = gaussFit->GetParameter(0);
  grSigmaSourceLow->Fit(gaussFit, "Q", "RQ", 0, 12);
  const double sigmaRadLow = gaussFit->GetParameter(0);
  std::cout << "p-Sigma0 radii " << sigmaRadLow << " " << sigmaRadDefault << " " << sigmaRadUp << "\n";


  auto gaussFitPP = new TF1(
      "PPgaus",
      [&](double *x, double *p) {return
        4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
        std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
      },
      0, 10, 1);
  gaussFitPP->SetParameter(0, 1.25);
  gaussFitPP->SetNpx(1000);
  gaussFitPP->SetLineColorAlpha(kTeal + 2, 0.8);
  gaussFitPP->SetLineWidth(2);
  gaussFitPP->SetLineStyle(2);

  auto d = new TCanvas();
  d->SetTopMargin(0.025);
  d->SetRightMargin(0.025);
  grSigmaSource->Draw("AL");
  grSigmaSource->GetXaxis()->SetNdivisions(505);
  grSigmaSource->GetYaxis()->SetNdivisions(505);
  grSigmaSource->GetXaxis()->SetRangeUser(0, 10);
  grSigmaSource->GetYaxis()->SetRangeUser(0, 0.5);
  grSigmaSource->Fit(gaussFit, "", "R", 0, 12);

  grPPSource->Draw("L same");
  grPPSource->Fit(gaussFitPP, "", "R", 0, 12);
  auto leg2 = new TLegend(0.45, 0.625, 0.75, 0.925);
  leg2->SetTextFont(43);
  leg2->SetTextSize(22);
  leg2->AddEntry(grPPSource, Form("p#minus#kern[-0.95]{ }p, #it{r}_{core} = (%.2f #pm 0.03) fm", rCorepp), "l");
  leg2->AddEntry(gaussFitPP, Form("p#minus#kern[-0.95]{ }p, #it{r}_{eff} = (%.2f #pm 0.03) fm", gaussFitPP->GetParameter(0)) , "l");
  leg2->AddEntry(grSigmaSource,  Form("p#minus#kern[-0.95]{ }#Sigma^{0}, #it{r}_{core} = (%.2f #pm 0.06) fm", rCore), "l");
  leg2->AddEntry(gaussFit, Form("p#minus#kern[-0.95]{ }#Sigma^{0}, #it{r}_{eff} = (%.2f #pm 0.06) fm", gaussFit->GetParameter(0)), "l");
  leg2->Draw("same");
  d->Print("Source.pdf");

  app->Run();
}
