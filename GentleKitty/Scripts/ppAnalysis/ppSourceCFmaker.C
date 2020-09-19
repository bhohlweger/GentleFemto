#include "TApplication.h"
#include "TidyCats.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "DLM_CkDecomposition.h"
#include "TGraph.h"

void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

void FillSourceGraph(CATS &kitty, TGraph *gr) {
  for (double i = 0; i < 150; ++i) {
    gr->SetPoint(i, i * 0.1, kitty.EvaluateTheSource(0, i * 0.1, 0));
  }
}

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();

  TApplication *app = new TApplication("app", 0, 0);

  const float radiusGaussian = 1.25;
  const float radiusLevy = radiusGaussian;
  const float radiusCauchy = radiusGaussian;

  const int momBins = 500;
  const double kMin = 0;
  const double kMax = 2.f * momBins;

  TidyCats *tidy = new TidyCats();

  auto grGaussian = new TGraph();
  grGaussian->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  grGaussian->SetName("gauss");
  DreamPlot::SetStyleGraph(grGaussian, 20, kBlue+2);
  auto grLevy = new TGraph();
  grLevy->SetName("levy");
  grLevy->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grLevy, 20, kOrange+2);
  auto grCauchy = new TGraph();
  grCauchy->SetName("cauchy");
  grCauchy->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleGraph(grCauchy, 20, kGreen+2);

  auto grGaussianSource = new TGraph();
  grGaussianSource->SetName("gaussSource");
  grGaussianSource->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");
  DreamPlot::SetStyleGraph(grGaussianSource, 20, kBlue+3);
  auto grLevySource = new TGraph();
  grLevySource->SetName("levySource");
  grLevySource->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");
  DreamPlot::SetStyleGraph(grLevySource, 20, kCyan+1);
  auto grCauchySource = new TGraph();
  grCauchySource->SetName("cauchySource");
  grCauchySource->SetTitle(";#it{r} (fm); 4#pi#it{r}^{2} S(#it{r}) (fm^{-1})");
  DreamPlot::SetStyleGraph(grCauchySource, 20, kGreen+3);

  CATS ppGaussian;
  tidy->GetCatsProtonProton(&ppGaussian, momBins, kMin, kMax, TidyCats::sGaussian);
  ppGaussian.KillTheCat();
  auto CFppGaussian = new DLM_Ck(1, 0, ppGaussian);
  CFppGaussian->SetSourcePar(0, radiusGaussian);
  CFppGaussian->Update();
  FillWaveGraph(ppGaussian, grGaussian);
  FillSourceGraph(ppGaussian, grGaussianSource);

  CATS ppLevy;
  tidy->GetCatsProtonProton(&ppLevy, momBins, kMin, kMax, TidyCats::sLevy);
  ppLevy.KillTheCat();
  auto CFppLevy = new DLM_Ck(2, 0, ppLevy);
  CFppLevy->SetSourcePar(0, radiusLevy);
  CFppLevy->SetSourcePar(1, 1.5);
  CFppLevy->Update();
  FillWaveGraph(ppLevy, grLevy);
  FillSourceGraph(ppLevy, grLevySource);

  CATS ppCauchy;
  tidy->GetCatsProtonProton(&ppCauchy, momBins, kMin, kMax,
                            TidyCats::sCauchy);
  ppCauchy.KillTheCat();
  auto CFppCauchy = new DLM_Ck(1, 0, ppCauchy);
  CFppCauchy->SetSourcePar(0, radiusCauchy);
  CFppCauchy->Update();
  FillWaveGraph(ppCauchy, grCauchy);
  FillSourceGraph(ppCauchy, grCauchySource);


  auto cGaussian = new TCanvas();
  cGaussian->SetTopMargin(0.025);
  cGaussian->SetRightMargin(0.015);
  grGaussian->Draw("AL");
  grLevy->Draw("Lsame");
  grCauchy->Draw("Lsame");

  auto cGaussianSource = new TCanvas();
  grGaussianSource->Draw("AL");
  grGaussianSource->GetXaxis()->SetNdivisions(505);
  grGaussianSource->GetYaxis()->SetNdivisions(505);
  grGaussianSource->GetXaxis()->SetRangeUser(0, 10);
  grGaussianSource->GetYaxis()->SetRangeUser(0, 0.4);
  grLevySource->Draw("L same");
  grCauchySource->Draw("L same");

  auto leg2 = new TLegend(0.6, 0.625, 0.9, 0.925);
  leg2->SetTextFont(43);
  leg2->SetTextSize(22);
  leg2->SetHeader(Form("#it{r}_{source} = %.2f fm", radiusGaussian));
  leg2->AddEntry(grCauchySource, "#alpha = 1", "l");
  leg2->AddEntry(grLevySource, "#alpha = 1.5", "l");
  leg2->AddEntry(grGaussianSource, "#alpha = 2", "l");
  leg2->Draw("same");

  cGaussianSource->Print("SourceShapes.pdf");

  app->Run();

  return 0;
}
