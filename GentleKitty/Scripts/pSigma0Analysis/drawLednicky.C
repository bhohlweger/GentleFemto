#include "DLM_CkModels.h"
#include "DLM_CkDecomposition.h"
#include "TApplication.h"
#include "TMath.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  TApplication *app = new TApplication("app", 0, 0);

  DreamPlot::SetStyle();

  const int momBins = 100;
  const double kmin = 0;
  const double kmax = 300;
  const double radius1 = 1.25;
  const double radius2 = 4;
  const double scatLen1 = 1.5;
  const double effRange = 1;
  const double scatLen2 = -1.5;
  const double scatLen3 = 0.8;

  const int lineWidth = 3;


  auto grCk11 = new TGraph();
  DreamPlot::SetStyleGraph(grCk11, 20, kBlue + 3, 0.8);
  grCk11->SetLineWidth(lineWidth);

  auto grCk12 = new TGraph();
  DreamPlot::SetStyleGraph(grCk12, 20, kBlue + 3, 0.8);
  grCk12->SetLineWidth(lineWidth);
  grCk12->SetLineStyle(2);

  auto grCk21 = new TGraph();
  DreamPlot::SetStyleGraph(grCk21, 20, kGreen + 3, 0.8);
  grCk21->SetLineWidth(lineWidth);

  auto grCk22 = new TGraph();
  DreamPlot::SetStyleGraph(grCk22, 20, kGreen + 3, 0.8);
  grCk22->SetLineWidth(lineWidth);
  grCk22->SetLineStyle(2);

  auto grCk31 = new TGraph();
  DreamPlot::SetStyleGraph(grCk31, 20, kOrange + 2, 0.8);
  grCk31->SetLineWidth(lineWidth);

  auto grCk32 = new TGraph();
  DreamPlot::SetStyleGraph(grCk32, 20, kOrange + 2, 0.8);
  grCk32->SetLineWidth(lineWidth);
  grCk32->SetLineStyle(2);

  auto Ck11 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck11->SetPotPar(0, scatLen1);
  Ck11->SetPotPar(1, effRange);
  Ck11->SetSourcePar(0, radius1);
  Ck11->Update();
  FillCkGraph(Ck11, grCk11);

  auto Ck12 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck12->SetPotPar(0, scatLen1);
  Ck12->SetPotPar(1, effRange);
  Ck12->SetSourcePar(0, radius2);
  Ck12->Update();
  FillCkGraph(Ck12, grCk12);

  auto Ck21 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck21->SetPotPar(0, scatLen2);
  Ck21->SetPotPar(1, effRange);
  Ck21->SetSourcePar(0, radius1);
  Ck21->Update();
  FillCkGraph(Ck21, grCk21);

  auto Ck22 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck22->SetPotPar(0, scatLen2);
  Ck22->SetPotPar(1, effRange);
  Ck22->SetSourcePar(0, radius2);
  Ck22->Update();
  FillCkGraph(Ck22, grCk22);

  auto Ck31 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck31->SetPotPar(0, scatLen3);
  Ck31->SetPotPar(1, effRange);
  Ck31->SetSourcePar(0, radius1);
  Ck31->Update();
  FillCkGraph(Ck31, grCk31);

  auto Ck32 = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet);
  Ck32->SetPotPar(0, scatLen3);
  Ck32->SetPotPar(1, effRange);
  Ck32->SetSourcePar(0, radius2);
  Ck32->Update();
  FillCkGraph(Ck32, grCk32);

  auto dummyHist = new TH1F("dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, kmax);
  DreamPlot::SetStyleHisto(dummyHist, 20, kWhite);
  dummyHist->GetYaxis()->SetRangeUser(0., 3);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);

  auto c = new TCanvas("CFComp1", "CFComp1"); //, 0, 0, 650, 550);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.025);
  dummyHist->Draw();
  grCk11->Draw("L3");
  grCk21->Draw("L3");
  grCk31->Draw("L3");
  grCk12->Draw("L3");
  grCk22->Draw("L3");
  grCk32->Draw("L3");

  TLatex text2;
  text2.SetNDC(true);
  text2.SetTextSize(gStyle->GetTextSize() * 0.9);
  text2.DrawLatex(0.45, 0.855, "#bf{Strong Interaction}");

  auto leg1 = new TLegend(0.45, 0.56, 0.65, 0.84);
  leg1->SetBorderSize(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg1->SetHeader(Form("#it{r}_{0} = %.2f fm", radius1));
  leg1->AddEntry(grCk11, Form("#it{f}_{0} = %.1f fm", scatLen1), "l");
  leg1->AddEntry(grCk31, Form("#it{f}_{0} = %.1f fm", scatLen3), "l");
  leg1->AddEntry(grCk21, Form("#it{f}_{0} = %.1f fm", scatLen2), "l");
  leg1->Draw("same");

  auto leg2 = new TLegend(0.7, 0.56, 0.9, 0.84);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg2->SetHeader(Form("#it{r}_{0} = %.2f fm", radius2));
  leg2->AddEntry(grCk12, Form("#it{f}_{0} = %.1f fm", scatLen1), "l");
  leg2->AddEntry(grCk32, Form("#it{f}_{0} = %.1f fm", scatLen3), "l");
  leg2->AddEntry(grCk22, Form("#it{f}_{0} = %.1f fm", scatLen2), "l");
  leg2->Draw("same");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextSize(gStyle->GetTextSize() * 0.9);
  text.SetNDC(true);
  text.DrawLatex(0.7, 0.25, Form("#it{d}_{0} = %.1f fm", effRange));

  c->Print("SomeLednicky.pdf");

  app->Run();

}
