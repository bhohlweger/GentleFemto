#include "VariationAnalysisLAL.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
int main(int argc, char *argv[]) {
  const char* filename = argv[1];
  const char* SystFile = argv[2];
  const char* model = argv[3];
  int selector;
  TString convmodel = model;
  if(convmodel=="Haidenbauer") selector=0;
  if(convmodel=="Lednicky") selector=1;
  if(convmodel=="Coulomb") selector=2;

  TFile* systFile = TFile::Open(SystFile, "read");
  if (!systFile) {
    std::cout << "no syst file " << std::endl;
    return 0;
  }

  TF1* systematic = (TF1*) systFile->Get("SystError");
  VariationAnalysisLAL* analysis = new VariationAnalysisLAL("hCk_ReweightedLALVar");

  TCut chiSqCut = "chiSqNDF<30";
  analysis->AppendAndCut(chiSqCut);
  analysis->ReadFitFile(filename);


  DreamData *LambdaAntiLambda = new DreamData("LambdaAntiLambda");
  LambdaAntiLambda->SetUnitConversionData(1);
  LambdaAntiLambda->SetUnitConversionCATS(1);
  LambdaAntiLambda->SetCorrelationFunction(analysis->GetCorrelationFunction(0));
  LambdaAntiLambda->SetSystematics(systematic, 8);
  if(selector==0) LambdaAntiLambda->FemtoModelFitBands(analysis->GetModel(), 1, 1, 3., 0.45, true);
  if(selector==1) LambdaAntiLambda->FemtoModelFitBands(analysis->GetModel(), 2, 1, 3., 0.45, true);
  if(selector==2) LambdaAntiLambda->FemtoModelFitBands(analysis->GetModel(), 3, 1, 3., 0.45, true);
  LambdaAntiLambda->FemtoModelDeviations(analysis->GetDeviationByBin(), 8);

  TCanvas* c_LAL = new TCanvas("CFLAL", "CFLAL", 0, 0, 650, 650);
  DreamPlot::SetStyle();
  c_LAL->cd();
  TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
  p1->SetRightMargin(0.025);
  p1->SetTopMargin(0.025);
  p1->SetBottomMargin(0.12);
  p1->Draw();

  LambdaAntiLambda->SetLegendName("#Lambda-#bar{#Lambda}", "fpe");
  LambdaAntiLambda->SetLegendName("Lednicky-Lyuboshits (fit)", "l");
  LambdaAntiLambda->SetRangePlotting(0, 300, 0.6, 1.1);
  LambdaAntiLambda->SetNDivisions(505);
  LambdaAntiLambda->SetLegendCoordinates(
      0.40, 0.35 - 0.09 * LambdaAntiLambda->GetNumberOfModels(), 0.6, 0.425);
  // void DreamData::DrawCorrelationPlot(TPad* c, const int color,
  //                                     const int systematicsColor,
  //                                     const float legendTextScale, const float markersize)
  LambdaAntiLambda->DrawCorrelationPlot(p1,0,kGray+2,0.7,0.8);
  p1->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .55);
  BeamText.SetNDC(kTRUE);
//    BeamText.DrawLatex(0.5, 0.875, "ALICE");

  BeamText.DrawLatex(0.2, 0.91, Form("#bf{ALICE}"));
  BeamText.DrawLatex(0.2, 0.85,
                     Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.2, 0.79, "High Mult. (0-0.17% INEL > 0)");
  text.SetNDC();
  text.SetTextColor(1);
  text.SetTextSize(gStyle->GetTextSize() * 0.55);
  // text.DrawLatex(0.32, 0.73, "Gaussian + Resonance source");
//      Form("#it{r}_{#kern[-0.17]{core}} = %.3f#kern[-0.1]{#pm}%.3f(stat.)^{+%.3f}_{-%.3f}(syst.) fm",
//           analysis->GetRadMean(), analysis->GetRadStatErr(),
//           analysis->GetRadSystUp(), analysis->GetRadSystDown()));
//  c_PP->cd();
//  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., .3);
//  p2->SetRightMargin(0.025);
//  p2->SetTopMargin(0.);
//  p2->SetBottomMargin(0.3);
//  p2->Draw();
//  ProtonProton->DrawDeviationPerBin(p2);
  TFile* out = TFile::Open(Form("tmp_LAL_%s.root",model), "recreate");
  out->cd();
  c_LAL->Write();
  c_LAL->SaveAs(Form("CF_LAL_%s.pdf", model));
  systFile->Close();
  out->Write();
  out->Close();

  return 0;
}
