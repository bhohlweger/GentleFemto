#include "VariationAnalysispAL.h"
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
  if(convmodel=="Lednicky") selector=1;
  else selector=0;
  TFile* systFile = TFile::Open(SystFile, "read");
  if (!systFile) {
    std::cout << "no syst file " << std::endl;
    return 0;
  }

  TF1* systematic = (TF1*) systFile->Get("SystError");


  VariationAnalysispAL* analysis = new VariationAnalysispAL("hCk_ReweightedpALVar");

  TCut chiSqCut = "chiSqNDF<30";
  analysis->AppendAndCut(chiSqCut);
  analysis->ReadFitFile(filename);

  DreamData *ProtonAntiLambda = new DreamData("ProtonAntiLambda");
  ProtonAntiLambda->SetUnitConversionData(1);
  ProtonAntiLambda->SetUnitConversionCATS(1);
  ProtonAntiLambda->SetCorrelationFunction(analysis->GetCorrelationFunction(0));
  ProtonAntiLambda->SetSystematics(systematic, 6);
  ProtonAntiLambda->FemtoModelFitBands(analysis->GetModel(), 2, 1, 3, -3000, true);
  ProtonAntiLambda->FemtoModelDeviations(analysis->GetDeviationByBin(), 2);

  TCanvas* c_PAL = new TCanvas("CFpAL", "CFpAL", 0, 0, 650, 650);
  DreamPlot::SetStyle();
  c_PAL->cd();
  TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
  p1->SetRightMargin(0.025);
  p1->SetTopMargin(0.025);
  p1->SetBottomMargin(0.12);
  p1->Draw();

  ProtonAntiLambda->SetLegendName("p-#bar{#Lambda} #oplus #bar{p}-#Lambda", "fpe");
  ProtonAntiLambda->SetLegendName("Lednicky-Lyuboshits (fit)", "l");
  ProtonAntiLambda->SetRangePlotting(0, 300, 0.6, 1.1);
  ProtonAntiLambda->SetNDivisions(505);
  ProtonAntiLambda->SetLegendCoordinates(
      0.40, 0.35 - 0.09 * ProtonAntiLambda->GetNumberOfModels(), 0.6, 0.425);
  // void DreamData::DrawCorrelationPlot(TPad* c, const int color,
  //                                     const int systematicsColor,
  //                                     const float legendTextScale, const float markersize)
  ProtonAntiLambda->DrawCorrelationPlot(p1,0,kGray+2,0.7,0.8);
  p1->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .55);
  BeamText.SetNDC(kTRUE);

  BeamText.DrawLatex(0.2, 0.91, Form("#bf{ALICE}"));
  BeamText.DrawLatex(0.2, 0.85,
                     Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.2, 0.79, "High Mult. (0-0.072% INEL)");
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
  TFile* out = TFile::Open(Form("tmp_pAl_%s.root",model), "recreate");
  out->cd();
  c_PAL->Write();
  c_PAL->SaveAs(Form("CF_pAL_%s.pdf", model));
  systFile->Close();
  out->Write();
  out->Close();

  return 0;
}
