#include "VariationAnalysis.h"
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
  TFile* systFile = TFile::Open(SystFile, "read");
  if (!systFile) {
    std::cout << "no syst file " << std::endl;
    return 0;
  }
  TF1* systematic = (TF1*) systFile->Get("SystError");
//  VariationAnalysis* analysis = new VariationAnalysis("hCk_FixShiftedppVar", 26,
//   81);
  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar");
  //  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar", 26,
//                                                      54);
  TCut chiSqCut = "chiSqNDF<30";
  analysis->AppendAndCut(chiSqCut);
//  TCut PolCut = "PolBL==1";
//  analysis->AppendAndCut(PolCut);

  analysis->ReadFitFile(filename);
  analysis->EvalRadius();
  std::cout << "Radius Mean: " << analysis->GetRadMean() << " Radius stat Err: "
            << analysis->GetRadStatErr() << " Radius SystErr Down: "
            << analysis->GetRadSystDown() << " Radius Syst Err Up: "
            << analysis->GetRadSystUp() << std::endl;
  DreamData *ProtonProton = new DreamData("ProtonProton");
  ProtonProton->SetUnitConversionData(1);
  ProtonProton->SetUnitConversionCATS(1);
  ProtonProton->SetCorrelationFunction(analysis->GetCorrelationFunction(0));
  ProtonProton->SetSystematics(systematic, 2);
  ProtonProton->FemtoModelFitBands(analysis->GetModel(), 2, 1, 3, -3000, true);
  ProtonProton->FemtoModelDeviations(analysis->GetDeviationByBin(), 2);

  TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 650);
//  TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 687.5);
  DreamPlot::SetStyle();
  c_PP->cd();
  TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
//  TPad *p1 = new TPad("p1", "p1", 0., 0.3, 1., 1.);
  p1->SetRightMargin(0.025);
  p1->SetTopMargin(0.025);
//  p1->SetBottomMargin(0.0);
  p1->SetBottomMargin(0.12);
  p1->Draw();
  ProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p}", "fpe");
  ProtonProton->SetLegendName("Coulomb + Argonne #nu_{18} (fit)", "l");
  ProtonProton->SetRangePlotting(0, 200, 0.725, 3.5);
  ProtonProton->SetInletRangePlotting(50, 375, 0.94, 1.06);
  ProtonProton->SetInletCoordinates(0.27, 0.23, 0.95, 0.55);
//  ProtonProton->SetInletCoordinates(0.27, 0.15, 0.95, 0.55);
  ProtonProton->SetNDivisions(505);
  ProtonProton->SetLegendCoordinates(
      0.30, 0.65 - 0.09 * ProtonProton->GetNumberOfModels(), 0.7, 0.725);
  ProtonProton->DrawCorrelationPlot(p1);
  p1->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .85);
  BeamText.SetNDC(kTRUE);
//    BeamText.DrawLatex(0.5, 0.875, "ALICE");

  BeamText.DrawLatex(0.32, 0.91, Form("#bf{ALICE}"));
  BeamText.DrawLatex(0.32, 0.85,
                     Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.32, 0.79, "High Mult. (0-0.17% INEL > 0)");
  text.SetNDC();
  text.SetTextColor(1);
  text.SetTextSize(gStyle->GetTextSize() * 0.85);
  text.DrawLatex(0.32, 0.73, "Gaussian + Resonance source");
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
  TFile* out = TFile::Open("tmp.root", "update");
  out->cd();
  c_PP->Write();
  c_PP->SaveAs("CF_Gauss.pdf");
  systFile->Close();
  out->Write();
  out->Close();
  return 0;
}
