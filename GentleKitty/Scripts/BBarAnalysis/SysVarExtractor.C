#include "VariationAnalysispAp.h"
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
  printf("--- debug 1 ---\n");
  TFile* systFile = TFile::Open(SystFile, "read");
  if (!systFile) {
    std::cout << "no syst file " << std::endl;
    return 0;
  }
  printf("--- debug 2 ---\n");

  TF1* systematic = (TF1*) systFile->Get("SystError");
//  VariationAnalysis* analysis = new VariationAnalysis("hCk_FixShiftedppVar", 26,
//   81);
printf("--- debug 3 ---\n");

  VariationAnalysispAp* analysis = new VariationAnalysispAp("hCk_ReweightedpApVar");
  //  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar", 26,
//                                                      54);
  TCut chiSqCut = "chiSqNDF<30";
  analysis->AppendAndCut(chiSqCut);
//  TCut PolCut = "PolBL==1";
//  analysis->AppendAndCut(PolCut);
printf("--- debug 4 ---\n");
  analysis->ReadFitFile(filename);
  printf("--- debug 5 ---\n");


  DreamData *ProtonAntiProton = new DreamData("ProtonAntiProton");
  ProtonAntiProton->SetUnitConversionData(1);
  ProtonAntiProton->SetUnitConversionCATS(1);
  ProtonAntiProton->SetCorrelationFunction(analysis->GetCorrelationFunction(0));
  ProtonAntiProton->SetSystematics(systematic, 2);
  ProtonAntiProton->FemtoModelFitBands(analysis->GetModel(), 2, 1, 3, -3000, true);
  ProtonAntiProton->FemtoModelDeviations(analysis->GetDeviationByBin(), 2);
  printf("--- debug 6 ---\n");

  TCanvas* c_PAP = new TCanvas("CFpAp", "CFpAp", 0, 0, 650, 650);
//  TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 687.5);
  DreamPlot::SetStyle();
  c_PAP->cd();
  TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
//  TPad *p1 = new TPad("p1", "p1", 0., 0.3, 1., 1.);
  p1->SetRightMargin(0.025);
  p1->SetTopMargin(0.025);
//  p1->SetBottomMargin(0.0);
  p1->SetBottomMargin(0.12);
  p1->Draw();
  printf("--- debug 7 ---\n");

  ProtonAntiProton->SetLegendName("p-#bar{p}", "fpe");
  ProtonAntiProton->SetLegendName("Coulomb + #chi EFT (fit)", "l");
  ProtonAntiProton->SetRangePlotting(0, 300, 0.725, 3.);
  ProtonAntiProton->SetNDivisions(505);
  ProtonAntiProton->SetLegendCoordinates(
      0.30, 0.65 - 0.09 * ProtonAntiProton->GetNumberOfModels(), 0.7, 0.725);
  // void DreamData::DrawCorrelationPlot(TPad* c, const int color,
  //                                     const int systematicsColor,
  //                                     const float legendTextScale, const float markersize)
  ProtonAntiProton->DrawCorrelationPlot(p1,0,kGray+2,0.7,0.8);
  p1->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .55);
  BeamText.SetNDC(kTRUE);
//    BeamText.DrawLatex(0.5, 0.875, "ALICE");

  BeamText.DrawLatex(0.32, 0.91, Form("#bf{ALICE}"));
  BeamText.DrawLatex(0.32, 0.85,
                     Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.32, 0.79, "High Mult. (0-0.072% INEL)");
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
  TFile* out = TFile::Open("tmp.root", "update");
  out->cd();
  c_PAP->Write();
  c_PAP->SaveAs("CF_pAp_Haide.pdf");
  printf("--- debug 9---\n");

  systFile->Close();
  printf("--- debug 10---\n");

  out->Write();
  printf("--- debug 11---\n");

  out->Close();
  printf("--- debug 12---\n");

  return 0;
}
