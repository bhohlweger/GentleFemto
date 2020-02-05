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
  const char* filename = argv[1];//root file where TTree from sysfit is saved (should be original binning of default data)
  const char* SystFile = argv[2];//systematic root file where err.relative fit is stored from sys data (!= from original binning)
  const char* model = argv[3];
  int selector;
  TString convmodel = model;
  if(convmodel=="Haidenbauer") selector=0;
  if(convmodel=="Lednicky") selector=1;
  if(convmodel=="Coulomb") selector=2;

  std::cout<<"selector = "<<selector<< "----" << "Model = " << convmodel <<std::endl;

  TFile* systFile = TFile::Open(SystFile, "read");
  if (!systFile) {
    std::cout << "no syst file " << std::endl;
    return 0;
  }

  TF1* systematic = (TF1*) systFile->Get("SystError");

  VariationAnalysispAp* analysis = new VariationAnalysispAp("hCk_ReweightedpApVar");

  TCut chiSqCut = "chiSqNDF<30";
  // TCut polOrdCut = "iPolOrd==3";
  analysis->AppendAndCut(chiSqCut);
  // analysis->AppendAndCut(polOrdCut);
  analysis->ReadFitFile(filename);

  DreamData *ProtonAntiProton = new DreamData("ProtonAntiProton");
  ProtonAntiProton->SetUnitConversionData(1);
  ProtonAntiProton->SetUnitConversionCATS(1);
  ProtonAntiProton->SetCorrelationFunction(analysis->GetCorrelationFunction(0));
  ProtonAntiProton->SetSystematics(systematic, 2);
  if(selector==0) ProtonAntiProton->FemtoModelFitBands(analysis->GetModel(), 1, 1, 0., 0.45, true);
  if(selector==1) ProtonAntiProton->FemtoModelFitBands(analysis->GetModel(), 2, 1, 0., 0.45, true);
  if(selector==2) ProtonAntiProton->FemtoModelFitBands(analysis->GetModel(), 3, 1, 0., 0.45, true);

  ProtonAntiProton->FemtoModelDeviations(analysis->GetDeviationByBin(), 2);


  TCanvas* c_PAP = new TCanvas("CFpAp", "CFpAp", 0, 0, 650, 650);
  DreamPlot::SetStyle();
  c_PAP->cd();
  TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
  p1->SetRightMargin(0.025);
  p1->SetTopMargin(0.025);
  p1->SetBottomMargin(0.12);
  p1->Draw();

  ProtonAntiProton->SetLegendName("p-#bar{p}", "fpe");
  if(selector==0){
	  ProtonAntiProton->SetLegendName("Coulomb + #chi EFT (fit)", "l");
  }else if(selector==2){
 	  ProtonAntiProton->SetLegendName("Coulomb", "l");
  }else if(selector==1){
      ProtonAntiProton->SetLegendName("Coulomb + Lednicky-Lyuboshits (fit)", "l");
  }
  ProtonAntiProton->SetRangePlotting(0, 300, 0.725, 3.);
  ProtonAntiProton->SetNDivisions(505);
  ProtonAntiProton->SetLegendCoordinates(
      0.30, 0.65 - 0.09 * ProtonAntiProton->GetNumberOfModels(), 0.7, 0.725);
  ProtonAntiProton->DrawCorrelationPlot(p1,0,kGray+2,0.7,0.8);
  p1->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .55);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.32, 0.91, Form("#bf{ALICE}"));
  BeamText.DrawLatex(0.32, 0.85,
                     Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(0.32, 0.79, "High Mult. (0-0.17% INEL > 0)");
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
  TFile* out = TFile::Open(Form("tmp_pAp_%s.root",model), "recreate");
  out->cd();
  c_PAP->Write();
  c_PAP->SaveAs(Form("CF_pAp_%s.pdf", model));

  out->Write();
  out->Close();
  systFile->Close();

delete analysis;
delete ProtonAntiProton;
delete c_PAP;

  return 0;
}
