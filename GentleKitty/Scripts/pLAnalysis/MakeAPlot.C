#include "TFile.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "DreamData.h"
#include "DreamPlot.h"

int main(int argc, char* argv[]) {
  DreamPlot::SetStyle();
  TFile* inFileNLO = TFile::Open(argv[1], "read");
  TFile* inFileLO = TFile::Open(argv[2], "read");
  TH1F* data = (TH1F*) inFileNLO->Get("Data_1");
  TF1* sysErr = (TF1*) inFileNLO->Get("Systematic_1");
  TGraphErrors* modelNLO = (TGraphErrors*) inFileNLO->Get("Model_0_imT_1");
  TGraphErrors* modelLO = (TGraphErrors*) inFileLO->Get("Model_0_imT_1");

  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("p#minus#kern[0.2]{#Lambda} #oplus #bar{p}#minus#kern[0.2]{#bar{#Lambda}}");
  LegNames.push_back("#chi EFT NLO (fit)");
  LegNames.push_back("#chi EFT LO (fit)");

  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  LegOptions.push_back("l");
  LegOptions.push_back("l");

  c1 = new TCanvas("c2", "c2", 0, 0, 500, 800);
  int counter = 1;
  TFile* out = TFile::Open("tmp.root", "recreate");
  TPad* pad;
  float LatexX = 0.;
  c1 = new TCanvas(Form("c_%u", counter), Form("c_%u", counter), 0, 0, 650,
                   650);
  pad = (TPad*) c1->cd(counter);
  pad->SetRightMargin(0.025);
  pad->SetTopMargin(0.025);
  pad->SetBottomMargin(0.12);
  pad->Draw();
  pad->cd();
  DreamData *Data = new DreamData(Form("Data_%i", counter));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
  Data->SetCorrelationFunction(data);
  Data->SetSystematics(sysErr, 2);

  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(true);
  Data->SetRangePlotting(0, 230, 0.9, data->GetMaximum() * 1.08);  //ranges
  Data->SetNDivisions(505);
  double lineWidth = 3;
  Data->FemtoModelFitBands(modelNLO, 1, 1, lineWidth, -1, true);  //Model colors
  Data->FemtoModelFitBands(modelLO, 3, 1, lineWidth, -1, true);  //Model colors
  float legXmin = 0.33;
  float fTextXMin = 0.35;
  Data->SetLegendCoordinates(legXmin, 0.685 - 0.09 * Data->GetNumberOfModels(),
                             legXmin + 0.4, 0.72);
  Data->DrawCorrelationPlot(pad);
  pad->cd();
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize() * .85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(fTextXMin, 0.91,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(
      fTextXMin,
      0.85,
      "High-mult. (0#kern[-0.65]{ }#minus#kern[-0.65]{ }0.17#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");

  TLatex text;
  text.SetNDC();
  text.SetTextColor(1);
  text.SetTextSize(gStyle->GetTextSize() * 0.85);
  text.DrawLatex(
      fTextXMin, 0.79,
      TString::Format("m_{T} #in [%.2f, %.2f) GeV/#it{c}^{2}", 1.26, 1.32));
  text.DrawLatex(fTextXMin, 0.73, TString::Format("%s", "Gaussian Source"));
  out->cd();
  c1->Write();
  c1->SaveAs(Form("mTBin_%u.pdf", counter));
  counter++;
  out->Close();
  return 0;
}
