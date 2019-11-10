#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TApplication.h"

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  const char* WorkDir = argv[1];
  TApplication app("TheApp",&argc,argv);
  TString cfName = TString::Format("%s/debug_Var0.root", WorkDir).Data();
  TString sysDataName =
      TString::Format("%s/Systematics_pXi.root", WorkDir).Data();
  TString sysNormName = TString::Format("%s/Systematics_pXiNorm.root", WorkDir)
      .Data();
  TString sysLamName = TString::Format("%s/Systematics_pXiLam.root",
                                       WorkDir).Data();

  TFile* cfFile = TFile::Open(cfName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
      "InputCF_ResGami_woBL_ResGami_GenuineGami");
  if (!cf_default) {
    std::cout << "Default  data not found \n";
    cfFile->ls();
    return 0;
  }
  TGraphErrors* coulomb = (TGraphErrors*) cfFile->FindObjectAny("Coulomb");
  TGraphErrors* hal = (TGraphErrors*) cfFile->FindObjectAny("HalQCD");
  TGraphErrors* esc = (TGraphErrors*) cfFile->FindObjectAny("ESC");

  TFile* sysDataFile = TFile::Open(sysDataName, "read");
  TF1* systDataErr = (TF1*) sysDataFile->FindObjectAny("SystError");
  if (!systDataErr) {
    std::cout << "systDataErr not found \n";
    sysDataFile->ls();
    return 0;
  }
  TFile* sysNormFile = TFile::Open(sysNormName, "read");
  TH1F* systNormErr = (TH1F*) sysNormFile->FindObjectAny("SystErrRel");
  if (!systNormErr) {
    std::cout << "systNormErr not found \n";
    sysNormFile->ls();
    return 0;
  }
  TFile* sysLamFile = TFile::Open(sysLamName, "read");
  TF1* systLamErr = (TF1*) sysLamFile->FindObjectAny("SystError");
  if (!systLamErr) {
    std::cout << "systLamErr not found \n";
    sysLamFile->ls();
    return 0;
  }

  TH1F* SystError = (TH1F*)cf_default->Clone("Sytematics");
  SystError->Reset();

  for (int iBin = 1; iBin <= SystError->GetNbinsX(); ++iBin) {
    double kStar = cf_default->GetBinCenter(iBin);
    double Ck = cf_default->GetBinContent(iBin);
    double errSystData = systDataErr->Eval(kStar) * Ck;
    double errSystNorm = systNormErr->GetBinContent(iBin) * Ck;
    double errSystLam = systLamErr->Eval(kStar) * Ck;
    double totErr = TMath::Sqrt(
        errSystData * errSystData + errSystNorm * errSystNorm
            + errSystLam * errSystLam);
    SystError->SetBinContent(iBin, totErr);
  }

  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("p-#Xi^{-} #bf{ALICE} data");
  LegNames.push_back("Coulomb");
  LegNames.push_back("HAL-QCD");
  LegNames.push_back("ESC 16");
  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  LegOptions.push_back("f");
  LegOptions.push_back("f");
  LegOptions.push_back("f");

  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.025);
  pad->SetTopMargin(0.025);
  pad->SetBottomMargin(0.12);
  pad->Draw();
  pad->cd();
  float fXmin = 0;
  float fXmax = 300;
  float fTextXMin = 0.35;
  float ymaxL = 0.81;
  DreamData *Data = new DreamData(Form("Data"));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
  Data->SetCorrelationFunction(cf_default);
  Data->SetSystematics(SystError, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(true);
  Data->FemtoModelFitBands(coulomb, 12, 1, 3, -4000, true);
  Data->FemtoModelFitBands(hal, 10, 10, 0, -4000, true);
  Data->FemtoModelFitBands(esc, 11, 8, 0, -4000, true);
  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_default->GetMaximum() * 1.5);  //ranges
  Data->SetNDivisions(803);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.45,  0.6,  0.8, ymaxL+0.03);
  Data->DrawCorrelationPlot(pad);
  pad->cd();
  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize() * .85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(fTextXMin, 0.91,
                     Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
  BeamText.DrawLatex(fTextXMin, 0.85, "High-mult.");
  BeamText.DrawLatex(
      fTextXMin,
      0.79,
      "(0#kern[-0.95]{ }#minus#kern[-0.05]{ }0.072#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
  out->cd();
  c1->Write();
  c1->SaveAs(Form("WithSys.pdf"));
  SystError->Write();
  out->Close();
  app.Run();
}
