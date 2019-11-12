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
  cf_default->SetName("DefaultMeV");
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
    double errSystData = systDataErr->Eval(kStar);
    double errSystNorm = systNormErr->GetBinContent(iBin);
    double errSystLam = systLamErr->Eval(kStar);
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
  float xmin = 0.;
  float xmax = 300.;
  float ymin = 0.6;
  float ymax = 10.2;

  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = c1->DrawFrame(0,0.6,300,4.5);
  const char *  texPtY="#it{C}(#it{k}*)";
  const char *  texPtX="#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetNdivisions(803);
  TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.1);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.1);
  pad->SetBottomMargin(0.1);
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
  Data->SetDrawAxis(false);
  Data->FemtoModelFitBands(coulomb, kGreen+1, 1,  3, -4000, true, false);
  Data->FemtoModelFitBands(hal,     kOrange+1, 10, 0, -4000, true, false);
  Data->FemtoModelFitBands(esc,     11, 8,  0, -4000, true);
  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_default->GetMaximum() * 1.5);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.45,  0.6,  0.8, ymaxL+0.03);
  Data->DrawCorrelationPlot(pad);
  pad->cd();
  out->cd();
  c1->Write();
  c1->SaveAs(Form("WithSys.pdf"));
  SystError->Write();
  out->Close();
  app.Run();
}
