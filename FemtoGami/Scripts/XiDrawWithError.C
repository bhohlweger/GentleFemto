#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TApplication.h"

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  gStyle->SetEndErrorSize(5);
  const char* WorkDir = argv[1];
  TApplication app("TheApp", &argc, argv);
  TString cfName = TString::Format("%s/debug_Var0.root", WorkDir).Data();
  TString cfgraphName = TString::Format("%s/CFinput_Var0.root", WorkDir).Data();
  TString sysDataName =
      TString::Format("%s/Systematics_pXi.root", WorkDir).Data();
  TString sysNormName = TString::Format("%s/Systematics_pXiNorm.root", WorkDir)
      .Data();
  TString sysLamName = TString::Format("%s/Systematics_pXiLam.root", WorkDir)
      .Data();
  TString sysResName = TString::Format("%s/Systematics_pXiRes.root", WorkDir)
      .Data();

  TFile* cfFile = TFile::Open(cfName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
      "InputCF_ResGami_woBL_ResGami_GenuineGami");
  TFile* cfgraphFile = TFile::Open(cfgraphName, "read");
  TGraphAsymmErrors* cf_graph = (TGraphAsymmErrors*) cfgraphFile->FindObjectAny(
      "Graph_from_hCk_RebinnedpXiVar0_0MeV");
  if (!cf_default) {
    std::cout << "Default  data not found \n";
    cfFile->ls();
    return 0;
  }
  cf_default->SetName("DefaultMeV");
  if (!cf_graph) {
    std::cout << "Default  graph not found \n";
    cfFile->ls();
    return 0;
  }
  TGraphAsymmErrors* cf_graphWidth = new TGraphAsymmErrors(*cf_graph);
  TGraphErrors* coulomb = (TGraphErrors*) cfFile->FindObjectAny("Coulomb");
  TGraphErrors* halRad = (TGraphErrors*) cfFile->FindObjectAny("HalAndRad");
  TGraphErrors* halOnly = (TGraphErrors*) cfFile->FindObjectAny("HalOnly");
//  TGraphErrors* esc = (TGraphErrors*) cfFile->FindObjectAny("ESC");

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
  TFile* sysResFile = TFile::Open(sysResName, "read");
  TH1F* systResErr = (TH1F*) sysResFile->FindObjectAny("SystErrRel");
  if (!systResErr) {
    std::cout << "systResErr not found \n";
    sysResFile->ls();
    return 0;
  }

  TH1F* SystError = (TH1F*) cf_default->Clone("Sytematics");
  SystError->Reset();

  for (int iBin = 1; iBin <= SystError->GetNbinsX(); ++iBin) {
    double kStar = cf_default->GetBinCenter(iBin);
    double Ck = cf_default->GetBinContent(iBin);
    double x,y;
    cf_graph->GetPoint(iBin-1, x, y);
    cf_graph->SetPoint(iBin-1, x, Ck);

    double xErrLeft = (x - kStar + cf_default->GetBinWidth(iBin) / 2.)*0.95;
    double xErrRight = (kStar - x + cf_default->GetBinWidth(iBin) / 2.)*0.95;
    cf_graphWidth->SetPoint(iBin-1, x, Ck);
    cf_graphWidth->SetPointError(iBin-1, xErrLeft, xErrRight, 0., 0.);

    double errSystData = systDataErr->Eval(kStar);
    double errSystNorm = systNormErr->GetBinContent(iBin);
    double errSystLam = systLamErr->Eval(kStar);
    double errSystMom = systResErr->GetBinContent(iBin);
    double totErr = TMath::Sqrt(
        errSystData * errSystData + errSystNorm * errSystNorm
            + errSystLam * errSystLam + errSystMom * errSystMom);
    SystError->SetBinContent(iBin, totErr);
  }

  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("p-#Xi^{-} #bf{ALICE} data");
  LegNames.push_back("Coulomb + HAL QCD");
  LegNames.push_back("Coulomb");
  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  LegOptions.push_back("f");
  LegOptions.push_back("f");

  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = c1->DrawFrame(0, 0.7, 305, 3.7);
  const char * texPtY = "#it{C}(#it{k}*)";
  const char * texPtX = "#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetNdivisions(806);
  TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.1);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.1);
  pad->SetBottomMargin(0.115);
  pad->Draw();
  pad->cd();
  float fXmin = 0;
  float fXmax = 305;
  float fTextXMin = 0.35;
  float ymaxL = 0.81;
  DreamData *Data = new DreamData(Form("Data"));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
//  Data->SetCorrelationFunction(cf_default);
  Data->SetCorrelationGraph(cf_graph);
  Data->SetSystematics(SystError, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(false);
  Data->FemtoModelFitBands(halOnly, kPink + 1, 10, 0, -4000, true, false);
  Data->FemtoModelFitBands(halRad, kGray + 1, 0.5, false);
  Data->FemtoModelFitBands(coulomb, kGreen + 1, 1, 3, -4000, true, false);
//  Data->FemtoModelFitBands(esc,     11, 8,  0, -4000, true);
  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_default->GetMaximum() * 1.5);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.45, 0.6, 0.8, ymaxL + 0.03);
  Data->DrawCorrelationPlot(pad);

  pad->cd();
  cf_graphWidth->SetLineWidth(1);
  cf_graphWidth->SetLineColorAlpha(kBlack, 0.9);
  cf_graphWidth->Draw("same []");

  out->cd();
  c1->Write();
  c1->SaveAs(Form("CF_pXi.pdf"));
  SystError->Write();
  out->Close();
  app.Run();
}
