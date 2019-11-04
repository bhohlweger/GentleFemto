#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "TLatex.h"
#include "TStyle.h"

int main(int argc, char *argv[]) {
  const char* cf_fileName = argv[1];
  const char* sysData_fileName = argv[2];
  const char* varName = argv[3];

  TFile* cfFile = TFile::Open(cf_fileName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
      "InputCF_ResGami_woBL_ResGami_GenuineGami");
  if (!cf_default) {
    std::cout << "Default  data not found \n";
    cfFile->ls();
    return 0;
  }
  TFile* sysFile = TFile::Open(sysData_fileName, "read");
  TF1* systErr = (TF1*) sysFile->FindObjectAny("SystError");
  if (!systErr) {
    std::cout << "systErr not found \n";
    systErr->ls();
    return 0;
  }
  TColor myColor1;
  TFile* shittymodel = TFile::Open("~/cernbox/debug.root", "read");
  TGraph* coulomb = (TGraph*)shittymodel->Get("coulomb");
  coulomb->SetMarkerColor( myColor1.GetColor(178,223,138) );
  coulomb->SetLineColor( myColor1.GetColor(178,223,138) );
  coulomb->SetLineWidth(2);
  TGraph* halqcd = (TGraph*)shittymodel->Get("halqcd");
  halqcd->SetMarkerColor(myColor1.GetColor(255,127,0) );
  halqcd->SetLineColor( myColor1.GetColor(255,127,0) );
  halqcd->SetLineWidth(2);
  TGraph* esc = (TGraph*)shittymodel->Get("esc");
  esc->SetMarkerColor(    myColor1.GetColor(31,120,180)  );
  esc->SetLineColor(     myColor1.GetColor(31,120,180) );
  esc->SetLineWidth(2);
  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("p-#Xi^{-} #oplus #bar{p}-#bar{#Xi}^{+}");
  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");

  c1 = new TCanvas("c2", "c2", 0, 0, 500, 800);
  TFile* out = TFile::Open(Form("out_%s.root",varName), "recreate");
  TPad* pad;
  float LatexX = 0.;
  c1 = new TCanvas(Form("c"), Form("c"), 0, 0, 650, 650);
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.025);
  pad->SetTopMargin(0.025);
  pad->SetBottomMargin(0.12);
  pad->Draw();
  pad->cd();
  float fXmin = 0;
  float fXmax = 230;
  float fTextXMin = 0.35;
  DreamData *Data = new DreamData(Form("Data"));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
  Data->SetCorrelationFunction(cf_default);
  Data->SetSystematics(systErr, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(true);
  Data->SetRangePlotting(fXmin, fXmax, 0.7, cf_default->GetMaximum() * 1.17);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(legXmin, 0.625 - 0.09 * Data->GetNumberOfModels(),
                             legXmin + 0.4, 0.66);
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

  TLatex text;
  text.SetNDC();
  text.SetTextColor(1);
  text.SetTextSize(gStyle->GetTextSize() * 0.85);
  c1->cd();
  coulomb->Draw("lsame");
  halqcd->Draw("lsame");
  esc->Draw("lsame");
  out->cd();
  c1->Write();
  c1->SaveAs(Form("WithSys%s.pdf",varName));
}
