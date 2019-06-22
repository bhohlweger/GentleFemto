/*
 * VariationmTAnalysis.cxx
 *
 *  Created on: Jun 18, 2019
 *      Author: schmollweger
 */
#include "VariationmTAnalysis.h"
#include "TSystemDirectory.h"
#include "DreamPlot.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include <iostream>

VariationmTAnalysis::VariationmTAnalysis()
    : fAnalysis(),
      fSystematic(),
      fHistname(),
      fmTAverage(),
      fmTRadiusSyst(new TGraphErrors()),
      fmTRadiusStat(new TGraphErrors()) {
  // TODO Auto-generated constructor stub

}

VariationmTAnalysis::~VariationmTAnalysis() {
  // TODO Auto-generated destructor stub
}

void VariationmTAnalysis::SetSystematic(const char* DataDir) {
  static int outputCounter = 0;
  TSystemDirectory *workdir = new TSystemDirectory("workdir", DataDir);
  TList *RootList = workdir->GetListOfFiles();
  RootList->Sort();
  TIter next(RootList);
  TObject* obj = nullptr;
  DreamSystematics Systematics(DreamSystematics::pp);
  Systematics.SetUpperFitRange(150);
  Systematics.SetBarlowUpperRange(150);
  while (obj = next()) {
    TString FileName = obj->GetName();
    if (FileName.Contains(".root")) {
      TH1F* histo = nullptr;
      TFile* File = TFile::Open(
          TString::Format("%s/%s", DataDir, FileName.Data()).Data(), "read");
      if (!File) {
        Warning(
            "VariationmTAnalysis::SetSystematic",
            TString::Format("File %s does not exist, exiting \n",
                            FileName.Data()));
        return;
      }
      TList* FileKeys = File->GetListOfKeys();
      TIter FileIter(FileKeys);
      TObject* FileObj;
      while (FileObj = FileIter()) {
        TString FileObjName = FileObj->GetName();
        if (FileObjName.Contains(fHistname) && FileObjName.Contains("MeV")) {
          histo = (TH1F*) (File->FindObjectAny(FileObjName.Data()))->Clone(
              TString::Format("%sClone", FileObjName.Data()));
          histo->SetDirectory(0);
          break;
        }
      }
      if (!histo) {
        Warning(
            "VariationmTAnalysis::SetSystematic",
            TString::Format("No Histogram found for %s in file %s. Exiting \n",
                            fHistname, FileName.Data()).Data());
        return;
      }
      if (FileName.Contains("Var0")) {
        //that's the default
        Systematics.SetDefaultHist(histo);
      } else {
        Systematics.SetVarHist(histo);
      }
      File->Close();
    }
  }
  Systematics.EvalSystematics();
  Systematics.WriteOutput(Form("%u", outputCounter));
  outputCounter++;
  fSystematic.push_back(Systematics);
  return;
}

void VariationmTAnalysis::SetVariation(const char* VarDir) {
  VariationAnalysis analysis = VariationAnalysis(fHistname, 26, 81);
  TString filename = Form("%s/OutFileVarpp.root", VarDir);
  analysis.ReadFitFile(filename.Data());
  analysis.EvalRadius();
  float radius = analysis.GetRadMean();
  float radiusErrStat = analysis.GetRadStatErr();
  float radiusErrSyst = (analysis.GetRadSystDown() + analysis.GetRadSystUp())
      / 2.;
  if (!fmTAverage) {
    Warning("SetVariation", "No Average mT histo set, exiting \n");
    return;
  }
  const int iPoint = fmTRadiusSyst->GetN();
  double mT, dummy;
  fmTAverage->GetPoint(iPoint, dummy, mT);
  fmTRadiusSyst->SetPoint(iPoint, mT, radius);
  fmTRadiusStat->SetPoint(iPoint, mT, radius);
  fmTRadiusSyst->SetPointError(iPoint, fmTAverage->GetErrorY(iPoint),
                               radiusErrSyst);
  fmTRadiusStat->SetPointError(iPoint, fmTAverage->GetErrorY(iPoint),
                               radiusErrStat);

  fAnalysis.push_back(analysis);
  return;
}

void VariationmTAnalysis::MakePlots() {
  DreamPlot::SetStyle();
  gStyle->SetLabelSize(14, "xyz");
  gStyle->SetTitleSize(14, "xyz");
  gStyle->SetTitleOffset(2., "x");
  gStyle->SetTitleOffset(3., "y");
  auto c1 = new TCanvas("c2", "c2", 0, 0, 500, 800);
//  c1->Divide(4, 2);
  int counter = 1;
  TFile* out = TFile::Open("tmp.root", "recreate");
  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5 };
  std::vector<float> xMinPad = { 0., 0.5, 0., 0.5, 0., 0.5, 0., 0.5 };
  std::vector<float> xMaxPad = { 0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0 };
  std::vector<float> yMinPad = { 0.75, 0.75, 0.5, 0.5, 0.25, 0.25, 0., 0. };
  std::vector<float> yMaxPad = { 1.0, 1.0, 0.75, 0.75, 0.5, 0.5, 0.25, 0.25 };
  for (auto it : fSystematic) {
    c1->cd();
    TPad* pad = new TPad(Form("p%u", counter), Form("p%u", counter),
                         xMinPad[counter - 1], yMinPad[counter - 1],
                         xMaxPad[counter - 1], yMaxPad[counter - 1]);
    pad->SetTopMargin(0.);
    if (counter % 2 == 0 ) {
      pad->SetRightMargin(0.07);
      pad->SetLeftMargin(0.);
      if (counter < 6) {
        pad->SetBottomMargin(0.);
      }
    } else {
      pad->SetRightMargin(0.);
      pad->SetLeftMargin(0.2);
      if (counter < 6) {
        pad->SetBottomMargin(0.);
      }
    }
    pad->Draw();
    pad->cd();
    DreamData *ProtonProton = new DreamData(Form("ProtonProton%i", counter));
    ProtonProton->SetMultiHisto(true);
    ProtonProton->SetUnitConversionData(1);
    ProtonProton->SetUnitConversionCATS(1);
    ProtonProton->SetCorrelationFunction(it.GetDefault());
    ProtonProton->SetSystematics(it.GetSystematicError(), 2);
    ProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p}", "fpe");
    ProtonProton->SetLegendName("#splitline{Coulomb +}{Argonne #nu_{18} (fit)}",
                                "l");
    ProtonProton->SetRangePlotting(0, 200, 0.725, 4.3);
    ProtonProton->SetNDivisions(505);
    ProtonProton->FemtoModelFitBands(fAnalysis[counter - 1].GetModel(), 2, 1, 3,
                                     -3000, true);
    ProtonProton->SetLegendCoordinates(0., 0.2, 1.0, 0.8, false);
    ProtonProton->DrawCorrelationPlot(pad, 0, kBlack, 1.2);
    TLatex text;
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSize(gStyle->GetTextSize() * 1.2);
    text.DrawLatex(
        0.35,
        0.85,
        TString::Format("#splitline{m_{T} #in}{[%.2f, %.2f] (GeV/#it{c}^{2})}",
                        mTppBins[counter - 1], mTppBins[counter]));
//    if (counter == 1) {
//      TPad* tmp2 = (TPad*) c1->cd(8);
//      ProtonProton->DrawLegendExternal(tmp2);
//    }
    counter++;
  }
  out->cd();
  c1->Write();
  c1->SaveAs("mTPlots.pdf");

  auto c4 = new TCanvas("c8", "c8");
  c4->cd();
  fmTRadiusSyst->SetLineColor(kBlack);
  fmTRadiusSyst->SetTitle("; < m_{T} >  (MeV/#it{c}^{2}); r_{Core} (fm)");

  fmTRadiusSyst->GetXaxis()->SetTitleSize(22);
  fmTRadiusSyst->GetYaxis()->SetTitleSize(22);
  fmTRadiusSyst->GetXaxis()->SetTitleOffset(1.5);
  fmTRadiusSyst->GetYaxis()->SetTitleOffset(1.5);

  fmTRadiusSyst->GetXaxis()->SetLabelSize(22);
  fmTRadiusSyst->GetYaxis()->SetLabelSize(22);
  fmTRadiusSyst->GetXaxis()->SetLabelOffset(.02);
  fmTRadiusSyst->GetYaxis()->SetLabelOffset(.02);

  fmTRadiusSyst->GetXaxis()->SetRangeUser(0.95, 2.7);
  fmTRadiusSyst->GetYaxis()->SetRangeUser(0.65, 1.2);
  //  fmTRadiusSyst->GetXaxis()->SetRangeUser(0.95, 2.7);
  //  fmTRadiusSyst->GetYaxis()->SetRangeUser(0.95, 1.55);

  fmTRadiusSyst->SetMarkerColorAlpha(kBlack, 0.);
  fmTRadiusSyst->SetLineWidth(0);
  fmTRadiusSyst->Draw("APZ");
  fmTRadiusSyst->SetFillColorAlpha(kBlack, 0.4);
  fmTRadiusSyst->Draw("2Z same");
  TGraphErrors fakeGraph;
  fakeGraph.SetMarkerColor(kBlack);
  fakeGraph.SetLineWidth(3);
  fakeGraph.SetDrawOption("z");
  fakeGraph.SetFillColorAlpha(kBlack, 0.4);
  TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->SetFillStyle(4000);
  leg->AddEntry(&fakeGraph, "p#minus p (AV18)", "lef");

  fmTRadiusStat->SetMarkerColor(kBlack);
  fmTRadiusStat->SetLineWidth(3);
  fmTRadiusStat->Draw("pez same");
  leg->Draw("same");
  c4->SaveAs("mTvsRad.pdf");
  c4->Write();
  fmTRadiusSyst->SetName("mTRadiusSyst");
  fmTRadiusSyst->Write();
  fmTRadiusStat->SetName("mTRadiusStat");
  fmTRadiusStat->Write();

  out->Write();
  out->Close();
}
