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

VariationmTAnalysis::VariationmTAnalysis(int nModels, int nData, int nVars)
    : fAnalysis(nModels),
      fnModel(nModels),
      fnData(nData),
      fnVars(nVars),
      fSystematic(),
      fHistname(),
      fFileName(),
      fmTAverage(),
      fmTRadiusSyst(),
      fmTRadiusStat() {
  for (int iMod = 0; iMod < fnModel; ++iMod) {
    fmTRadiusSyst.emplace_back(new TGraphErrors());
    fmTRadiusStat.emplace_back(new TGraphErrors());
  }
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

void VariationmTAnalysis::SetVariation(const char* VarDir, int iModel) {
  VariationAnalysis analysis = VariationAnalysis(fHistname, fnData, fnVars);
  TString filename = Form("%s/%s", VarDir, fFileName);
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
  const int iPoint = fmTRadiusSyst[iModel]->GetN();

  double mT, dummy;
  fmTAverage->GetPoint(iPoint, dummy, mT);
  fmTRadiusSyst[iModel]->SetPoint(iPoint, mT, radius);
  fmTRadiusStat[iModel]->SetPoint(iPoint, mT, radius);
  fmTRadiusSyst[iModel]->SetPointError(iPoint, fmTAverage->GetErrorY(iPoint),
                                       radiusErrSyst);
  fmTRadiusStat[iModel]->SetPointError(iPoint, fmTAverage->GetErrorY(iPoint),
                                       radiusErrStat);

  fAnalysis[iModel].push_back(analysis);
  return;
}

void VariationmTAnalysis::MakeCFPlotsPP() {
  DreamPlot::SetStyle();
  gStyle->SetLabelSize(16, "xyz");
  gStyle->SetTitleSize(16, "xyz");
  gStyle->SetTitleOffset(3.5, "x");
  gStyle->SetTitleOffset(3.5, "y");
  TGraphErrors* AxisGraph = new TGraphErrors();
  AxisGraph->SetPoint(0, 4, 1.);
  AxisGraph->SetPoint(1, 210, 1);
  AxisGraph->SetLineColor(kWhite);
//  AxisGraph->GetYaxis()->SetTitleOffset(1.5);
  AxisGraph->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  AxisGraph->GetXaxis()->SetRangeUser(4, 210);
  AxisGraph->GetYaxis()->SetRangeUser(0.725, 4.3);
  AxisGraph->GetXaxis()->SetNdivisions(505);
  auto c1 = new TCanvas("c2", "c2", 0, 0, 500, 800);
//  c1->Divide(4, 2);
  int counter = 1;
  TFile* out = TFile::Open("tmp.root", "recreate");
  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5 };
  for (auto it : fSystematic) {
    c1->cd();
    TPad* pad = GetFormattedPad(counter);
    pad->SetTopMargin(0.);
    float LatexX = 0.;
    //left sided pads
    if (counter % 2 == 0) {
      LatexX = 0.35;
    } else {  //right sided pads
      LatexX = 0.25;
    }
    pad->Draw();
    pad->cd();
    AxisGraph->Draw("Ap");
    DreamData *ProtonProton = new DreamData(Form("ProtonProton%i", counter));
    ProtonProton->SetMultiHisto(true);
    ProtonProton->SetUnitConversionData(1);
    ProtonProton->SetUnitConversionCATS(1);
    ProtonProton->SetCorrelationFunction(it.GetDefault());
    ProtonProton->SetSystematics(it.GetSystematicError(), 2);
    ProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p}", "fpe");
    ProtonProton->SetLegendName("#splitline{Coulomb +}{Argonne #nu_{18} (fit)}",
                                "l");
    ProtonProton->SetDrawAxis(false);
    ProtonProton->SetRangePlotting(4, 208, 0.725, 4.3);
    ProtonProton->SetNDivisions(505);
    ProtonProton->FemtoModelFitBands(fAnalysis[0][counter - 1].GetModel(), 2, 1,
                                     3, -3000, true);
    ProtonProton->SetLegendCoordinates(0., 0.2, 1.0, 0.8, false);
    ProtonProton->DrawCorrelationPlot(pad, 0, kBlack, 1.8);
    TLatex text;
    text.SetTextFont(43);
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSizePixels(14);
    text.DrawLatex(
        LatexX,
        0.8,
        TString::Format("m_{T} #in [%.2f, %.2f] (GeV/#it{c}^{2})",
                        mTppBins[counter - 1], mTppBins[counter]));
    if (counter == 1) {
      TPad* tmp2 = GetFormattedPad(0);
      c1->cd();
      tmp2->SetFillStyle(4000);
      tmp2->Draw();
      tmp2->cd();
      ProtonProton->DrawLegendExternal(tmp2);
    }

    counter++;
  }
  out->cd();
  c1->Write();
  c1->SaveAs("mTPlots.pdf");
  out->Close();
}
void VariationmTAnalysis::MakeRadPlotsPP() {
  TFile* out = TFile::Open("tmp.root", "update");
  auto c4 = new TCanvas("c8", "c8");
  c4->cd();
  fmTRadiusSyst[0]->SetLineColor(kBlack);
  fmTRadiusSyst[0]->SetTitle("; < m_{T} >  (MeV/#it{c}^{2}); r_{Core} (fm)");

  fmTRadiusSyst[0]->GetXaxis()->SetTitleSize(22);
  fmTRadiusSyst[0]->GetYaxis()->SetTitleSize(22);
  fmTRadiusSyst[0]->GetXaxis()->SetTitleOffset(1.5);
  fmTRadiusSyst[0]->GetYaxis()->SetTitleOffset(1.5);

  fmTRadiusSyst[0]->GetXaxis()->SetLabelSize(22);
  fmTRadiusSyst[0]->GetYaxis()->SetLabelSize(22);
  fmTRadiusSyst[0]->GetXaxis()->SetLabelOffset(.02);
  fmTRadiusSyst[0]->GetYaxis()->SetLabelOffset(.02);

  fmTRadiusSyst[0]->GetXaxis()->SetRangeUser(0.95, 2.7);
  fmTRadiusSyst[0]->GetYaxis()->SetRangeUser(0.65, 1.2);
  //  fmTRadiusSyst->GetXaxis()->SetRangeUser(0.95, 2.7);
  //  fmTRadiusSyst->GetYaxis()->SetRangeUser(0.95, 1.55);

  fmTRadiusSyst[0]->SetMarkerColorAlpha(kBlack, 0.);
  fmTRadiusSyst[0]->SetLineWidth(0);
  fmTRadiusSyst[0]->Draw("APZ");
  fmTRadiusSyst[0]->SetFillColorAlpha(kBlack, 0.4);
  fmTRadiusSyst[0]->Draw("2Z same");
  TGraphErrors fakeGraph;
  fakeGraph.SetMarkerColor(kBlack);
  fakeGraph.SetLineWidth(3);
  fakeGraph.SetDrawOption("z");
  fakeGraph.SetFillColorAlpha(kBlack, 0.4);
  TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->SetFillStyle(4000);
  leg->AddEntry(&fakeGraph, "p#minus p (AV18)", "lef");

  fmTRadiusStat[0]->SetMarkerColor(kBlack);
  fmTRadiusStat[0]->SetLineWidth(3);
  fmTRadiusStat[0]->Draw("pez same");
  leg->Draw("same");
  c4->SaveAs("mTvsRad.pdf");
  c4->Write();
  fmTRadiusSyst[0]->SetName("mTRadiusSyst");
  fmTRadiusSyst[0]->Write();
  fmTRadiusStat[0]->SetName("mTRadiusStat");
  fmTRadiusStat[0]->Write();

  out->Write();
  out->Close();
}

void VariationmTAnalysis::MakeRadPlotsPL(const char* ppFilePath) {
  TFile* ppFile = TFile::Open(ppFilePath, "read");
  TGraphErrors* mTppSys = (TGraphErrors*) ppFile->Get("mTRadiusSyst");
  TGraphErrors* mTppStat = (TGraphErrors*) ppFile->Get("mTRadiusStat");
  TFile* out = TFile::Open("tmp.root", "update");
  out->cd();
  auto c4 = new TCanvas("c8", "c8");
  c4->cd();
  TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->SetFillStyle(4000);
  TGraphErrors fakeGraphUsmani;
  fakeGraphUsmani.SetLineWidth(3);
  fakeGraphUsmani.SetDrawOption("z");
  TGraphErrors fakeGraphNLO;
  fakeGraphNLO.SetLineWidth(3);
  fakeGraphNLO.SetDrawOption("z");
  TGraphErrors fakeGraphLO;
  fakeGraphLO.SetLineWidth(3);
  fakeGraphLO.SetDrawOption("z");

  mTppSys->SetTitle("; < m_{T} >  (MeV/#it{c}^{2}); r_{Core} (fm)");
  mTppSys->GetXaxis()->SetTitleSize(22);
  mTppSys->GetYaxis()->SetTitleSize(22);
  mTppSys->GetXaxis()->SetTitleOffset(1.5);
  mTppSys->GetYaxis()->SetTitleOffset(1.5);
  mTppSys->GetXaxis()->SetLabelSize(22);
  mTppSys->GetYaxis()->SetLabelSize(22);
  mTppSys->GetXaxis()->SetLabelOffset(.02);
  mTppSys->GetYaxis()->SetLabelOffset(.02);
  mTppSys->GetXaxis()->SetRangeUser(0.95, 2.7);
  mTppSys->GetYaxis()->SetRangeUser(0.3, 1.35);

  mTppSys->SetMarkerColorAlpha(kBlack, 0.);
  mTppSys->SetLineWidth(0);
  mTppSys->Draw("APZ");
  mTppSys->SetFillColorAlpha(kBlack, 0.4);
  mTppSys->Draw("2Z same");
  TGraphErrors fakeGraph;
  fakeGraph.SetMarkerColor(kBlack);
  fakeGraph.SetLineWidth(3);
  fakeGraph.SetDrawOption("z");
  fakeGraph.SetFillColorAlpha(kBlack, 0.4);

  leg->AddEntry(&fakeGraph, "p#minus p (AV18)", "lef");

  mTppStat->SetMarkerColor(kBlack);
  mTppStat->SetLineWidth(3);
  mTppStat->Draw("pez same");

  for (int iMod = 0; iMod < fnModel; ++iMod) {

    //  fmTRadiusSyst->GetXaxis()->SetRangeUser(0.95, 2.7);
    //  fmTRadiusSyst->GetYaxis()->SetRangeUser(0.95, 1.55);

    fmTRadiusSyst[iMod]->SetLineWidth(0);
    if (iMod == 0) {
      fmTRadiusSyst[iMod]->SetLineColor(kCyan + 2);

      fmTRadiusSyst[iMod]->SetMarkerColorAlpha(kCyan + 2, 0.);
      fmTRadiusSyst[iMod]->Draw("PZSame");

      fmTRadiusSyst[iMod]->SetFillColorAlpha(kCyan + 2, 0.4);
      fmTRadiusSyst[iMod]->Draw("2Z same");

      fakeGraphUsmani.SetMarkerColor(kCyan + 2);
      fakeGraphUsmani.SetFillColorAlpha(kCyan + 2, 0.4);
      fmTRadiusStat[iMod]->SetMarkerColor(kCyan + 2);
      leg->AddEntry(&fakeGraphUsmani, "p#minus #Lambda (Usmani)", "lef");
      fmTRadiusStat[iMod]->SetLineWidth(3);
    } else if (iMod == 1) {
      fmTRadiusSyst[iMod]->SetLineColor(kRed + 1);

      fmTRadiusSyst[iMod]->SetMarkerColorAlpha(kRed + 1, 0.);
      fmTRadiusSyst[iMod]->Draw("PZSame");

      fmTRadiusSyst[iMod]->SetFillColorAlpha(kRed + 1, 0.4);
      fmTRadiusSyst[iMod]->Draw("2Z same");

      fakeGraphNLO.SetMarkerColor(kRed + 1);
      fakeGraphNLO.SetFillColorAlpha(kRed + 1, 0.4);
      fmTRadiusStat[iMod]->SetMarkerColor(kRed + 1);
      leg->AddEntry(&fakeGraphNLO, "p#minus #Lambda (#chi EFT NLO)", "lef");
      fmTRadiusStat[iMod]->SetLineWidth(3);
    } else {
      fmTRadiusSyst[iMod]->SetLineColor(kGreen + 3);

      fmTRadiusSyst[iMod]->SetMarkerColorAlpha(kGreen + 3, 0.);
      fmTRadiusSyst[iMod]->Draw("PZSame");

      fmTRadiusSyst[iMod]->SetFillColorAlpha(kGreen + 3, 0.4);
      fmTRadiusSyst[iMod]->Draw("2Z same");

      fakeGraphLO.SetMarkerColor(kGreen + 3);
      fakeGraphLO.SetFillColorAlpha(kGreen + 3, 0.4);
      fmTRadiusStat[iMod]->SetMarkerColor(kGreen + 3);
      leg->AddEntry(&fakeGraphLO, "p#minus #Lambda (#chi EFT LO)", "lef");
      fmTRadiusStat[iMod]->SetLineWidth(3);
    }
    fmTRadiusStat[iMod]->Draw("pez same");
    fmTRadiusSyst[iMod]->SetName(TString::Format("mTRadiusSystMod_%u", iMod));
    fmTRadiusSyst[iMod]->Write();
    fmTRadiusStat[iMod]->SetName(TString::Format("mTRadiusStatMod_%u", iMod));
    fmTRadiusStat[iMod]->Write();
  }
  leg->Draw("same");
  c4->SaveAs("mTvsRad.pdf");
  c4->Write();
  out->Write();
  out->Close();
  ppFile->Close();
}

void VariationmTAnalysis::MakeCFPlotsPL() {
  DreamPlot::SetStyle();
//  gStyle->SetLabelSize(16, "xyz");
//  gStyle->SetTitleSize(16, "xyz");
//  gStyle->SetTitleOffset(3.5, "x");
//  gStyle->SetTitleOffset(3.5, "y");
  TFile* out = TFile::Open("tmp.root", "recreate");
  std::vector<float> mTppBins = { 1.08, 1.26, 1.32, 1.44, 1.65, 1.9, 4.5 };

  int counter = 1;
  for (auto it : fSystematic) {
    int iModCounter = 0;
    for (auto iModel : fAnalysis) {
      auto c1 = new TCanvas(TString::Format("c%u_%u", counter, iModCounter),
                            TString::Format("c%u_%u", counter, iModCounter), 0,
                            0, 650, 650);
      c1->cd();
      TPad *p1 = new TPad("p1", "p1", 0., 0., 1., 1.);
      p1->SetRightMargin(0.025);
      p1->SetTopMargin(0.025);
      p1->SetBottomMargin(0.12);
      p1->Draw();
      p1->cd();
      DreamData *ProtonLambda = new DreamData(
          Form("ProtonLambda%u_%u", counter, iModCounter));
      ProtonLambda->SetMultiHisto(false);
      ProtonLambda->SetUnitConversionData(1);
      ProtonLambda->SetUnitConversionCATS(1);
      ProtonLambda->SetCorrelationFunction(it.GetDefault());
      ProtonLambda->SetSystematics(it.GetSystematicError(), 2);
      ProtonLambda->SetLegendName("p-#Lambda #oplus #bar{p}-#bar{#Lambda}",
                                  "fpe");
      if (iModCounter == 0) {
        ProtonLambda->SetLegendName("Usmani", "fl");
        ProtonLambda->FemtoModelFitBands(
            fAnalysis[iModCounter][counter - 1].GetModel(), 11, 7, 0, 3244,
            true);
      } else if (iModCounter == 1) {
        ProtonLambda->SetLegendName("#chi_{EFT} NLO", "fl");
        ProtonLambda->FemtoModelFitBands(
            fAnalysis[iModCounter][counter - 1].GetModel(), 1, 7, 0, 3244,
            true);
      } else if (iModCounter == 2) {
        ProtonLambda->SetLegendName("#chi_{EFT} LO", "fl");
        ProtonLambda->FemtoModelFitBands(
            fAnalysis[iModCounter][counter - 1].GetModel(), 3, 7, 0, 3244,
            true);
      }
      ProtonLambda->SetDrawAxis(true);
      ProtonLambda->SetRangePlotting(4, 208, 0.87, 2.4);
      ProtonLambda->SetNDivisions(505);
      ProtonLambda->SetLegendCoordinates(
          0.4, 0.6 - 0.09 * ProtonLambda->GetNumberOfModels(), 0.7, 0.665);
      ProtonLambda->DrawCorrelationPlot(p1);

      p1->cd();
      TLatex BeamText;
      TLatex text;
      BeamText.SetTextSize(gStyle->GetTextSize() * .85);
      BeamText.SetNDC(kTRUE);
      BeamText.DrawLatex(0.42, 0.91, Form("#bf{ALICE}"));
      BeamText.DrawLatex(0.42, 0.85,
                         Form("%s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
      BeamText.DrawLatex(0.42, 0.79, "High Mult. (0-0.072% INEL)");
      text.SetNDC();
      text.SetTextColor(1);
      text.SetTextSize(gStyle->GetTextSize() * 0.85);
      text.DrawLatex(0.42, 0.73, "Gaussian + Resonance source");
      text.DrawLatex(
          0.42,
          0.67,
          TString::Format("m_{T} #in [%.2f, %.2f] (GeV/#it{c}^{2})",
                          mTppBins[counter - 1], mTppBins[counter]));
      out->cd();
      c1->Write();
      c1->SaveAs(Form("mTPlots%u_%u.pdf", counter, iModCounter));
      iModCounter++;
    }
    counter++;
  }
  out->Close();
}

TPad* VariationmTAnalysis::GetFormattedPad(int counter) {
  std::vector<float> xMinPad = { 0.1, 0.4, 0., 0.5, 0., 0.5, 0., 0.5 };
  std::vector<float> xMaxPad = { 0.4, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0 };
  std::vector<float> yMinPad = { 0.75, 0.75, 0.52, 0.52, 0.29, 0.29, 0., 0. };
  std::vector<float> yMaxPad = { .98, 0.98, 0.75, 0.75, 0.52, 0.52, 0.29, 0.29 };
  TPad* pad = new TPad(Form("p%u", counter), Form("p%u", counter),
                       xMinPad[counter], yMinPad[counter], xMaxPad[counter],
                       yMaxPad[counter]);
  pad->SetTopMargin(0.);
  //left sided pads
  if (counter % 2 == 0) {
    pad->SetRightMargin(0.);
    pad->SetLeftMargin(0.2);
    if (counter < 5) {
      pad->SetBottomMargin(0.);
    } else {
      pad->SetBottomMargin(0.06 / 0.29);
    }
  } else {  //right sided pads
    if (counter != 1) {
      pad->SetLeftMargin(0.);
      pad->SetRightMargin(0.07);
    } else {
      pad->SetLeftMargin(0.1 / 0.6);
      pad->SetRightMargin(0.035 / 0.6);
    }
    if (counter < 6) {
      pad->SetBottomMargin(0.);
    } else {
      pad->SetBottomMargin(0.06 / 0.29);
    }
  }
  return pad;
}
