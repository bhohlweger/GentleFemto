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
      fmTBins(),
      fHistname(),
      fFileName(),
      fDataName(),
      fDataOption(),
      fModelName(),
      fModelOption(),
      fSourceName(),
      fColor(),
      fTextXMin(0.32),
      fXmin(4),
      fXmax(250),
      fSystematic(),
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
void VariationmTAnalysis::MakeCFPlotsSingleBand() {
  DreamPlot::SetStyle();
//  gStyle->SetLabelSize(16, "xyz");
//  gStyle->SetTitleSize(16, "xyz");
//  gStyle->SetTitleOffset(3.5, "x");
//  gStyle->SetTitleOffset(3.5, "y");
  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back(fDataName);
  LegNames.insert(LegNames.end(), fModelName.begin(), fModelName.end());
  std::vector<const char*> LegOptions;
  LegOptions.push_back(fDataOption);
  LegOptions.insert(LegOptions.end(), fModelOption.begin(), fModelOption.end());

  c1 = new TCanvas("c2", "c2", 0, 0, 500, 800);
  int counter = 1;
  TFile* out = TFile::Open("tmp.root", "recreate");
  for (auto it : fSystematic) {
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
    Data->SetCorrelationFunction(it.GetDefault());
    Data->SetSystematics(it.GetSystematicError(), 2);
    Data->SetLegendName(LegNames, LegOptions);
    Data->SetDrawAxis(true);
    Data->SetRangePlotting(fXmin, fXmax, 0.9, it.GetDefault()->GetMaximum()*1.2); //ranges
    Data->SetNDivisions(505);
    for (int iMod = 0; iMod < fnModel; ++iMod) {
      Data->FemtoModelFitBands(fAnalysis[iMod][counter - 1].GetModel(), fColor[iMod], 1, 3,
                               -3000, true); //Model colors
    }
    float legXmin = fTextXMin-0.02;
    Data->SetLegendCoordinates(legXmin, 0.67 - 0.09 * Data->GetNumberOfModels(),
                               legXmin + 0.4, 0.725);
    Data->DrawCorrelationPlot(pad);
    pad->cd();
    TLatex BeamText;
    BeamText.SetTextSize(gStyle->GetTextSize() * .85);
    BeamText.SetNDC(kTRUE);
    BeamText.DrawLatex(fTextXMin, 0.91,
                       Form("ALICE %s #sqrt{#it{s}} = %i TeV", "pp", (int) 13));
    BeamText.DrawLatex(fTextXMin, 0.85, "High Mult. (0-0.072% INEL)");

    TLatex text;
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSize(gStyle->GetTextSize() * 0.85);
    text.DrawLatex(fTextXMin, 0.79, fSourceName);
    text.DrawLatex(
        fTextXMin,
        0.73,
        TString::Format("m_{T} #in [%.2f, %.2f] (GeV/#it{c}^{2})",
                        fmTBins[counter - 1], fmTBins[counter]));
    out->cd();
    c1->Write();
    c1->SaveAs(Form("mTBin_%u.pdf", counter));
    counter++;
  }
  out->Close();
}

void VariationmTAnalysis::MakeOnePanelPlots() {
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

void VariationmTAnalysis::StoreRadvsmT(const char* fileName) {
  TFile* out = TFile::Open(fileName, "recreate");
  fmTRadiusSyst[0]->SetName("mTRadiusSyst");
  fmTRadiusSyst[0]->Write();
  fmTRadiusStat[0]->SetName("mTRadiusStat");
  fmTRadiusStat[0]->Write();
  out->Write();
  out->Close();
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
