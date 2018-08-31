/*
 * DreamData.cxx
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#include "DreamData.h"
#include "TLegend.h"
#include "TStyle.h"
DreamData::DreamData(const char* particlePair)
    : fName(particlePair),
      fCorrelationFunction(nullptr),
      fCorrelationFunctionSimulation(nullptr),
      fSystematics(nullptr),
      fSysError(nullptr),
      fBaseLine(new TF1(Form("%sBaseLine", particlePair), "pol1", 0, 1)),
      fXMin(0),
      fXMax(0.5),
      fYMin(0),
      fYMax(0.5),
      fLegendName(),
      fFemtoModdeled(),
      fFakeGraph() {
  fFillColors = {kGray + 1, kRed - 10, kBlue - 9, kGreen - 8, kMagenta - 9,
    kOrange - 9, kCyan - 3, kYellow - 7};
  fColors = {kBlack, kRed + 1, kBlue + 2, kGreen + 3, kMagenta + 1, kOrange - 1,
    kCyan + 2, kYellow + 2};
  fMarkers = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond,
    kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};
  fBaseLine->SetLineStyle(2);
  fBaseLine->SetLineColor(fFillColors[6]);

}

DreamData::~DreamData() {
  // TODO Auto-generated destructor stub
}

void DreamData::SetSystematics(TF1* parameters, int UnitConv,
                               float errorwidth) {
  if (parameters) {
//    TF1* UnbinnedSys = new TF1(Form("%sUnbinnedSyst", fName), "pol2", 0, 1);
//    UnbinnedSys->SetParameter(0, parameters->GetBinContent(1));
//    UnbinnedSys->SetParameter(1, parameters->GetBinContent(2));
//    UnbinnedSys->SetParameter(2, parameters->GetBinContent(3));

    if (fCorrelationFunction) {
      int nBinsX = fCorrelationFunction->GetNbinsX();
      float minX = fCorrelationFunction->GetXaxis()->GetXmin();
      float maxX = fCorrelationFunction->GetXaxis()->GetXmax();
      fSystematics = new TH1F(Form("Systematics%s", fName),
                              Form("Systematics%s", fName), nBinsX, minX, maxX);
      for (int iBin = 1; iBin < nBinsX; iBin++) {
        const float x = fCorrelationFunction->GetBinCenter(iBin);
        const float y = fCorrelationFunction->GetBinContent(iBin);
        fSystematics->SetBinContent(iBin,
                                    y * parameters->Eval(x / (float) UnitConv));
      }
      fSystematics->SetLineWidth(2.0);
      fSysError = new TGraphErrors();

      for (int i = 0; i < nBinsX; i++) {
        if (fCorrelationFunction->GetBinCenter(i + 1) > 0.2)
          continue;
        fSysError->SetPoint(i, fCorrelationFunction->GetBinCenter(i + 1),
                            fCorrelationFunction->GetBinContent(i + 1));
        fSysError->SetPointError(i, errorwidth,
                                 fSystematics->GetBinContent(i + 1));
      }
      TGraph *grFakeSys = new TGraph();
      grFakeSys->SetFillColor(fFillColors[0]);
      grFakeSys->SetLineColor(fFillColors[0]);
      SetStyleGraph(grFakeSys, 0, 0);
      fFakeGraph.push_back(grFakeSys);
    } else {
      std::cout << "For " << fName
                << " set the CF before adding the systematics \n";
    }
  } else {
    std::cout << "Paramters input missing for " << fName << std::endl;
  }
  return;
}

void DreamData::FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower,
                                   TGraph *grUpper, int UnitConv, int color) {
  if (fSystematics) {
    TGraphErrors *grFemtoModel = new TGraphErrors();
    grFemtoModel->SetName(grMedian1->GetName());
    double x, yM1, yLo, yUp;
    int count = 0;
    for (int i = 0; i < grMedian1->GetN(); ++i) {
      grMedian1->GetPoint(i, x, yM1);
      grLower->GetPoint(i, x, yLo);
      grUpper->GetPoint(i, x, yUp);
      std::vector<float> yAll;
      yAll.push_back(yM1);
      yAll.push_back(yLo);
      yAll.push_back(yUp);
      std::sort(yAll.begin(), yAll.end());
      grFemtoModel->SetPoint(count, x / (float) UnitConv,
                             (yAll[2] + yAll[0]) / 2.f);
      grFemtoModel->SetPointError(count++, 0,
                                  (yAll[2] + yAll[0]) / 2.f - yAll[0]);
    }
    grFemtoModel->SetFillColor(fColors[color]);
    grFemtoModel->SetLineColor(fColors[color]);
    grFemtoModel->SetLineWidth(3);
    fFemtoModdeled.push_back(grFemtoModel);
    TGraph *grFakeModel = new TGraph();
    grFakeModel->SetLineColor(fColors[color]);
    grFakeModel->SetLineWidth(4);
    fFakeGraph.push_back(grFakeModel);
  } else {
    std::cout << "Set Systematics first for " << fName << "\n";
  }
  return;
}

void DreamData::SetStyleHisto(TH1 *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

void DreamData::DrawCorrelationPlot(TCanvas* c) {
  c->cd();
  SetStyleHisto(fCorrelationFunction, 0, 0);
  fCorrelationFunction->GetXaxis()->SetRangeUser(0, 0.4);
  fCorrelationFunction->GetXaxis()->SetRangeUser(0, 4);
  fSysError->SetLineColor(kWhite);
  fSysError->Draw("Ap");
  fBaseLine->Draw("same");
  fSysError->SetTitle("; #it{k}* (GeV/#it{c}); #it{C}(#it{k}*)");
  fSysError->GetXaxis()->SetRangeUser(fXMin, fXMax);
  fSysError->GetYaxis()->SetRangeUser(fYMin, fYMax);
  TLegend *leg = new TLegend(0.48,0.535,0.85,0.75);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize() * 0.75);
  int legendCounter = 1;
//  leg->AddEntry(fCorrelationFunction, fLegendName[0], "pe");
  leg->AddEntry(fFakeGraph[0], fLegendName[0], "fpe");
  leg->AddEntry(fBaseLine, "BaseLine", "l");
  leg->Draw("same");
  for (auto &it : fFemtoModdeled) {
    it->Draw("L3 same");
    leg->AddEntry(fFakeGraph[legendCounter], fLegendName[legendCounter], "l");
    legendCounter++;
  }
  fSysError->SetFillColorAlpha(kBlack, 0.4);
  fSysError->Draw("2 same");
  fCorrelationFunction->Draw("pe same");
  leg->Draw("same");
}

void DreamData::SetStyleGraph(TGraph *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(1.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}
