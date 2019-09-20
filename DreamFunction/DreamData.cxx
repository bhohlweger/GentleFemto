/*
 * DreamData.cxx
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#include "DreamData.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TRandom.h"
DreamData::DreamData(const char* particlePair)
    : fName(particlePair),
      fCorrelationFunction(nullptr),
      fCorrelationFunctionSimulation(nullptr),
      fCorrelationGraph(nullptr),
      fSystematics(nullptr),
      fSysError(nullptr),
      fCorrelatedError(nullptr),
      fDummyHist(nullptr),
      fBaseLine(new TF1(Form("%sBaseLine", particlePair), "pol1", 0, 1000)),
      fDrawAxis(true),
      fForceAxis(false),
      fXMin(0),
      fXMax(0.5),
      fYMin(0),
      fYMax(0.5),
      fInlet(false),
      fXMinZoom(0),
      fXMaxZoom(0.5),
      fYMinZoom(0),
      fYMaxZoom(0.5),
      fXMinInlet(0),
      fXMaxInlet(0.5),
      fYMinInlet(0),
      fYMaxInlet(0.5),
      fLegend(nullptr),
      fXMinLegend(0),
      fXMaxLegend(0.5),
      fYMinLegend(0),
      fYMaxLegend(0.5),
      fDrawLegend(true),
      fUnitConversionData(1),
      fUnitConversionCATS(1),
      fMultiHisto(false),
      fLegendName(),
      fLegendOption(),
      fFemtoModdeled(),
      fFakeGraph() {
  fFillColors = {
    kGray + 1,
    kRed - 10,
    kBlue - 9,
    kGreen - 8,
    kMagenta - 9,
    kOrange - 9,
    kCyan - 3,
    kYellow - 7,
    kBlue + 3
  };
  TColor myColor1;
  fColors = {
    kBlack,         //0
    kRed + 1,//1
    kBlue + 2,//2
    kGreen + 3,//3
    kMagenta - 8,//4
    kOrange - 7,//5
    kCyan + 2,//6
    kYellow + 2,//7
    kWhite,//8
    kGreen - 5,//9
    myColor1.GetColor(255,127,0),//10
    myColor1.GetColor(31,120,180),//11
    myColor1.GetColor(178,223,138),//12
    kBlue + 3//13
  };
  fMarkers = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond,
    kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};
  fBaseLine->SetLineStyle(7);
  fBaseLine->SetLineWidth(1);
  fBaseLine->SetLineColor(fColors[0]);

}

DreamData::~DreamData() {
  // TODO Auto-generated destructor stub
}

void DreamData::SetSystematics(TF1* parameters, float errorwidth) {
  if (parameters) {
    if (fCorrelationFunction) {
      int nBinsX = fCorrelationFunction->GetNbinsX();
      float minX = fCorrelationFunction->GetXaxis()->GetXmin();
      float maxX = fCorrelationFunction->GetXaxis()->GetXmax();
      fSystematics = new TH1F(Form("Systematics%s", fName),
                              Form("Systematics%s", fName), nBinsX, minX, maxX);
      for (int iBin = 1; iBin < nBinsX; iBin++) {
        const float x = fCorrelationFunction->GetBinCenter(iBin);
        const float y = fCorrelationFunction->GetBinContent(iBin);
//        std::cout << "x = " << x << '\t' << " y = " << y << '\n'
//                  << " rel error "
//                  << parameters->Eval(x / (float) fUnitConversionData) << '\t'
//                  << " Error = "
//                  << y * parameters->Eval(x / (float) fUnitConversionData)
//                  << std::endl;
        fSystematics->SetBinContent(
            iBin, y * parameters->Eval(x / (float) fUnitConversionData));
      }
      fSystematics->SetLineWidth(2.0);
      fSysError = new TGraphAsymmErrors();

      for (int i = 0; i < nBinsX; i++) {
        fSysError->SetPoint(i, fCorrelationFunction->GetBinCenter(i + 1),
                            fCorrelationFunction->GetBinContent(i + 1));
        fSysError->SetPointError(i, errorwidth, errorwidth,
                                 fSystematics->GetBinContent(i + 1),
                                 fSystematics->GetBinContent(i + 1));
      }
      TGraph *grFakeSys = new TGraph();
      SetStyleGraph(grFakeSys, 2, 0);
      grFakeSys->SetFillColor(fFillColors[0]);
      grFakeSys->SetLineColor(fFillColors[0]);
      grFakeSys->SetLineWidth(0);
      fFakeGraph.push_back(grFakeSys);
    } else if (fCorrelationGraph) {
      // TODO inplement
      fSysError = new TGraphAsymmErrors();
      double x, y;
      for (int i = 0; i < fCorrelationGraph->GetN(); ++i) {
        fCorrelationGraph->GetPoint(i, x, y);
        fSysError->SetPoint(i, x, y);
        fSysError->SetPointError(
            i, errorwidth, errorwidth,
            y * parameters->Eval(x / (float) fUnitConversionData),
            y * parameters->Eval(x / (float) fUnitConversionData));
      }

      TGraph *grFakeSys = new TGraph();
      SetStyleGraph(grFakeSys, 2, 0);
      grFakeSys->SetFillColor(fFillColors[0]);
      grFakeSys->SetLineColor(fFillColors[0]);
      grFakeSys->SetLineWidth(0);
      fFakeGraph.push_back(grFakeSys);
    } else {
      Warning("DreamData", "For %s set the CF before adding the systematics",
              fName);
    }
  } else {
    Warning("DreamData", "Parameters input missing for %s", fName);
  }
  return;
}

void DreamData::SetSystematics(TH1* parameters, float errorwidth) {
  if (parameters) {
    if (fCorrelationFunction) {
      int nBinsX = fCorrelationFunction->GetNbinsX();
      float minX = fCorrelationFunction->GetXaxis()->GetXmin();
      float maxX = fCorrelationFunction->GetXaxis()->GetXmax();
      fSystematics = new TH1F(Form("Systematics%s", fName),
                              Form("Systematics%s", fName), nBinsX, minX, maxX);
      for (int iBin = 1; iBin < nBinsX; iBin++) {
        const float x = fCorrelationFunction->GetBinCenter(iBin);
        const float y = fCorrelationFunction->GetBinContent(iBin);
        fSystematics->SetBinContent(
            iBin, y * parameters->GetBinContent(parameters->FindBin(x)));
      }
      fSystematics->SetLineWidth(2.0);
      fSysError = new TGraphAsymmErrors();
      for (int i = 0; i < nBinsX; i++) {
        fSysError->SetPoint(i, fCorrelationFunction->GetBinCenter(i + 1),
                            fCorrelationFunction->GetBinContent(i + 1));
        fSysError->SetPointError(i, errorwidth, errorwidth,
                                 fSystematics->GetBinContent(i + 1),
                                 fSystematics->GetBinContent(i + 1));
      }
      TGraph *grFakeSys = new TGraph();
      SetStyleGraph(grFakeSys, 2, 0);
      grFakeSys->SetFillColor(fFillColors[0]);
      grFakeSys->SetLineColor(fFillColors[0]);
      grFakeSys->SetLineWidth(0);
      fFakeGraph.push_back(grFakeSys);
    } else {
      Warning("DreamData", "For %s set the CF before adding the systematics",
              fName);
    }
  } else {
    Warning("DreamData", "Parameters input missing for %s", fName);
  }
  return;
}

void DreamData::SetCorrelatedError(TGraphErrors *grError, int color, float colorAlpha,
                          bool useDefaultColors) {
  fCorrelatedError = grError;
  fCorrelatedError->SetLineColorAlpha((useDefaultColors ? fColors[color] : color), 0);
  fCorrelatedError->SetFillColorAlpha((useDefaultColors ? fColors[color] : color), colorAlpha);
}

void DreamData::SetCorrelatedError(TGraphErrors *grError, int color, int fillstyle,
                          bool useDefaultColors) {
  fCorrelatedError = grError;
  fCorrelatedError->SetLineWidth(1);
  fCorrelatedError->SetLineColor(useDefaultColors ? fColors[color] : 0);
  fCorrelatedError->SetFillColor(useDefaultColors ? fColors[color] : color);
  fCorrelatedError->SetFillStyle(fillstyle);
}

void DreamData::FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower,
                                   TGraph *grUpper, int color, int lineStyle,
                                   double lineWidth, int fillStyle,
                                   bool addtoLegend) {
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
      grFemtoModel->SetPoint(count, x / (float) fUnitConversionCATS,
                             (yAll[2] + yAll[0]) / 2.f);
      grFemtoModel->SetPointError(count++, 0,
                                  (yAll[2] + yAll[0]) / 2.f - yAll[0]);
    }
    grFemtoModel->SetLineColor(fColors[color]);
    grFemtoModel->SetFillColor(fColors[color]);
    grFemtoModel->SetLineWidth(lineWidth);
    grFemtoModel->SetLineStyle(lineStyle);
    if (fillStyle > 0)
      grFemtoModel->SetFillStyle(fillStyle);
    fFemtoModdeled.push_back(grFemtoModel);
    if (addtoLegend) {
      TGraph *grFakeModel = new TGraph();
      grFakeModel->SetLineColor(fColors[color]);
      grFakeModel->SetFillColor(fColors[color]);
      grFakeModel->SetLineWidth(lineWidth * 1.8);
      grFakeModel->SetLineStyle(lineStyle);
      if (fillStyle > 0) {
        grFakeModel->SetFillStyle(fillStyle);
      }
      fFakeGraph.push_back(grFakeModel);
    }
  } else {
    Warning("DreamData", "Set Systematics first for %s", fName);
  }
  return;
}

void DreamData::FemtoModelFitBands(TGraphErrors *grFemtoModel, int color,
                                   int lineStyle, double lineWidth,
                                   int fillStyle, bool addtoLegend,
                                   bool useDefaultColors) {
  grFemtoModel->SetLineColor(useDefaultColors ? fColors[color] : color);
  grFemtoModel->SetFillColor(useDefaultColors ? fColors[color] : color);
  grFemtoModel->SetLineWidth(lineWidth);
  grFemtoModel->SetLineStyle(lineStyle);
  if (fillStyle > 0)
    grFemtoModel->SetFillStyle(fillStyle);
  fFemtoModdeled.push_back(grFemtoModel);
  if (addtoLegend) {
    TGraph *grFakeModel = new TGraph();
    grFakeModel->SetLineColor(useDefaultColors ? fColors[color] : color);
    grFakeModel->SetFillColor(useDefaultColors ? fColors[color] : color);
    grFakeModel->SetLineWidth(lineWidth * 1.8);
    grFakeModel->SetLineStyle(lineStyle);
    if (fillStyle > 0) {
      grFakeModel->SetFillStyle(fillStyle);
    }
    fFakeGraph.push_back(grFakeModel);
  }
}

void DreamData::FemtoModelFitBands(TGraphErrors *grFemtoModel, int color,
                                   float colorAlpha, bool addtoLegend) {
  grFemtoModel->SetLineColorAlpha(color, colorAlpha);
  grFemtoModel->SetFillColorAlpha(color, colorAlpha);
  grFemtoModel->SetLineWidth(0);
  fFemtoModdeled.push_back(grFemtoModel);
  if (addtoLegend) {
    TGraph *grFakeModel = new TGraph();
    grFakeModel->SetLineColorAlpha(color, colorAlpha);
    grFakeModel->SetFillColorAlpha(color, colorAlpha);
    grFakeModel->SetLineWidth(5);
    fFakeGraph.push_back(grFakeModel);
  }
}

void DreamData::FemtoModelDeviations(TGraphErrors* grDeviation, int color, bool useDefaultColors) {
//  SetStyleGraph(grDeviation, 2, color);
  grDeviation->SetMarkerStyle(fMarkers[0]);
  grDeviation->SetMarkerColor(useDefaultColors ? fColors[color] : color);
  grDeviation->SetLineColor(useDefaultColors ? fColors[color] : color);
  grDeviation->SetFillColor(useDefaultColors ? fColors[color] : color);
  fFemtoDeviation.push_back(grDeviation);
}

void DreamData::FemtoModelDeviations(TGraphErrors *grDeviation, int color,
                                     int lineStyle, double lineWidth,
                                     int fillStyle, bool useDefaultColors) {
  grDeviation->SetLineColor(useDefaultColors ? fColors[color] : color);
  grDeviation->SetFillColor(useDefaultColors ? fColors[color] : color);
  grDeviation->SetLineWidth(lineWidth);
  grDeviation->SetLineStyle(lineStyle);
  if (fillStyle > 0)
    grDeviation->SetFillStyle(fillStyle);
  fFemtoDeviation.push_back(grDeviation);
}


void DreamData::SetStyleHisto(TH1 *histo, int marker, int color) {
  if (fMultiHisto) {
    SetStyleMultiHisto(histo, marker, color);
  } else {
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetTitleOffset(1.25);
    histo->SetMarkerSize(1.4);
    histo->SetLineWidth(2);
    histo->SetMarkerStyle(fMarkers[marker]);
    histo->SetMarkerColor(fColors[color]);
    histo->SetLineColor(fColors[color]);
  }
}

void DreamData::SetStyleMultiHisto(TH1 *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.015);
  histo->GetXaxis()->SetTitleSize(0.015);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.015);
  histo->GetYaxis()->SetTitleSize(0.015);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(0.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
}

void DreamData::DrawCorrelationPlot(TPad* c, const int color,
                                    const int systematicsColor,
                                    const float legendTextScale) {
  c->cd();
  TString CFName;
  Color_t markerColor;
  int markerStyle;
  if (fCorrelationFunction) {
    SetStyleHisto(fCorrelationFunction, 2, color);
    fCorrelationFunction->GetXaxis()->SetRangeUser(fXMin, fXMax);
    fCorrelationFunction->GetYaxis()->SetRangeUser(fYMin, fYMax);
    CFName = fCorrelationFunction->GetName();
    markerColor = fCorrelationFunction->GetMarkerColor();
    markerStyle = fCorrelationFunction->GetMarkerStyle();
  } else if (fCorrelationGraph) {
    SetStyleGraph(fCorrelationGraph, 2, color);
    fCorrelationGraph->GetXaxis()->SetRangeUser(fXMin, fXMax);
    fCorrelationGraph->GetYaxis()->SetRangeUser(fYMin, fYMax);
    fCorrelationGraph->SetMarkerSize(1);
    CFName = fCorrelationGraph->GetName();
    markerColor = fCorrelationGraph->GetMarkerColor();
    markerStyle = fCorrelationGraph->GetMarkerStyle();
  }
  fSysError->SetLineColor(kWhite);
  if (!fMultiHisto)
    fSysError->GetYaxis()->SetTitleOffset(1.5);
  fSysError->GetXaxis()->SetRangeUser(fXMin, fXMax);
  fSysError->GetYaxis()->SetRangeUser(fYMin, fYMax);

  if (CFName.Contains("MeV")) {
    fSysError->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  } else {
    fSysError->SetTitle("; #it{k}* (GeV/#it{c}); #it{C}(#it{k}*)");
  }

  if (fDrawAxis) {
    if (fForceAxis) {
      fDummyHist = new TH1F(Form("dummyHist_%.3f", gRandom->Uniform()), fSysError->GetTitle(), 100, fXMin, fXMax);
      fDummyHist->SetMinimum(fYMin);
      fDummyHist->SetMaximum(fYMax);
      fDummyHist->GetXaxis()->SetNdivisions(fSysError->GetXaxis()->GetNdivisions());
      fDummyHist->SetLineColor(kWhite);
      fDummyHist->Draw();
      fSysError->Draw("same");
    } else {
      fSysError->Draw("Ap");
    }
  }
  fBaseLine->Draw("same");

  fLegend = new TLegend(fXMinLegend, fYMinLegend, fXMaxLegend, fYMaxLegend);
//  TLegend *leg = new TLegend(0.5, 0.55, 0.62, 0.875);
  fLegend->SetBorderSize(0);
  fLegend->SetTextFont(42);
  fLegend->SetTextSize(gStyle->GetTextSize() * legendTextScale);
  int legendCounter = 1;
//  leg->AddEntry(fCorrelationFunction, fLegendName[0], "pe");
  fFakeGraph[0]->SetMarkerStyle(markerStyle);
  fFakeGraph[0]->SetMarkerColor(markerColor);
  fFakeGraph[0]->SetFillColorAlpha(systematicsColor, 0.4);
  fLegend->AddEntry(fFakeGraph[0], fLegendName[0], fLegendOption[0]);
//  leg->AddEntry(fBaseLine, "Baseline", "l");
  if (fDrawLegend)
    fLegend->Draw("same");
  for (auto &it : fFemtoModdeled) {
    if (legendCounter < fFakeGraph.size()) {
      fLegend->AddEntry(fFakeGraph[legendCounter], fLegendName[legendCounter],
                        fLegendOption[legendCounter]);
    }
    legendCounter++;
  }
  auto it = fFemtoModdeled.rbegin();
  while (it != fFemtoModdeled.rend()) {
    (*it)->Draw("L3 same");
    it++;
  }
  fSysError->SetFillColorAlpha(systematicsColor, 0.4);
  fSysError->Draw("2 same");
  if(fCorrelationFunction) {
    fCorrelationFunction->DrawCopy("pe same");
  } else if (fCorrelationGraph) {
    fCorrelationGraph->Draw("pez same");
  }
  if(fCorrelatedError) {
    fCorrelatedError->Draw("L3same");
  }
  if (fDrawLegend)
    fLegend->Draw("same");
  if (fInlet) {
    DrawInlet(c);
  }
}

void DreamData::DrawInlet(TPad *c) {
  c->cd();
  TPad *inset_pad = new TPad("insert", "insertPad", fXMinInlet, fYMinInlet,
                             fXMaxInlet, fYMaxInlet);
  inset_pad->SetTopMargin(0.01);
  inset_pad->SetRightMargin(0.05);
  inset_pad->SetBottomMargin(0.28);
  inset_pad->SetLeftMargin(0.28);
  inset_pad->SetFillStyle(4000);
  inset_pad->Draw();
  inset_pad->cd();
  TGraphErrors* SysErrCopy = (TGraphErrors*) fSysError->Clone(
      Form("%s_clone", fSysError->GetName()));
  TH1F* CFCopy = (TH1F*) fCorrelationFunction->Clone(
      Form("%s_Cloned", fCorrelationFunction->GetName()));
  SetStyleHisto(CFCopy, 2, 0);
  CFCopy->GetXaxis()->SetRangeUser(fXMinZoom, fXMaxZoom);
  CFCopy->GetYaxis()->SetRangeUser(fYMinZoom, fYMaxZoom);
  SysErrCopy->GetYaxis()->SetNdivisions(203);
  SysErrCopy->GetXaxis()->SetNdivisions(204);
  SysErrCopy->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  SysErrCopy->GetXaxis()->SetTitleOffset(3.0);
  SysErrCopy->GetYaxis()->CenterTitle(true);
  SysErrCopy->GetYaxis()->SetTitleOffset(1.8);
  SysErrCopy->SetLineColor(kWhite);
  SysErrCopy->Draw("Ap");
  fBaseLine->Draw("same");
//  SysErrCopy->SetTitle(" ; ; ");
  SysErrCopy->GetXaxis()->SetRangeUser(fXMinZoom, fXMaxZoom);
  SysErrCopy->GetYaxis()->SetRangeUser(fYMinZoom, fYMaxZoom);
  for (auto &it : fFemtoModdeled) {
    it->Draw("L3 same");
  }
  SysErrCopy->SetFillColorAlpha(kBlack, 0.4);
  SysErrCopy->Draw("2 same");
  CFCopy->SetMarkerSize(0.6);
  CFCopy->DrawCopy("pe same");
  return;
}

void DreamData::SetStyleGraph(TGraph *histo, int marker, int color) {
  if (fMultiHisto) {
    SetStyleGraphMulti(histo, marker, color);
  } else {
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetXaxis()->SetTitleOffset(1.2);
    histo->GetXaxis()->SetLabelFont(42);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetTitleOffset(1.25);
    histo->SetMarkerSize(1.4);
    histo->SetLineWidth(2);
    histo->SetMarkerStyle(fMarkers[marker]);
    histo->SetMarkerColor(fColors[color]);
    histo->SetLineColor(fColors[color]);
    histo->SetFillColor(fColors[color]);
  }
}

void DreamData::SetStyleGraphMulti(TGraph *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerSize(0.5);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
  histo->SetFillColor(fColors[color]);
}

void DreamData::DrawDeviationPerBin(TPad* c) {
  c->cd();
  TString CFName;
  if(fCorrelationFunction) {
    CFName = fCorrelationFunction->GetName();
  } else if ( fCorrelationGraph) {
    CFName = fCorrelationGraph->GetName();
  }
  TGraphErrors* GraphAxis = (TGraphErrors*) fSysError->Clone("Dummy");
  GraphAxis->Clear();
  GraphAxis->GetYaxis()->SetTitle("n#sigma_{local}");
  GraphAxis->GetYaxis()->CenterTitle(true);
  GraphAxis->GetYaxis()->SetNdivisions(203);
  GraphAxis->GetXaxis()->SetTitleOffset(3.);
  TLine lineOne = TLine(fXMin, 1, fXMax, 1);
  lineOne.SetLineWidth(3);
  lineOne.SetLineColor(kBlack);
  lineOne.SetLineStyle(2);
  for (auto it : fFemtoDeviation) {
    GraphAxis->GetYaxis()->SetRangeUser(it->GetYaxis()->GetXmin(),
                                        it->GetYaxis()->GetXmax());
    GraphAxis->DrawClone("AP");
    it->Draw("L3 same");
    lineOne.DrawLine(fXMin, 0, fXMax, 0);
  }
}

void DreamData::DrawDeviationPerBin(TCanvas* c, float ylow, float yup, float nSigmaMax) {
  const float ylowDraw = ylow;
  const float yupDraw = yup;
  TH1F* histDummy;
  TGraphErrors* grDummy;
  if (fForceAxis) {
    histDummy = (TH1F*) fDummyHist->Clone("Dummy2");
    histDummy->Reset();
    histDummy->GetXaxis()->SetNdivisions(fSysError->GetXaxis()->GetNdivisions());
    histDummy->GetYaxis()->SetTitle("");
    histDummy->GetYaxis()->CenterTitle(true);
    histDummy->GetYaxis()->SetTickLength(0);
    histDummy->GetXaxis()->SetTitleOffset(6.);
    histDummy->GetYaxis()->SetTitleOffset(2.5);
    histDummy->GetYaxis()->SetLabelSize(0.);
    fDummyHist->GetYaxis()->SetTitleOffset(2.5);
  } else {
    grDummy = (TGraphErrors*) fSysError->Clone("Dummy2");
    grDummy->Clear();
    grDummy->GetXaxis()->SetNdivisions(fSysError->GetXaxis()->GetNdivisions());
    grDummy->GetYaxis()->SetTitle("");
    grDummy->GetYaxis()->CenterTitle(true);
    grDummy->GetYaxis()->SetTickLength(0);
    grDummy->GetXaxis()->SetTitleOffset(6.);
    grDummy->GetYaxis()->SetTitleOffset(2.5);
    grDummy->GetYaxis()->SetLabelSize(0.);
    fSysError->GetYaxis()->SetTitleOffset(2.5);
  }
  c->Update();
  const float ticklength =
      (fForceAxis) ?
          histDummy->GetXaxis()->GetTickLength() :
          grDummy->GetXaxis()->GetTickLength();
  const float ticklengthY = fSysError->GetYaxis()->GetTickLength();
  float ymin = 100;
  float ymax = -100.;
  for (auto it : fFemtoDeviation) {
    if(ymin > it->GetYaxis()->GetXmin()) {
      ymin = it->GetYaxis()->GetXmin();
    }
    if(ymax < it->GetYaxis()->GetXmax()) {
      ymax = it->GetYaxis()->GetXmax();
    }
  }
  ymin = std::round(ymin);
  ymax = std::round(ymax);
  if(ymax > std::abs(ymin)) {
    ymin = -1 * ymax;
  } else {
    ymax = std::abs(ymin);
  }
  float xMaxVal;
  if (fForceAxis) {
    histDummy->GetYaxis()->SetRangeUser(ymin - 0.5, ymax + 0.5);
    xMaxVal = fXMax;
  } else {
    grDummy->GetYaxis()->SetRangeUser(ymin - 0.5, ymax + 0.5);
    xMaxVal = grDummy->GetXaxis()->GetXmax();
  }
  TLine lineOne = TLine(fXMin, 1, fXMax, 1);
  lineOne.SetLineWidth(2);
  lineOne.SetLineColor(kBlack);
  lineOne.SetLineStyle(2);

  TGaxis *axis1 = new TGaxis( fXMin, -1. * nSigmaMax, fXMin, nSigmaMax, -1. * nSigmaMax, nSigmaMax, 202,"", ticklengthY); // redraw the axis to have control over the labels
  axis1->SetName("weloveroot");
  axis1->SetLabelFont((fForceAxis) ?
      histDummy->GetXaxis()->GetLabelFont() :
      grDummy->GetXaxis()->GetLabelFont());
  axis1->SetLabelSize((fForceAxis) ?
      histDummy->GetXaxis()->GetLabelSize() :
      grDummy->GetXaxis()->GetLabelSize());
  axis1->SetLabelOffset((fForceAxis) ?
      histDummy->GetXaxis()->GetLabelOffset() :
      grDummy->GetXaxis()->GetLabelOffset());

  TGaxis *axis2 = new TGaxis( fXMax, -1. * nSigmaMax, fXMax, nSigmaMax, -1. * nSigmaMax, nSigmaMax, 202,"+L", ticklengthY); // redraw the axis to have control over the labels
  axis2->SetName("rootisthebest");
  axis2->SetLabelSize(0);

  const float otherPadHeight = 1. - yup;

  const float bottomMargin = 0.4 * 0.5/(yup - ylow);
  const float scalePadHeight = float(fFemtoDeviation.size() - 1. - bottomMargin) / float(fFemtoDeviation.size() - 1);
  const float padHeight = (yup - ylow) / (fFemtoDeviation.size());
  const int nPads = fFemtoDeviation.size();
  TPad* pad[nPads];
  int counter = 0;
  for (auto it : fFemtoDeviation) {
    float currentPadHeight = padHeight;
    if (counter != fFemtoDeviation.size() - 1) {
      currentPadHeight *= scalePadHeight;
      ylow = yup - currentPadHeight;
    } else {
      currentPadHeight *= (1. + bottomMargin);
      ylow = yup - currentPadHeight;
      ylow = (ylow < 0) ? 0. : ylow;
    }
    pad[counter] = new TPad(Form("pad_%i", counter), Form("pad_%i", counter), 0., ylow, 1., yup);
    pad[counter]->cd();
    pad[counter]->SetRightMargin(c->GetRightMargin());
    pad[counter]->SetTopMargin(0.);
    if(counter != fFemtoDeviation.size() - 1) {
      pad[counter]->SetBottomMargin(0.);
    } else {
      pad[counter]->SetBottomMargin(bottomMargin);
    }

    if (fForceAxis) {
      histDummy->GetXaxis()->SetTickLength(
          ticklength * otherPadHeight / currentPadHeight);  // magic scaling - same length as on the large pad!
      histDummy->DrawClone();
    } else {
      grDummy->GetXaxis()->SetTickLength(
          ticklength * otherPadHeight / currentPadHeight);  // magic scaling - same length as on the large pad!
      grDummy->DrawClone("AP");
    }
    axis1->SetTickLength(ticklengthY);
    axis1->SetTickSize(ticklengthY);
    axis1->DrawClone("same");
    axis2->DrawClone("same");
    it->DrawClone("f2");
    lineOne.DrawLine(fXMin, 0, xMaxVal, 0);
    c->cd();
    pad[counter]->DrawClone();
    yup = ylow;
    counter++;
  }
  c->cd();
  TLatex label;
  label.SetNDC(kTRUE);
  label.SetTextFont((fForceAxis) ?
      histDummy->GetYaxis()->GetTitleFont() :
      grDummy->GetYaxis()->GetTitleFont());
  label.SetTextSize((fForceAxis) ?
      histDummy->GetYaxis()->GetTitleSize() :
      grDummy->GetYaxis()->GetTitleSize());
  label.SetTextAngle(90);
  label.DrawLatex(0.04, (yupDraw - ylowDraw) / 2., "#it{n}_{#sigma, local}");
}

void DreamData::DrawLegendExternal(TPad* LegPad) {
  if (!fLegend) {
    Error("DreamData::DrawLegendExternal",
          "No Legend Created yet, call after DrawCorrelationPlot. Exiting \n");
    return;
  }
  LegPad->cd();
  fLegend->Draw();
}
