/*
 * DreamData.h
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#ifndef DREAMFUNCTION_DREAMDATA_H_
#define DREAMFUNCTION_DREAMDATA_H_
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>

class DreamData {
 public:
  DreamData(const char* particlePair);
  void SetSystematics(TF1* parameters, float errorwidth);
  void SetSystematics(TH1* parameters, float errorwidth);
  void SetCorrelationFunction(TH1F* CF) {
    fCorrelationFunction = CF;
  }
  ;
  void SetCorrelationGraph(TGraphAsymmErrors* graph) {
    fCorrelationGraph = graph;
  }
  const TH1F* GetCorrelationFunction() {
    return fCorrelationFunction;
  }
  void SetCorrelationFunctionSimulation(TH1F* CF) {
    fCorrelationFunctionSimulation = CF;
  }
  ;
  void SetBaseline(float BL0, float BL1) {
    fBaseLine->SetParameter(0, BL0);
    fBaseLine->SetParameter(1, BL1);
  }
  ;
  void SetCorrelatedError(TGraphErrors *geError, int color, float colorAlpha,
                          bool useDefaultColors = true);
  void SetCorrelatedError(TGraphErrors *geError, int color, int fillstyle,
                          bool useDefaultColors = true);
  void FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper,
                          int color, int lineStyle, double lineWidth,
                          int fillStyle, bool addtoLegend = true);
  void FemtoModelFitBands(TGraphErrors *gr, int color, int lineStyle,
                          double lineWidth, int fillStyle, bool addtoLegend,
                          bool useDefaultColors = true);
  void FemtoModelFitBands(TGraphErrors *grFemtoModel, int color,
                          float colorAlpha, bool addtoLegend);
  void FemtoModelFitBands(TGraphErrors *graphMod, int color, int lineStyle,
                          double lineWidth, double opacity, bool addtoLegend,
                          bool useDefaultColors = true);
  void FemtoModelDeviations(TGraphErrors* grDeviation, int color, bool useDefaultColors = true);
  void FemtoModelDeviations(TGraphErrors *grDeviation, int color, int lineStyle,
                            double lineWidth, int fillStyle,
                            bool useDefaultColors = true);
  void SetStyleHisto(TH1 *histo, int marker, int color);
  void SetStyleMultiHisto(TH1 *histo, int marker, int color);
  void DrawCorrelationPlot(TPad* c, const int color = 0,
                           const int systematicsColor = kBlack,
                           const float legendTextScale = 0.96, const float markersize = 1.4);
  void DrawDeviationPerBin(TPad* c);
  void DrawDeviationPerBin(TCanvas* c, float ylow, float yup, float nSigmaMax);
  void SetDrawAxis(bool drawAxis) {
    fDrawAxis = drawAxis;
  }
  void SetForceAxisRanges(bool force) { fForceAxis = force; }
  void SetRangePlotting(float xMin, float xMax, float yMin, float yMax) {
    fXMin = xMin;
    fXMax = xMax;
    fYMin = yMin;
    fYMax = yMax;
  }
  ;
  void SetInletRangePlotting(float xRangeMin, float xRangeMax, float yRangeMin,
                             float yRangeMax) {
    fInlet = true;
    fXMinZoom = xRangeMin;
    fXMaxZoom = xRangeMax;
    fYMinZoom = yRangeMin;
    fYMaxZoom = yRangeMax;
  }
  ;
  void SetInletCoordinates(float xMin, float yMin, float xMax, float yMax) {
    fInlet = true;
    fXMinInlet = xMin;
    fXMaxInlet = xMax;
    fYMinInlet = yMin;
    fYMaxInlet = yMax;
  }
  ;
  void SetLegendCoordinates(float xMin, float yMin, float xMax, float yMax,
                            bool drawLegend = true) {
    fXMinLegend = xMin;
    fXMaxLegend = xMax;
    fYMinLegend = yMin;
    fYMaxLegend = yMax;
    fDrawLegend = drawLegend;
  }
  ;
  void SetAxisOffsetInlet(float xaxis, float yaxis) {
    fXAxisOffsetInlet = xaxis;
    fYAxisOffsetInlet = yaxis;
  }
  void SetTextSizeLegend(float size) {
    fTextSizeInlet = size;
  }
  void DrawLegendExternal(TPad* LegPad);
  void SetNDivisions(int nDiv) {
    if (fSysError) {
      fSysError->GetXaxis()->SetNdivisions(nDiv);
    } else {
      Warning("DreamData", "No sys err for %s", fName);
    }
  }
  ;
  void SetLegendName(const char* name, const char* option) {
    fLegendName.push_back(name);
    fLegendOption.push_back(option);
  }
  ;
  void SetLegendName(std::vector<const char*> name, std::vector<const char*> option) {
    fLegendName = name;
    fLegendOption = option;
  }
  ;
  void SetUnitConversionData(int unit) {
    fUnitConversionData = unit;
  }
  ;
  void SetUnitConversionCATS(int unit) {
    fUnitConversionCATS = unit;
  }
  ;
  void SetMultiHisto(bool multiHisto) {
    fMultiHisto = multiHisto;
  }
  void SetStyleGraph(TGraph *histo, int marker, int color);
  void SetStyleGraphMulti(TGraph *histo, int marker, int color);
  void DrawInlet(TPad *c);
  int GetNumberOfModels() const {
    return fFemtoModdeled.size();
  }
  virtual ~DreamData();
  const char* fName;
  TH1F* fCorrelationFunction;
  TH1F* fCorrelationFunctionSimulation;
  TH1F* fSystematics;
  TGraphAsymmErrors* fCorrelationGraph;
  TGraphAsymmErrors* fSysError;
  TGraphErrors *fCorrelatedError;
  TH1F* fDummyHist;
  TF1* fBaseLine;
  bool fDrawAxis;
  bool fForceAxis;
  float fXMin;
  float fXMax;
  float fYMin;
  float fYMax;
  bool fInlet;
  float fXMinZoom;
  float fXMaxZoom;
  float fYMinZoom;
  float fYMaxZoom;
  float fXMinInlet;
  float fXMaxInlet;
  float fYMinInlet;
  float fYMaxInlet;
  float fXAxisOffsetInlet;
  float fYAxisOffsetInlet;
  float fTextSizeInlet;
  TLegend* fLegend;
  float fXMinLegend;
  float fXMaxLegend;
  float fYMinLegend;
  float fYMaxLegend;
  bool fDrawLegend;
  int fUnitConversionData;
  int fUnitConversionCATS;
  bool fMultiHisto;
  std::vector<const char*> fLegendName;
  std::vector<const char*> fLegendOption;
  std::vector<TGraphErrors*> fFemtoModdeled;
  std::vector<TGraphErrors*> fFemtoDeviation;
  std::vector<TGraph*> fFakeGraph;
  std::vector<int> fFillColors;
  std::vector<int> fColors;
  std::vector<int> fMarkers;
};

#endif /* DREAMFUNCTION_DREAMDATA_H_ */
