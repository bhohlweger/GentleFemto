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
#include <vector>

class DreamData {
 public:
  DreamData(const char* particlePair);
  void SetSystematics(TF1* parameters, float errorwidth);
  void SetCorrelationFunction(TH1F* CF) {
    fCorrelationFunction = CF;
  }
  ;
  const TH1F* GetCorrelationFunction() { return fCorrelationFunction; }
  void SetCorrelationFunctionSimulation(TH1F* CF) {
    fCorrelationFunctionSimulation = CF;
  }
  ;
  void SetBaseline(float BL0, float BL1) {
    fBaseLine->SetParameter(0, BL0);
    fBaseLine->SetParameter(1, BL1);
  }
  ;
  void FemtoModelFitBands(TGraph *grMedian1, TGraph *grLower, TGraph *grUpper,
                          int color, int lineStyle, double lineWidth,
                          int fillStyle, bool addtoLegend = true);
  void FemtoModelFitBands(TGraphErrors *gr, int color, int lineStyle,
                          double lineWidth, int fillStyle, bool addtoLegend);
  void SetStyleHisto(TH1 *histo, int marker, int color);
  void DrawCorrelationPlot(TCanvas* c, const int color = 0);
  void SetRangePlotting(float xMin, float xMax, float yMin, float yMax) {
    fXMin = xMin;
    fXMax = xMax;
    fYMin = yMin;
    fYMax = yMax;
  }
  ;
  void SetInletRangePlotting(float xRangeMin, float xRangeMax, float yRangeMin, float yRangeMax) {
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
  void SetLegendCoordinates(float xMin, float yMin, float xMax, float yMax) {
    fXMinLegend = xMin;
    fXMaxLegend = xMax;
    fYMinLegend = yMin;
    fYMaxLegend = yMax;
  }
  ;
  void SetNDivisions(int nDiv) {
    if (fSysError) {
      fSysError->GetXaxis()->SetNdivisions(nDiv);
    } else {
      std::cout << "No sys err for " << fName;
    }
  }
  ;
  void SetLegendName(const char* name, const char* option) {
    fLegendName.push_back(name);
    fLegendOption.push_back(option);
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
  void SetStyleGraph(TGraph *histo, int marker, int color);
  void DrawInlet(TCanvas *c);
  int GetNumberOfModels() const { return fFemtoModdeled.size(); }
  virtual ~DreamData();
  const char* fName;
  TH1F* fCorrelationFunction;
  TH1F* fCorrelationFunctionSimulation;
  TH1F* fSystematics;
  TGraphErrors* fSysError;
  TF1* fBaseLine;
  float fXMin;
  float fXMax;
  float fYMin;
  float fYMax;
  bool  fInlet;
  float fXMinZoom;
  float fXMaxZoom;
  float fYMinZoom;
  float fYMaxZoom;
  float fXMinInlet;
  float fXMaxInlet;
  float fYMinInlet;
  float fYMaxInlet;
  float fXMinLegend;
  float fXMaxLegend;
  float fYMinLegend;
  float fYMaxLegend;
  int fUnitConversionData;
  int fUnitConversionCATS;
  std::vector<const char*> fLegendName;
  std::vector<const char*> fLegendOption;
  std::vector<TGraphErrors*> fFemtoModdeled;
  std::vector<TGraph*> fFakeGraph;
  std::vector<int> fFillColors;
  std::vector<int> fColors;
  std::vector<int> fMarkers;
};

#endif /* DREAMFUNCTION_DREAMDATA_H_ */
