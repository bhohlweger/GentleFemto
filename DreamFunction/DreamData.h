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
  void SetSystematics(TF1* parameters, int UnitConv, float errorwidth);
  void SetCorrelationFunction(TH1F* CF) {
    fCorrelationFunction = CF;
  }
  ;
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
                          int UnitConv, int color);
  void SetStyleHisto(TH1 *histo, int marker, int color);
  void DrawCorrelationPlot(TCanvas* c);
  void SetRangePlotting(float xMin, float xMax, float yMin, float yMax) {
    fXMin = xMin;
    fXMax = xMax;
    fYMin = yMin;
    fYMax = yMax;
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
  void SetLegendName(const char* name) {
    fLegendName.push_back(name);
  }
  ;
  void SetStyleGraph(TGraph *histo, int marker, int color);
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
  std::vector<const char*> fLegendName;
  std::vector<TGraphErrors*> fFemtoModdeled;
  std::vector<TGraph*> fFakeGraph;
  std::vector<int> fFillColors;
  std::vector<int> fColors;
  std::vector<int> fMarkers;
};

#endif /* DREAMFUNCTION_DREAMDATA_H_ */
