/*
 * DrawPP.h
 *
 *  Created on: May 8, 2019
 *      Author: schmollweger
 *      all units in MeV
 */
#ifndef GENTLEKITTY_VARIATIONANALYSISPAL_H_
#define GENTLEKITTY_VARIATIONANALYSISPAL_H_
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include <vector>
#include "DreamPlot.h"

class VariationAnalysispAL {
 public:
  VariationAnalysispAL(const char* histName);
  virtual ~VariationAnalysispAL();
  void ReadFitFile(TString FileName);
  TGraphErrors* EvaluateCurves(TNtuple* tuple, TGraph* RefGraph);
  TGraphErrors* DeviationByBin(TH1F* RefHist, TGraphErrors* model);
  TH1F* GetCorrelationFunction(int i) const {
    return fCk.at(i);
  }
  TH1D* GetRadDist() const {
    return fRadiusDist;
  }
  TGraphErrors* GetModel() {
    return fModel;
  }
  TGraphErrors* GetDeviationByBin() {
    return fDeviationByBin;
  }
  float GetRadStatErr() const {
    return fRadStat;
  }
  float GetRadMean() const {
    return fRadMean;
  }
  float GetRadSystDown() const {
    return fRadSystDown;
  }
  float GetRadSystUp() const {
    return fRadSystUp;
  }
  void AppendAndCut(TCut anotherCut) {
    fSelector = fSelector && anotherCut;
  }
  void AppendOrCut(TCut anotherCut) {
    fSelector = fSelector || anotherCut;
  }
 private:
  TFile* fInFile;
  TCut fSelector;
  const char* fHistname;
  int fnDataVarsStart;
  int fnDataVarsEnd;
  int fnFitVarsStart;
  int fnFitVarsEnd;
  float fRadMean;
  float fRadSystUp;
  float fRadSystDown;
  float fRadStat;
  std::vector<TH1F*> fCk;
  TGraphErrors* fModel;
  TGraphErrors* fDeviationByBin;
  TH1D* fRadiusDist;
  TH1D* fRadiusErrDist;
};

#endif /* GENTLEKITTY_VARIATIONANALYSIS_H_ */
