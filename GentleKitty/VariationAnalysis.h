/*
 * DrawPP.h
 *
 *  Created on: May 8, 2019
 *      Author: schmollweger
 *      all units in MeV
 */
#ifndef GENTLEKITTY_VARIATIONANALYSIS_H_
#define GENTLEKITTY_VARIATIONANALYSIS_H_
#include "TROOT.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TGraphErrors.h"
#include <vector>

class VariationAnalysis {
 public:
  VariationAnalysis(const char* histName, const int nDataVars,
                    const int nFitVars);
  virtual ~VariationAnalysis();
  void ReadFitFile(TString FileName);
  TGraphErrors* EvaluateCurves(TNtuple* tuple, TGraph* RefGraph);
  void EvalRadius();
  void SetRadiusRanges(int nRadBins, float radMin, float radMax) {
    fnRadBins = nRadBins;
    fRadMin = radMin;
    fRadMax = radMax;
  }
  TH1F* GetCorrelationFunctio(int i) const {
    return fCk.at(i);
  }
  TH1F* GetRadDist() const {
    return fRadiusDist;
  }
  TGraphErrors* GetModel() {
    return fModel;
  }
  float GetRadStatErr() const {
    return fRadStat;
  }
 private:
  TFile* fInFile;
  const char* fHistname;
  const int fnDataVars;
  const int fnFitVars;
  int fnRadBins;
  float fRadMin;
  float fRadMax;
  float fRadStat;
  std::vector<TH1F*> fCk;
  TGraphErrors* fModel;
  TH1F* fRadiusDist;
};

#endif /* GENTLEKITTY_VARIATIONANALYSIS_H_ */
