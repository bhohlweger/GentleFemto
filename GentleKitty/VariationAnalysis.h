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
  TGraphErrors* ModelFitBands();
  TGraphErrors* EvaluateCurves(TNtuple* tuple, const int nBins, const int kMin,
                               const int dkStar);
  TH1F* GetCorrelationFunctio(int i) const {
    return fCk.at(i);
  }
 private:
  TFile* fInFile;
  const char* fHistname;
  const int fnDataVars;
  const int fnFitVars;
  float fnkMin;
  int fnModelBins;
  int fdkstar;
  std::vector<TH1F*> fCk;
  TNtuple* fFits;
};

#endif /* GENTLEKITTY_VARIATIONANALYSIS_H_ */
