/*
 * VariationmTAnalysis.h
 *
 *  Created on: Jun 18, 2019
 *      Author: schmollweger
 */

#ifndef GENTLEKITTY_VARIATIONMTANALYSIS_H_
#define GENTLEKITTY_VARIATIONMTANALYSIS_H_

#include "VariationAnalysis.h"
#include "DreamSystematics.h"
#include <vector>

class VariationmTAnalysis {
 public:
  VariationmTAnalysis();
  virtual ~VariationmTAnalysis();
  void SetSystematic(const char* DataDir);
  void SetVariation(const char* VarDir);
  void SetHistName(const char* Histname) {
    fHistname = Histname;
  }
  void MakePlots();
  void SetmTAverage(TGraphErrors* mTAvg) {
    fmTAverage = mTAvg;
  }
  TPad* GetFormattedPad(int counter);
 private:
  std::vector<VariationAnalysis> fAnalysis;
  std::vector<DreamSystematics> fSystematic;
  const char* fHistname;
  TGraphErrors* fmTAverage;
  TGraphErrors* fmTRadiusSyst;
  TGraphErrors* fmTRadiusStat;
};

#endif /* GENTLEKITTY_VARIATIONMTANALYSIS_H_ */
