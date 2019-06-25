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
  VariationmTAnalysis(int nModels, int nData, int nVars);
  virtual ~VariationmTAnalysis();
  void SetSystematic(const char* DataDir);
  void SetVariation(const char* VarDir, int iModel);
  void SetHistName(const char* Histname) {
    fHistname = Histname;
  }
  void SetFileName(const char* Filename) {
    fFileName= Filename;
  }
  void MakeCFPlotsPP();
  void MakeCFPlotsPL();
  void MakeRadPlotsPP();
  void MakeRadPlotsPL(const char* ppFilePath);
  void SetmTAverage(TGraphErrors* mTAvg) {
    fmTAverage = mTAvg;
  }
  TPad* GetFormattedPad(int counter);
 private:
  std::vector<std::vector<VariationAnalysis>> fAnalysis;
  const int fnModel;
  const int fnData;
  const int fnVars;
  std::vector<DreamSystematics> fSystematic;
  const char* fHistname;
  const char* fFileName;
  TGraphErrors* fmTAverage;
  std::vector<TGraphErrors*> fmTRadiusSyst;
  std::vector<TGraphErrors*> fmTRadiusStat;
};
#endif /* GENTLEKITTY_VARIATIONMTANALYSIS_H_ */
