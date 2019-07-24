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
  void SetLegData(const char* DataName, const char* DataOption) {
    fDataName= DataName;
    fDataOption = DataOption;
  }
  void SetLegModel(const char* ModelName, const char* ModelOption, int Color) {
    fModelName.push_back(ModelName);
    fModelOption.push_back(ModelOption);
    fColor.push_back(Color);
  }
  void SourceName(const char* sourceName) {
    fSourceName = sourceName;
  }
  void MakeCFPlotsSingleBand();
  void MakeOnePanelPlots();
  void MakeCFPlotsPL();
  void SetmTAverage(TGraphErrors* mTAvg) {
    fmTAverage = mTAvg;
  }
  TPad* GetFormattedPad(int counter);
  void SetmTBins(std::vector<float> mTBins) {
    fmTBins = mTBins;
  }
  void StoreRadvsmT(const char* fileName);
  void SetTextXMin(float texMin) {
    fTextXMin = texMin;
  }
  void SetPlottingRange(float xMin, float xMax) {
    fXmin = xMin;
    fXmax = xMax;
  }
 private:
  std::vector<std::vector<VariationAnalysis>> fAnalysis;
  const int fnModel;
  const int fnData;
  const int fnVars;
  std::vector<float> fmTBins;

  const char* fHistname;
  const char* fFileName;
  const char* fDataName;
  const char* fDataOption;
  std::vector<const char*> fModelName;
  std::vector<const char*> fModelOption;
  const char* fSourceName;
  std::vector<int> fColor;
  float fTextXMin;
  float fXmin;
  float fXmax;

  std::vector<DreamSystematics> fSystematic;
  TGraphErrors* fmTAverage;
  std::vector<TGraphErrors*> fmTRadiusSyst;
  std::vector<TGraphErrors*> fmTRadiusStat;
};
#endif /* GENTLEKITTY_VARIATIONMTANALYSIS_H_ */
