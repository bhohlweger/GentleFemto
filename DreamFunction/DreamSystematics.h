#ifndef DREAMFUNCTION_DREAMPSYSTEMATICS_H_
#define DREAMFUNCTION_DREAMPSYSTEMATICS_H_

/// \file DreamSystematics.h
/// \brief Definition of a class to compute systematic uncertainties
/// \author A. Mathis

#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TNtuple.h"

#include "DreamPlot.h"
#include "DreamCF.h"
/// \class DreamSystematics
/// This class takes care of computing the systematic uncertainties for different particle pairs
class DreamSystematics {
 public:
  enum Pair {
    pp = 0,       ///< Proton-Proton correlation function
    pSigma0 = 1,       ///< Proton-Sigma0 correlation function
    pXi = 2,
    pL = 3,
    pAp = 4,       ///< Proton-AntiProton correlation function
    pAL = 5,       ///< Proton-AntiLambda correlation function
    LAL = 6,       ///< Lambda-AntiLambda correlation function
    pXiNorm = 7,
    pXiLam = 8,
    pXiRes = 9,
    pD = 10
  };

  enum Estimator {
    Uniform = 0,
    StdDev = 1,
    WeightedStdDev = 2
  };

  DreamSystematics();
  DreamSystematics(Pair pair);
  virtual ~DreamSystematics() = default;

  void SetUpperFitRange(float fitRange) {
    fSystematicFitRangeUp = fitRange;
  }
  void SetBarlowUpperRange(float upperRange) {
    fBarlowUpperRange = upperRange;
  }

  void SetDefaultHist(TH1F* histDef) {
    if (!histDef) {
      Warning("SetDefaultHist", "No Histogram exists, nothing set \n");
    } else {
      fHistDefault = histDef;
    }
  }
  void SetVarHist(TH1F* histVar) {
    if (!histVar) {
      Warning("SetDefaultHist", "No Histogram exists, nothing set \n");
    } else {
      fHistVar.emplace_back(histVar);
    }
  }
  void SetDefaultHist(DreamCF* CFDef, const char* CFName) {
    SetDefaultHist(CFDef->FindCorrelationFunction(CFName));
  }
  void SetVarHist(DreamCF* CFVar, const char* CFName) {
    SetVarHist(CFVar->FindCorrelationFunction(CFName));
  }
  void SetPair(unsigned int nDef, unsigned int nVar) {
    fnPairsDefault.push_back(nDef);
    fnPairsVar.push_back(nVar);
  }
  void SetParticles(unsigned int nDef1, unsigned int nDef2, unsigned int nPart1,
                    unsigned int nPart2) {
    fNPartOneDefault.push_back(nDef1);
    fNPartTwoDefault.push_back(nDef2);
    fNPartOneVariations.push_back(nPart1);
    fNPartTwoVariations.push_back(nPart2);
  }
  void SetPurity(float purDef1, float purDef2, float purPart1, float purPart2) {
    fPurityOneDefault.push_back(purDef1);
    fPurityTwoDefault.push_back(purDef2);
    fPurityOneVariations.push_back(purPart1);
    fPurityTwoVariations.push_back(purPart2);
  }
  int GetNumberOfVars() const {
    return vars[fParticlePairMode];
  }

  Pair GetPair() const {
    return fParticlePairMode;
  }
  TString GetPairName() const {
    return pairName[fParticlePairMode];
  }

  void SetEstimator(Estimator est) {
    fErrorEstimator = est;
  }

  TH1F* GetAbsError(TH1F* histDefault, TH1F* histVar) const;
  TH1F* GetErrorBudget(TH1F* histDefault, TH1F* histVar) const;
  TH1F* GetBarlow(TH1F* histDefault, TH1F* histVar);
  void FillTuple(TH1F* histVar);

  void EvalSystematics();
  void EvalSystematicsBBar(int doRebin);
  template<typename T>
  void EvalDifference(std::vector<T> &CountsDefault, std::vector<T> &CountsVar,
                      std::vector<float> &AbsDiff, std::vector<float> &RelDiff);
  TH1F* FillHisto(std::vector<float> Diff, const char* name);
  void EvalDifferenceInPairs();
  void EvalDifferenceInParticles();
  void EvalDifferenceInPurity();
  void ComputeUncertainty();
  void ComputeUncertaintyBBar(int doRebin);
  void WriteOutput(const char* outname = "");
  template<typename T>
  void WriteOutput(TFile* file, std::vector<T*>& histvec, const TString name);
  void FixStyle(TH1F* histCF) const;
  TF1* GetSystematicError() const {
    return fRatio;
  }
  ;
  TH1F* GetDefault() const {
    return fHistDefault;
  }
 private:
  float fSystematicFitRangeLow;
  float fSystematicFitRangeUp;
  float fBarlowUpperRange;
  Pair fParticlePairMode;
  Estimator fErrorEstimator;
  TH1F *fHistDefault;
  TH1F *fHistSystErrAbs;
  TH1F *fHistSystErrRel;
  TGraphErrors* fGrFinalError;
  TF1* fRatio;
  std::vector<TH1F*> fHistVar;
  std::vector<TH1F*> fHistAbsErr;
  std::vector<TH1F*> fHistErrBudget;
  std::vector<TH1F*> fHistBarlow;
  std::vector<TH1D*> fHistKstar;
  TNtuple *fCutTuple;
  std::vector<TString> fBarlowLabel;

  std::vector<unsigned int> fnPairsDefault;
  std::vector<unsigned int> fNPartOneDefault;
  std::vector<unsigned int> fNPartTwoDefault;

  std::vector<unsigned int> fnPairsVar;
  std::vector<unsigned int> fNPartOneVariations;
  std::vector<unsigned int> fNPartTwoVariations;

  std::vector<float> fPurityDefault;
  std::vector<float> fPurityOneDefault;
  std::vector<float> fPurityTwoDefault;

  std::vector<float> fPurityVar;
  std::vector<float> fPurityOneVariations;
  std::vector<float> fPurityTwoVariations;

  std::vector<float> fnPairsAbsDiff;
  std::vector<float> fnPartOneAbsDiff;
  std::vector<float> fnPartTwoAbsDiff;
  std::vector<float> fPurityOneAbsDiff;
  std::vector<float> fPurityTwoAbsDiff;

  std::vector<float> fnPairsRelDiff;
  std::vector<float> fnPartOneRelDiff;
  std::vector<float> fnPartTwoRelDiff;
  std::vector<float> fPurityOneRelDiff;
  std::vector<float> fPurityTwoRelDiff;

  TH1F* fHistPairsAbsDiff;
  TH1F* fHistPartOneAbsDiff;
  TH1F* fHistPartTwoAbsDiff;
  TH1F* fHistPurityOneAbsDiff;
  TH1F* fHistPurityTwoAbsDiff;

  TH1F* fHistPairsRelDiff;
  TH1F* fHistPartOneRelDiff;
  TH1F* fHistPartTwoRelDiff;
  TH1F* fHistPurityOneRelDiff;
  TH1F* fHistPurityTwoRelDiff;

  const int ppVariations = 44;
  const int pSigma0Variations = 25;
  const int pXiVariations = 44;
  const int pLVariations = 44;

  const std::vector<int> vars = { { ppVariations, pSigma0Variations,
      pXiVariations, pLVariations } };
  const std::vector<TString> pairName = { { "pp", "pSigma0", "pXi", "pL", "pAp",
      "pAL", "LAL", "pXiNorm", "pXiLam", "pXiRes", "pD" } };
};

inline
void DreamSystematics::FixStyle(TH1F* histCF) const {
  histCF->Sumw2();
  histCF->GetXaxis()->SetRangeUser(fSystematicFitRangeLow,
                                   fSystematicFitRangeUp);
  DreamPlot::SetStyleHisto(histCF, 22, 1, -1);
  histCF->GetXaxis()->SetTitle("#it{k}* (MeV/#it{c})");
  histCF->GetYaxis()->SetTitle("C(#it{k}*)");
}

inline void DreamSystematics::FillTuple(TH1F *histVar) {

  for (int i = 1; i < histVar->GetNbinsX() + 1; ++i) {
    fCutTuple->Fill(histVar->GetBinCenter(i), histVar->GetBinContent(i),
                    histVar->GetBinError(i));
  }
}

#endif // DREAMFUNCTION_DREAMPSYSTEMATICS_H_
