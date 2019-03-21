#ifndef DREAMFUNCTION_DREAMPSYSTEMATICS_H_
#define DREAMFUNCTION_DREAMPSYSTEMATICS_H_

/// \file DreamSystematics.h
/// \brief Definition of a class to compute systematic uncertainties
/// \author A. Mathis

#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"

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
    pL = 3
  };

  DreamSystematics();
  DreamSystematics(Pair pair);
  virtual ~DreamSystematics() = default;

  void SetUpperFitRange(float fitRange) {
    fSystematicFitRangeUp = fitRange;
  }

  void SetDefaultHist(TH1F* histDef) {
    fHistDefault = histDef;
  }
  void SetVarHist(TH1F* histVar) {
    fHistVar.emplace_back(histVar);
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
  ;
  int GetNumberOfVars() const {
    return vars[fParticlePairMode].size();
  }

  TString GetVariation(int count) const {
    return vars[fParticlePairMode][count];
  }

  Pair GetPair() const {
    return fParticlePairMode;
  }
  TString GetPairName() const {
    return pairName[fParticlePairMode];
  }

  TH1F* GetAbsError(TH1F* histDefault, TH1F* histVar) const;
  TH1F* GetErrorBudget(TH1F* histDefault, TH1F* histVar) const;
  TH1F* GetBarlow(TH1F* histDefault, TH1F* histVar) const;

  void EvalSystematics();
  void EvalDifference(std::vector<unsigned int> CountsDefault,
                      std::vector<unsigned int> CountsVar,
                      std::vector<float> *AbsDiff, std::vector<float> *RelDiff);
  TH1F* FillHisto(std::vector<float> Diff, const char* name);
  void EvalDifferenceInPairs();
//  void CountPairs();
  void EvalDifferenceInParticles();
//  void CountParticles();
  void ComputeUncertainty();
  void EvalProtonProton(const int kstar);
  void EvalProtonSigma(const int kstar);
  void EvalProtonXi(const int kstar);
  void EvalProtonLambda(const int kstar);
  void WriteOutput();
  void WriteOutput(TFile* file, std::vector<TH1F*>& histvec,
                   const TString name);
  void FixStyle(TH1F* histCF) const;

 private:
  float fSystematicFitRangeLow;
  float fSystematicFitRangeUp;
  float fFemtoRangeLow;
  float fFemtoRangeUp;
  Pair fParticlePairMode;
  TH1F *fHistDefault;
  TH1F *fHistSystErrAbs;
  TH1F *fHistSystErrRel;
  TGraphErrors* fGrFinalError;
  TF1* fRatio;
  std::vector<TH1F*> fHistVar;
  std::vector<TH1F*> fHistAbsErr;
  std::vector<TH1F*> fHistErrBudget;
  std::vector<TH1F*> fHistBarlow;

  std::vector<unsigned int> fnPairsDefault;
  std::vector<unsigned int> fNPartOneDefault;
  std::vector<unsigned int> fNPartTwoDefault;

  std::vector<unsigned int> fnPairsVar;
  std::vector<unsigned int> fNPartOneVariations;
  std::vector<unsigned int> fNPartTwoVariations;

  std::vector<float> fnPairsAbsDiff;
  std::vector<float> fnPartOneAbsDiff;
  std::vector<float> fnPartTwoAbsDiff;

  std::vector<float> fnPairsRelDiff;
  std::vector<float> fnPartOneRelDiff;
  std::vector<float> fnPartTwoRelDiff;

  TH1F* fHistPairsAbsDiff;
  TH1F* fHistPartOneAbsDiff;
  TH1F* fHistPartTwoAbsDiff;

  TH1F* fHistPairsRelDiff;
  TH1F* fHistPartOneRelDiff;
  TH1F* fHistPartTwoRelDiff;

  const std::vector<TString> ppVariations = { { "Proton #it{p}_{T} down",
      "Proton #it{p}_{T} up", "Proton #eta up", "Proton #eta down",
      "Proton n#sigma up", "Proton n#sigma down", "Proton FilterBit",
      "Proton TPC cluster down", "Proton TPC cluster up", "Proton CPR 0.012",
      "Proton CPR 0.014" } };

  const std::vector<TString> pSigma0Variations = { { "Proton #it{p}_{T} down",
      "Proton #it{p}_{T} up", "Proton #eta down", "Proton #eta up",
      "Proton n#sigma up", "Proton n#sigma down", "Proton FilterBit",
      "Proton TPC cluster down", "Proton TPC cluster up",
      "Lambda #it{p}_{T} down", "Lambda #it{p}_{T} up", "Lambda CPA up",
      "Lambda CPA down", "Lambda Daug n#sigma up", "Lambda Daug n#sigma down",
      "Lambda Armenteros-Podolandski", "Lambda Daug TPC ncls up",
      "Lambda Daug TPC ncls down", "Lambda Daug #eta down", "Lambda mass",
      "Photon #eta down", "Photon #it{p}_{T} down", "Photon #it{p}_{T} up",
      "Photon Daug TPC ncls finable up", "Photon Daug n#sigma down",
      "Photon Daug n#sigma up", "Photon q_{T} 1-D up", "Photon q_{T} 2-D down",
      "Photon #Psi_{Pair} 1-D up", "Photon #Psi_{Pair} 2-D down",
      "Photon CPA down", "Photon CPA up"} };

  const std::vector<TString> pXiVariations = {
      { "Proton #it{p}_{T} down", "Proton #it{p}_{T} up", "Proton #eta up",
          "Proton #eta down", "Proton n#sigma up", "Proton n#sigma down",
          "Proton FilterBit", "Proton TPC cluster down",
          "Proton TPC cluster up", "Casc Daug d_{track} up",
          "Casc Daug d_{track} down", "Casc Bach d_{Track, PV} down",
          "Casc Bach d_{Track, PV} up", "Casc CPA up 1", "Casc CPA up 2",
          "Casc Transverse Radius down", "Casc Transverse Radius up",
          "Casc V0 Daug d_{track} up 1", "Casc V0 Daug d_{track} up 2",
          "Casc V0 CPA down", "Casc V0 CPA up",
          "Casc V0 Transverse Radius down", "Casc V0 Transverse Radius up",
          "Casc V0 d_{V0, PV} up", "Casc V0 d_{V0, PV} down",
          "Casc V0 Daug d_{Track, PV} down", "Casc V0 Daug d_{Track, PV} up",
          "Casc Track #eta up", "Casc Track #eta down", "Casc Track n#sigma up",
          "Casc Track n#sigma down", "Casc #it{p}_{T} up 1",
          "Casc #it{p}_{T} up 2" } };

  const std::vector<TString> pLVariations = { { "Proton #it{p}_{T} down",
      "Proton #it{p}_{T} up", "Proton #eta up", "Proton #eta down",
      "Proton n#sigma up", "Proton n#sigma down", "Proton FilterBit",
      "Proton TPC cluster down", "Proton TPC cluster up", "V0 #it{p}_{T} down",
      "V0 #it{p}_{T} up", "V0 CPA up", "V0 Track n#sigma up",
      "V0 Track TPC cluster up", "V0 Track #eta up", "V0 Track #eta down",
      "V0 Track d_{track} up", "V0 Track d_{Track, PV}" } };

  const std::vector<std::vector<TString>> vars = { { ppVariations,
      pSigma0Variations, pXiVariations, pLVariations } };
  const std::vector<TString> pairName = { { "pp", "pSigma0", "pXi", "pL" } };
};

inline
void DreamSystematics::FixStyle(TH1F* histCF) const {
  histCF->Sumw2();
  histCF->GetXaxis()->SetRangeUser(fSystematicFitRangeLow,
                                   fSystematicFitRangeUp);
  DreamPlot::SetStyleHisto(histCF);
  histCF->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  histCF->GetYaxis()->SetTitle("C(#it{k}*)");
}

#endif // DREAMFUNCTION_DREAMPSYSTEMATICS_H_

