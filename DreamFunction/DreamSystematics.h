#ifndef DREAMFUNCTION_DREAMPSYSTEMATICS_H_
#define DREAMFUNCTION_DREAMPSYSTEMATICS_H_

/// \file DreamSystematics.h
/// \brief Definition of a class to compute systematic uncertainties
/// \author A. Mathis

#include <vector>
#include "TH1F.h"
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
    pXi = 2
  };

  DreamSystematics();
  DreamSystematics(Pair pair);
  virtual ~DreamSystematics() = default;

  void SetUpperFitRange(float fitRange) {
    fSystematicFitRangeUp = fitRange;
  }
  void SetDefaultPair(DreamCF* CFDef, TString CFName) {
    SetDefaultHist(CFDef->FindCorrelationFunction(CFName));
    fnPairsDefault = CFDef->GetFemtoPairs(fFemtoRangeLow, fFemtoRangeUp);
  }
  void SetVarPair(DreamCF* CFVar, TString CFName) {
    SetVarHist(CFVar->FindCorrelationFunction(CFName));
    fnPairsVar.push_back(CFVar->GetFemtoPairs(fFemtoRangeLow, fFemtoRangeUp));
  }
  void SetDefaultHist(TH1F* histDef) {
    fHistDefault = histDef;
  }
  void SetVarHist(TH1F* histVar) {
    fHistVar.emplace_back(histVar);
  }
  void SetNDefaultParticles(unsigned int nPart1, unsigned int nPart2) {
    fNPartOneDefault = nPart1;
    fNPartTwoDefault = nPart2;
  };
  void SetNVariationParticles(unsigned int nPart1, unsigned int nPart2) {
    fNPartOneVariations.push_back(nPart1);
    fNPartTwoVariations.push_back(nPart2);
  };
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
  void EvalDifferenceInPairs();
  void CountPairs();
  void EvalDifferenceInParticles();
  void CountParticles();
  void ComputeUncertainty();
  void EvalProtonProton(const int kstar);
  void EvalProtonSigma(const int kstar);
  void EvalProtonXi(const int kstar);
  void WriteOutput();
  void DrawAllCF();
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
  unsigned int fnPairsDefault;
  TF1* fRatio;
  std::vector<TH1F*> fHistVar;
  std::vector<TH1F*> fHistAbsErr;
  std::vector<TH1F*> fHistErrBudget;
  std::vector<TH1F*> fHistBarlow;
  std::vector<unsigned int> fnPairsVar;
  std::vector<unsigned int> fnPairsAbsDiff;
  std::vector<float> fnPairsRelDiff;
  unsigned int fNPartOneDefault;
  std::vector<unsigned int> fNPartOneVariations;
  unsigned int fNPartTwoDefault;
  std::vector<unsigned int> fNPartTwoVariations;
  std::vector<unsigned int> fnPartOneAbsDiff;
  std::vector<float> fnPartOneRelDiff;
  std::vector<unsigned int> fnPartTwoAbsDiff;
  std::vector<float> fnPartTwoRelDiff;
  TH1F* fHistPairsAbsDiff;
  TH1F* fHistPairsRelDiff;
  TH1F* fHistPartOneAbsDiff;
  TH1F* fHistPartOneRelDiff;
  TH1F* fHistPartTwoAbsDiff;
  TH1F* fHistPartTwoRelDiff;
  const std::vector<TString> ppVariations = { { "Proton #it{p}_{T} down",
      "Proton #it{p}_{T} up", "Proton #eta up", "Proton #eta down",
      "Proton n#sigma up", "Proton n#sigma down", "Proton FilterBit",
      "Proton TPC cluster down", "Proton TPC cluster up" } };

  const std::vector<TString> pSigma0Variations = { {
      "Lambda PileUp one daughter", "Lambda #it{p}_{T} 0.24",
      "Lambda #it{p}_{T} 0.36", "Lambda CPA 0.995", "Lambda CPA 0.98",
      "Lambda Daug n#sigma 3", "Lambda Daug n#sigma 3 only dE/dx",
      "Lambda Daug n#sigma 6", "Lambda Daug n#sigma 6 only dE/dx",
      "Lambda only dE/dx", "Lambda Daug TPC ncls 80", "Lambda Daug TPC ncls 60",
      "Lambda Daug #eta 0.8", "Lambda Daug #eta 0.75",
      "Lambda Daug d_{track} 1.2", "Lambda Daug d_{track} 0.9",
      "Lambda Daug d_{Track, PV} 0.06", "Lambda Daug d_{Track, PV} 0.04",
      "Lambda K0 Rejection", "Lambda Lambda selection 8",
      "Lambda radius up 120 decay vtx 120", "Lambda radius up 80",
      "Lambda radius low 0", "Lambda radius low 5", "Lambda radius low 5",
      "Lambda decay vtx 80", "Lambda decay vtx 120", "Photon #eta 0.8",
      "Photon #eta 0.75", "Photon r 0-180", "Photon r 10-180",
      "Photon #it{p}_{T,e} = 0, #it{p}_{T} = 0",
      "Photon #it{p}_{T,e} = 0.15, #it{p}_{T} = 0.02",
      "Photon #it{p}_{T,e} = 0.05, #it{p}_{T} = 0.15",
      "Photon Daug TPC ncls finable 0.6", "Photon Daug n#sigma 10",
      "Photon Daug n#sigma 3", "Photon pKpi rej.", "Photon q_{T} < 0.1 1-D",
      "Photon q_{T} < 0.02 2-D", "Photon #chi^{2} < 100",
      "Photon #chi^{2} < 10", "Photon #Psi_{Pair} < 0.2 1-D",
      "Photon #Psi_{Pair} < 0.1 2-D", "Photon CPA 0.98", "Photon CPA 0.995",
      "Photon Shared e", "Photon DCA_R < 5", "Photon DCA_Z < 5",
      "Photon Pile-up", "Proton #it{p}_{T} down", "Proton #it{p}_{T} up",
      "Proton #eta up", "Proton #eta down", "Proton n#sigma up",
      "Proton n#sigma down", "Proton FilterBit", "Proton TPC cluster down",
      "Proton TPC cluster up" } };

  const std::vector<TString> pXiVariations = { { "Proton #it{p}_{T} down",
      "Proton #it{p}_{T} up", "Proton #eta up", "Proton #eta down",
      "Proton n#sigma up", "Proton n#sigma down", "Proton FilterBit",
      "Proton TPC cluster down", "Proton TPC cluster up",
      "Casc Daug d_{track} up", "Casc Daug d_{track} down",
      "Casc Bach d_{Track, PV} down", "Casc Bach d_{Track, PV} up",
      "Casc CPA up 1", "Casc CPA up 2", "Casc Transverse Radius down",
      "Casc Transverse Radius up", "Casc V0 Daug d_{track} up 1",
      "Casc V0 Daug d_{track} up 2", "Casc V0 CPA down", "Casc V0 CPA up",
      "Casc V0 Transverse Radius down", "Casc V0 Transverse Radius up",
      "Casc V0 d_{V0, PV} up", "Casc V0 d_{V0, PV} down",
      "Casc V0 Daug d_{Track, PV} down", "Casc V0 Daug d_{Track, PV} up",
      "Casc Track #eta up", "Casc Track #eta down", "Casc Track n#sigma up",
      "Casc Track n#sigma down", "Casc #it{p}_{T} down" } };

  const std::vector<std::vector<TString>> vars = { { ppVariations,
      pSigma0Variations, pXiVariations } };
  const std::vector<TString> pairName = { { "pp", "pSigma0", "pXi" } };
};

inline void DreamSystematics::DrawAllCF() {
  auto c = new TCanvas();
  DreamPlot::SetStyleHisto(fHistDefault, 20, 0);
  fHistDefault->SetMarkerSize(0.5);
  fHistDefault->Draw();
  fHistDefault->GetXaxis()->SetRangeUser(fSystematicFitRangeLow,
                                         fSystematicFitRangeUp);
  fHistDefault->GetYaxis()->SetRangeUser(0, 4);
  int iCount = 0;
  for (auto &it : fHistVar) {
    DreamPlot::SetStyleHisto(it, 20 + iCount, ++iCount);
    it->SetMarkerSize(0.5);
    it->Draw("same");
  }
  c->Print(TString::Format("CF_var_%s.pdf", GetPairName().Data()));
}

inline
void DreamSystematics::FixStyle(TH1F* histCF) const {
  histCF->Sumw2();
  histCF->GetXaxis()->SetRangeUser(fSystematicFitRangeLow,
                                   fSystematicFitRangeUp);
  histCF->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
  histCF->GetYaxis()->SetTitle("C(#it{k}*)");
}

#endif // DREAMFUNCTION_DREAMPSYSTEMATICS_H_

