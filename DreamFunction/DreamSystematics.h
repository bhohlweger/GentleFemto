#ifndef DREAMFUNCTION_DREAMPSYSTEMATICS_H_
#define DREAMFUNCTION_DREAMPSYSTEMATICS_H_

/// \file DreamSystematics.h
/// \brief Definition of a class to compute systematic uncertainties
/// \author A. Mathis

#include <vector>
#include "TH1F.h"
#include "TString.h"

#include "DreamPlot.h"

/// \class DreamSystematics
/// This class takes care of computing the systematic uncertainties for different particle pairs
class DreamSystematics {
 public:
  enum Pair {
    pp = 0,       ///< Proton-Proton correlation function
    pSigma0 = 1,       ///< Proton-Sigma0 correlation function
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
  void ComputeUncertainty();
  void EvalProtonProton(const int kstar);
  void EvalProtonSigma(const int kstar);
  void WriteOutput();
  void DrawAllCF();
  void FixStyle(TH1F* histCF) const;

 private:
  float fSystematicFitRangeLow;
  float fSystematicFitRangeUp;
  Pair fParticlePairMode;
  TH1F *fHistDefault;
  TH1F *fHistSystErrAbs;
  TH1F *fHistSystErrRel;
  TF1* fRatio;
  std::vector<TH1F*> fHistVar;
  std::vector<TH1F*> fHistAbsErr;
  std::vector<TH1F*> fHistErrBudget;
  std::vector<TH1F*> fHistBarlow;

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

  const std::vector<std::vector<TString>> vars = { { ppVariations,
      pSigma0Variations } };
  const std::vector<TString> pairName = { { "pp", "pSigma0" } };
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

