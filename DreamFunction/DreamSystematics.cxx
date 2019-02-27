/// \file DreamSystematics.cxx
/// \brief Implementation of a class to compute systematic uncertainties
/// \author A. Mathis

#include "DreamSystematics.h"
#include "TFile.h"
#include <cmath>

DreamSystematics::DreamSystematics()
    : fSystematicFitRangeLow(0.f),
      fSystematicFitRangeUp(600.f),
      fFemtoRangeLow(0.f),
      fFemtoRangeUp(200.f),
      fParticlePairMode(DreamSystematics::pp),
      fHistDefault(nullptr),
      fHistSystErrAbs(nullptr),
      fHistSystErrRel(nullptr),
      fnPairsDefault(0),
      fRatio(nullptr),
      fHistVar(),
      fHistAbsErr(),
      fHistErrBudget(),
      fHistBarlow(),
      fnPairsVar() {
}

DreamSystematics::DreamSystematics(Pair pair)
    : fSystematicFitRangeLow(0.f),
      fSystematicFitRangeUp(600.f),
      fFemtoRangeLow(0.f),
      fFemtoRangeUp(200.f),
      fParticlePairMode(pair),
      fHistDefault(nullptr),
      fHistSystErrAbs(nullptr),
      fHistSystErrRel(nullptr),
      fnPairsDefault(0),
      fRatio(nullptr),
      fHistVar(),
      fHistAbsErr(),
      fHistErrBudget(),
      fHistBarlow(),
      fnPairsVar() {
}

TH1F* DreamSystematics::GetAbsError(TH1F* histDefault, TH1F* histVar) const {
  auto histNew = (TH1F*) histVar->Clone(Form("%s_AbsErr", histVar->GetName()));
  histNew->SetTitle(Form("AbsErr %s", histVar->GetTitle()));
  histNew->GetYaxis()->SetTitle("|CF_{default} - CF_{var}|");

  for (int i = 1; i < histVar->GetNbinsX() + 1; ++i) {
    histNew->SetBinContent(
        i, std::abs(histDefault->GetBinContent(i) - histVar->GetBinContent(i)));
    histNew->SetBinError(
        i,
        std::sqrt(
            std::pow(histDefault->GetBinError(i), 2)
                + std::pow(histVar->GetBinError(i), 2)));
  }
  return histNew;
}

TH1F* DreamSystematics::GetErrorBudget(TH1F* histDefault, TH1F* histVar) const {
  auto histNew = (TH1F*) histVar->Clone(
      Form("%s_ErrBudget", histVar->GetName()));
  histNew->SetTitle(Form("ErrBudget %s", histVar->GetTitle()));
  histNew->GetYaxis()->SetTitle("|1  - CF_{var} / CF_{default}|");

  float binContDef, binContVar, binErrDef, binErrVar;
  for (int i = 1; i < histVar->GetNbinsX() + 1; ++i) {
    binContDef = histDefault->GetBinContent(i);
    binContVar = histVar->GetBinContent(i);
    binErrDef = histDefault->GetBinError(i);
    binErrVar = histVar->GetBinError(i);

    histNew->SetBinContent(i, std::abs(1 - binContVar / binContDef));
    histNew->SetBinError(
        i,
        std::sqrt(
            std::pow(binErrVar, 2) / std::pow(binContDef, 2)
                + std::pow(binErrDef, 2) * std::pow(binContVar, 2)
                    / std::pow(binContDef, 4)

                    ));
  }
  return histNew;
}

TH1F* DreamSystematics::GetBarlow(TH1F* histDefault, TH1F* histVar) const {
  auto histNew = (TH1F*) histVar->Clone(Form("%s_Barlow", histVar->GetName()));
  histNew->SetTitle(Form("Barlow %s", histVar->GetTitle()));
  histNew->GetYaxis()->SetTitle("n#sigma");

  float binContDef, binContVar, binErrDef, binErrVar, statErrVariation;
  for (int i = 1; i < histVar->GetNbinsX() + 1; ++i) {
    binContDef = histDefault->GetBinContent(i);
    binErrDef = histDefault->GetBinError(i);
    binContVar = histVar->GetBinContent(i);
    binErrVar = histVar->GetBinError(i);
    statErrVariation = std::sqrt(binErrVar * binErrVar + binErrDef * binErrDef);

    histNew->SetBinContent(
        i, std::abs(binContVar - binContDef) / statErrVariation);
    histNew->SetBinError(i, 0);
  }
  return histNew;
}

void DreamSystematics::EvalSystematics() {
  FixStyle(fHistDefault);
  fSystematicFitRangeLow = fHistDefault->GetBinCenter(1)
      - fHistDefault->GetXaxis()->GetBinWidth(1) / 2.f;
  fHistDefault->GetXaxis()->SetRangeUser(fSystematicFitRangeLow,
                                         fSystematicFitRangeUp);
  int iVar = 0;
  for (auto &it : fHistVar) {
    it->SetTitle(GetVariation(iVar));
    FixStyle(it);
    fHistAbsErr.emplace_back(GetAbsError(fHistDefault, it));
    fHistErrBudget.emplace_back(GetErrorBudget(fHistDefault, it));
    fHistBarlow.emplace_back(GetBarlow(fHistDefault, it));
    ++iVar;
  }

  ComputeUncertainty();
}

void DreamSystematics::ComputeUncertainty() {

  fHistSystErrAbs = (TH1F*) fHistDefault->Clone(
      Form("%s_SystErrAbs", fHistDefault->GetName()));
  fHistSystErrAbs->GetYaxis()->SetTitle("Syst. error");
  fHistSystErrAbs->Clear();

  fHistSystErrRel = (TH1F*) fHistDefault->Clone(
      Form("%s_SystErrRel", fHistDefault->GetName()));
  fHistSystErrRel->GetYaxis()->SetTitle("Rel. syst. error");
  fHistSystErrRel->Clear();

  const int nBins = fHistDefault->GetXaxis()->FindBin(fSystematicFitRangeUp);
  for (int ikstar = 0; ikstar < nBins; ++ikstar) {
    switch (fParticlePairMode) {
      case Pair::pp:
        EvalProtonProton(ikstar);
        break;
      case Pair::pSigma0:
        EvalProtonSigma(ikstar);
        break;
    }
  }

  fRatio = new TF1(Form("Ratio_%s", fHistDefault->GetName()), "pol2", 0.,
                   fSystematicFitRangeUp);
  fHistSystErrRel->Fit(fRatio, "R", 0, fSystematicFitRangeUp);
}

void DreamSystematics::EvalProtonProton(const int kstar) {
  std::vector<float> addSyst;

  // pT
  addSyst.push_back(
      (fHistAbsErr[0]->GetBinContent(kstar)
          + fHistAbsErr[1]->GetBinContent(kstar)) / 2.);
  // Eta
  addSyst.push_back(
      (fHistAbsErr[2]->GetBinContent(kstar)
          + fHistAbsErr[3]->GetBinContent(kstar)) / 2.);
  // nSigma
  addSyst.push_back(
      (fHistAbsErr[4]->GetBinContent(kstar)
          + fHistAbsErr[5]->GetBinContent(kstar)) / 2.);
  // Filter Bit 96
  addSyst.push_back(fHistAbsErr[6]->GetBinContent(kstar));
  // TPC Cls
  addSyst.push_back(
      (fHistAbsErr[7]->GetBinContent(kstar)
          + fHistAbsErr[8]->GetBinContent(kstar)) / 2.);

  double sysErrTotal = 0.;
  for (size_t iAdd = 0; iAdd < addSyst.size(); ++iAdd) {
    sysErrTotal += (addSyst[iAdd] * addSyst[iAdd]);
  }
  sysErrTotal = std::sqrt(sysErrTotal);

  fHistSystErrAbs->SetBinContent(kstar, sysErrTotal);
  fHistSystErrRel->SetBinContent(
      kstar, sysErrTotal / fHistDefault->GetBinContent(kstar));
}

void DreamSystematics::EvalProtonSigma(const int kstar) {

}

void DreamSystematics::WriteOutput() {
  auto file = new TFile(Form("Systematics_%s.root", GetPairName().Data()),
                        "recreate");
  file->mkdir("Raw");
  file->cd("Raw");
  fHistDefault->Write("histDefault");
  int iVar = 0;
  for (auto &it : fHistVar) {
    it->Write(Form("histVar_%i", iVar++));
  }
  file->cd();
  file->mkdir("AbsErr");
  file->cd("AbsErr");
  iVar = 0;
  for (auto &it : fHistAbsErr) {
    it->Write(Form("histAbsErr_%i", iVar++));
  }

  file->cd();
  file->mkdir("ErrBudget");
  file->cd("ErrBudget");
  iVar = 0;
  for (auto &it : fHistErrBudget) {
    it->Write(Form("histErrBudget_%i", iVar++));
  }

  file->cd();
  file->mkdir("Barlow");
  file->cd("Barlow");
  iVar = 0;
  for (auto &it : fHistBarlow) {
    it->Write(Form("histBarlow_%i", iVar++));
  }

  file->cd();
  fHistSystErrAbs->Write("SystErrAbs");
  fHistSystErrRel->Write("SystErrRel");
  fRatio->Write("SystError");

  file->Close();

}
