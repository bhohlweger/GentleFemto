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
      fFemtoRangeUp(0.2),
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
      fnPairsVar(),
      fnPairsAbsDiff(),
      fnPairsRelDiff(),
      fNPartOneDefault(),
      fNPartOneVariations(),
      fNPartTwoDefault(),
      fNPartTwoVariations(),
      fHistPairsAbsDiff(nullptr),
      fHistPairsRelDiff(nullptr),
      fHistPartOneAbsDiff(nullptr),
      fHistPartOneRelDiff(nullptr),
      fHistPartTwoAbsDiff(nullptr),
      fHistPartTwoRelDiff(nullptr) {
}

DreamSystematics::DreamSystematics(Pair pair)
    : fSystematicFitRangeLow(0.f),
      fSystematicFitRangeUp(600.f),
      fFemtoRangeLow(0.f),
      fFemtoRangeUp(0.2),
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
      fnPairsVar(),
      fnPairsAbsDiff(),
      fnPairsRelDiff(),
      fNPartOneDefault(0),
      fNPartOneVariations(),
      fNPartTwoDefault(0),
      fNPartTwoVariations(),
      fHistPairsAbsDiff(nullptr),
      fHistPairsRelDiff(nullptr),
      fHistPartOneAbsDiff(nullptr),
      fHistPartOneRelDiff(nullptr),
      fHistPartTwoAbsDiff(nullptr),
      fHistPartTwoRelDiff(nullptr) {
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
void DreamSystematics::EvalDifference(std::vector<unsigned int> CountsDefault,
                                      std::vector<unsigned int> CountsVar,
                                      std::vector<float> *AbsDiff,
                                      std::vector<float> *RelDiff) {
  size_t nVars = CountsVar.size();
  for (unsigned int iVar = 0; iVar < nVars; ++iVar) {
    AbsDiff->push_back(
        std::abs(
            (float) ((float) CountsDefault[iVar] - (float) CountsVar[iVar])));
    RelDiff->push_back(
        CountsDefault[iVar] > 0 ?
            std::abs(
                (float) (1 - (CountsVar[iVar] / (float) CountsDefault[iVar]))) :
            0);
  }
  return;
}

TH1F* DreamSystematics::FillHisto(std::vector<float> Diff, const char* name) {
  std::cout << "" << Diff.size() << '\t'
            << TString::Format("%s%s", name, GetPairName().Data()) << std::endl;
  TH1F* outHisto = new TH1F(TString::Format("%s%s", name, GetPairName().Data()),
                            TString::Format("%s%s", name, GetPairName().Data()),
                            Diff.size(), 0, Diff.size());
  for (unsigned int iBin = 1; iBin < Diff.size(); ++iBin) {
    outHisto->GetXaxis()->SetBinLabel(iBin, GetVariation(iBin - 1).Data());
    outHisto->SetBinContent(iBin, Diff[iBin - 1]);
  }
  return outHisto;
}

void DreamSystematics::EvalDifferenceInPairs() {
  if (fnPairsDefault.size() == 0 || fnPairsVar.size()) {
    std::cout
        << "DreamSystematics::EvalDifferenceInPairs() : no variations set \n";
  } else {
    EvalDifference(fnPairsDefault, fnPairsVar, &fnPairsAbsDiff,
                   &fnPairsRelDiff);
    fHistPairsAbsDiff = FillHisto(fnPairsAbsDiff, "AbsDiffPair");
    fHistPairsRelDiff = FillHisto(fnPairsRelDiff, "RelDiffPair");
  }
}

void DreamSystematics::EvalDifferenceInParticles() {
  if (fNPartOneDefault.size() == 0 || fNPartOneVariations.size() == 0) {
    std::cout
        << "DreamSystematics::EvalDifferenceInParticles() : default or var not set for part one \n";
  } else {
    EvalDifference(fNPartOneDefault, fNPartOneVariations, &fnPartOneAbsDiff,
                   &fnPartOneRelDiff);
    fHistPartOneAbsDiff = FillHisto(fnPartOneAbsDiff, "AbsDiffPartOne");
    fHistPartOneRelDiff = FillHisto(fnPartOneRelDiff, "RelDiffPartOne");
  }
  if (fNPartTwoDefault.size() == 0 || fNPartTwoVariations.size() == 0) {
    std::cout
        << "DreamSystematics::EvalDifferenceInParticles() : default or var not set for part two \n";
  } else {
    EvalDifference(fNPartTwoDefault, fNPartTwoVariations, &fnPartTwoAbsDiff,
                   &fnPartTwoRelDiff);
    fHistPartTwoAbsDiff = FillHisto(fnPartTwoAbsDiff, "AbsDiffPartTwo");
    fHistPartTwoRelDiff = FillHisto(fnPartTwoRelDiff, "RelDiffPartTwo");
  }
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
      case Pair::pXi:
        EvalProtonXi(ikstar);
        break;
      default:
        std::cout << "Non implemented Particle mode \n";
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
  std::vector<float> addSyst;

  // === PROTON ===
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

  // === Lambda ===
  // pT
  addSyst.push_back(
      (fHistAbsErr[9]->GetBinContent(kstar)
          + fHistAbsErr[10]->GetBinContent(kstar)) / 2.);
  // CPA
  addSyst.push_back(
      (fHistAbsErr[11]->GetBinContent(kstar)
          + fHistAbsErr[12]->GetBinContent(kstar)) / 2.);
  // TPC nSigma PID
  addSyst.push_back(
      (fHistAbsErr[13]->GetBinContent(kstar)
          + fHistAbsErr[14]->GetBinContent(kstar)) / 2.);
  // Armenteros-Podolandski
  addSyst.push_back(fHistAbsErr[15]->GetBinContent(kstar));
  // TPC Cls
  addSyst.push_back(
      (fHistAbsErr[16]->GetBinContent(kstar)
          + fHistAbsErr[17]->GetBinContent(kstar)) / 2.);
  // Eta
  addSyst.push_back(fHistAbsErr[18]->GetBinContent(kstar));
  // Mass window
  addSyst.push_back(fHistAbsErr[19]->GetBinContent(kstar));

  // === Photon ===
  // Eta
  addSyst.push_back(fHistAbsErr[20]->GetBinContent(kstar));
  // pT
  addSyst.push_back(
      (fHistAbsErr[21]->GetBinContent(kstar)
          + fHistAbsErr[22]->GetBinContent(kstar)) / 2.);
  // TPC Cls
  addSyst.push_back(fHistAbsErr[23]->GetBinContent(kstar));
  // TPC nSigma PID
  addSyst.push_back(
      (fHistAbsErr[24]->GetBinContent(kstar)
          + fHistAbsErr[25]->GetBinContent(kstar)) / 2.);
  // Armenteros-Podolandski
  addSyst.push_back(
      (fHistAbsErr[26]->GetBinContent(kstar)
          + fHistAbsErr[27]->GetBinContent(kstar)) / 2.);
  // Psi_pair
  addSyst.push_back(
      (fHistAbsErr[28]->GetBinContent(kstar)
          + fHistAbsErr[29]->GetBinContent(kstar)) / 2.);
  // CPA
  addSyst.push_back(
      (fHistAbsErr[30]->GetBinContent(kstar)
          + fHistAbsErr[31]->GetBinContent(kstar)) / 2.);
  // DCA
  addSyst.push_back(
      (fHistAbsErr[32]->GetBinContent(kstar)
          + fHistAbsErr[33]->GetBinContent(kstar)) / 2.);

  double sysErrTotal = 0.;
  for (size_t iAdd = 0; iAdd < addSyst.size(); ++iAdd) {
    sysErrTotal += (addSyst[iAdd] * addSyst[iAdd]);
  }
  sysErrTotal = std::sqrt(sysErrTotal);

  fHistSystErrAbs->SetBinContent(kstar, sysErrTotal);
  fHistSystErrRel->SetBinContent(
      kstar, sysErrTotal / fHistDefault->GetBinContent(kstar));

}

void DreamSystematics::EvalProtonXi(const int kstar) {
  std::vector<float> addSyst;
  //pT
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
  //Daughter DCA to Casc Vtx
  addSyst.push_back(
      (fHistAbsErr[9]->GetBinContent(kstar)
          + fHistAbsErr[10]->GetBinContent(kstar)) / 2.);
  //Bach DCA to PV
  addSyst.push_back(
      (fHistAbsErr[11]->GetBinContent(kstar)
          + fHistAbsErr[12]->GetBinContent(kstar)) / 2.);
  //Xi CPA
  addSyst.push_back(
      (fHistAbsErr[13]->GetBinContent(kstar)
          + fHistAbsErr[14]->GetBinContent(kstar)) / 2.);
  //Xi Transverse Radius
  addSyst.push_back(
      (fHistAbsErr[15]->GetBinContent(kstar)
          + fHistAbsErr[16]->GetBinContent(kstar)) / 2.);
  //Daugh DCA to V0 Vtx
  addSyst.push_back(
      (fHistAbsErr[17]->GetBinContent(kstar)
          + fHistAbsErr[18]->GetBinContent(kstar)) / 2.);
  //V0 CPA
  addSyst.push_back(
      (fHistAbsErr[19]->GetBinContent(kstar)
          + fHistAbsErr[20]->GetBinContent(kstar)) / 2.);
  //V0 Transverse Radius
  addSyst.push_back(
      (fHistAbsErr[21]->GetBinContent(kstar)
          + fHistAbsErr[22]->GetBinContent(kstar)) / 2.);
  //V0 DCA to PV
  addSyst.push_back(
      (fHistAbsErr[23]->GetBinContent(kstar)
          + fHistAbsErr[24]->GetBinContent(kstar)) / 2.);
  //Daugh DCA to PV
  addSyst.push_back(
      (fHistAbsErr[25]->GetBinContent(kstar)
          + fHistAbsErr[26]->GetBinContent(kstar)) / 2.);
  //Xi Tracks Eta
  addSyst.push_back(
      (fHistAbsErr[27]->GetBinContent(kstar)
          + fHistAbsErr[28]->GetBinContent(kstar)) / 2.);
  //Xi Tracks PID
  addSyst.push_back(
      (fHistAbsErr[29]->GetBinContent(kstar)
          + fHistAbsErr[30]->GetBinContent(kstar)) / 2.);
  // Xi Pt
  addSyst.push_back(fHistAbsErr[31]->GetBinContent(kstar));

  double sysErrTotal = 0.;
  for (size_t iAdd = 0; iAdd < addSyst.size(); ++iAdd) {
    sysErrTotal += (addSyst[iAdd] * addSyst[iAdd]);
  }
  sysErrTotal = std::sqrt(sysErrTotal);

  fHistSystErrAbs->SetBinContent(kstar, sysErrTotal);
  fHistSystErrRel->SetBinContent(
      kstar, sysErrTotal / fHistDefault->GetBinContent(kstar));
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
  if (fnPairsDefault.size() > 0) {
    file->mkdir("nPairs");
    file->cd("nPairs");
    if (fHistPairsAbsDiff)
      fHistPairsAbsDiff->Write();
    if (fHistPairsRelDiff)
      fHistPairsRelDiff->Write();
  }
  if (fNPartOneDefault.size() > 0 || fNPartTwoDefault.size() > 0) {
    file->mkdir("nParticles");
    file->cd("nParticles");
    if (fHistPartOneAbsDiff)
      fHistPartOneAbsDiff->Write();
    if (fHistPartOneRelDiff)
      fHistPartOneRelDiff->Write();
    if (fHistPartTwoAbsDiff)
      fHistPartTwoAbsDiff->Write();
    if (fHistPartTwoRelDiff)
      fHistPartTwoRelDiff->Write();
  }
  file->cd();
  fHistSystErrAbs->Write("SystErrAbs");
  fHistSystErrRel->Write("SystErrRel");
  fRatio->Write("SystError");

  file->Close();

}
