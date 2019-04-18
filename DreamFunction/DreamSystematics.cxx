/// \file DreamSystematics.cxx
/// \brief Implementation of a class to compute systematic uncertainties
/// \author A. Mathis

#include "DreamSystematics.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"
#include <cmath>

DreamSystematics::DreamSystematics()
    : fSystematicFitRangeLow(0.f),
      fSystematicFitRangeUp(600.f),
      fBarlowUpperRange(600.f),
      fParticlePairMode(DreamSystematics::pp),
      fHistDefault(nullptr),
      fHistSystErrAbs(nullptr),
      fHistSystErrRel(nullptr),
      fGrFinalError(nullptr),
      fnPairsDefault(0),
      fRatio(nullptr),
      fHistVar(),
      fHistAbsErr(),
      fHistErrBudget(),
      fHistBarlow(),
      fBarlowLabel(),
      fnPairsVar(),
      fnPairsAbsDiff(),
      fnPairsRelDiff(),
      fNPartOneDefault(),
      fNPartOneVariations(),
      fNPartTwoDefault(),
      fNPartTwoVariations(),
      fPurityDefault(),
      fPurityOneDefault(),
      fPurityTwoDefault(),
      fPurityVar(),
      fPurityOneVariations(),
      fPurityTwoVariations(),
      fHistPairsAbsDiff(nullptr),
      fHistPartOneAbsDiff(nullptr),
      fHistPartTwoAbsDiff(nullptr),
      fHistPurityOneAbsDiff(nullptr),
      fHistPurityTwoAbsDiff(nullptr),
      fHistPairsRelDiff(nullptr),
      fHistPartOneRelDiff(nullptr),
      fHistPartTwoRelDiff(nullptr),
      fHistPurityOneRelDiff(nullptr),
      fHistPurityTwoRelDiff(nullptr) {
  DreamPlot::SetStyle();
}

DreamSystematics::DreamSystematics(Pair pair)
    : fSystematicFitRangeLow(0.f),
      fSystematicFitRangeUp(600.f),
      fBarlowUpperRange(600.f),
      fParticlePairMode(pair),
      fHistDefault(nullptr),
      fHistSystErrAbs(nullptr),
      fHistSystErrRel(nullptr),
      fGrFinalError(nullptr),
      fnPairsDefault(0),
      fRatio(nullptr),
      fHistVar(),
      fHistAbsErr(),
      fHistErrBudget(),
      fHistBarlow(),
      fBarlowLabel(),
      fnPairsVar(),
      fnPairsAbsDiff(),
      fnPairsRelDiff(),
      fNPartOneDefault(0),
      fNPartOneVariations(),
      fNPartTwoDefault(0),
      fNPartTwoVariations(),
      fPurityDefault(),
      fPurityOneDefault(),
      fPurityTwoDefault(),
      fPurityVar(),
      fPurityOneVariations(),
      fPurityTwoVariations(),
      fHistPairsAbsDiff(nullptr),
      fHistPartOneAbsDiff(nullptr),
      fHistPartTwoAbsDiff(nullptr),
      fHistPurityOneAbsDiff(nullptr),
      fHistPurityTwoAbsDiff(nullptr),
      fHistPairsRelDiff(nullptr),
      fHistPartOneRelDiff(nullptr),
      fHistPartTwoRelDiff(nullptr),
      fHistPurityOneRelDiff(nullptr),
      fHistPurityTwoRelDiff(nullptr) {
  DreamPlot::SetStyle();
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
  auto fit = new TF1(Form("%s_fit", histNew->GetName()), "pol0",
                     fSystematicFitRangeLow, fSystematicFitRangeUp);
  fit->SetLineColor(kGreen + 2);
  histNew->Fit(fit, "QR");
  return histNew;
}

TH1F* DreamSystematics::GetBarlow(TH1F* histDefault, TH1F* histVar) {
  auto histNew = (TH1F*) histVar->Clone(Form("%s_Barlow", histVar->GetName()));
  histNew->SetTitle(Form("Barlow %s", histVar->GetTitle()));
  histNew->GetYaxis()->SetTitle("n#sigma");
  int nSigBel2 = 0;
  int nSigBel3 = 0;
  int nSigAb3 = 0;
  int sigBins = 0;
  float binContDef, binContVar, binErrDef, binErrVar, statErrVariation;
  for (int i = 1; i < histVar->GetNbinsX() + 1; ++i) {
    binContDef = histDefault->GetBinContent(i);
    binErrDef = histDefault->GetBinError(i);
    binContVar = histVar->GetBinContent(i);
    binErrVar = histVar->GetBinError(i);
    statErrVariation = std::sqrt(
        TMath::Abs(binErrVar * binErrVar - binErrDef * binErrDef));
    float nSigma =0;
    if (statErrVariation > 1e-6) {
      nSigma = std::abs(binContVar - binContDef) / statErrVariation;
    } else {
      nSigma = 99;
    }
    histNew->SetBinContent(i, nSigma);
    histNew->SetBinError(i, 0);
  }
  for (int i = 1; i < histNew->GetNbinsX() + 1; ++i) {
    const float nSigma = histNew->GetBinContent(i);
    if(histNew->GetBinCenter(i) > fBarlowUpperRange) continue;
    if (nSigma < 2) {
      nSigBel2++;
    } else if (nSigma >= 2 && nSigma < 3) {
      nSigBel3++;
    } else if (nSigma >= 3) {
      nSigAb3++;
    }
    sigBins++;
  }

  TString label = TString::Format(
      "#splitline{#splitline{Below 2: %.2f (%%)}{Between 2 - 3: %.2f (%%)}}{Above 3: %.2f (%%)}",
      nSigBel2*100 / (float) sigBins, nSigBel3*100 / (float) sigBins,
      nSigAb3*100 / (float) sigBins);
  fBarlowLabel.push_back(label);
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

template<typename T>
void DreamSystematics::EvalDifference(std::vector<T> &CountsDefault,
                                      std::vector<T> &CountsVar,
                                      std::vector<float> &AbsDiff,
                                      std::vector<float> &RelDiff) {
  for (size_t iVar = 0; iVar < CountsVar.size(); ++iVar) {
    float def = CountsDefault[iVar];
    float var = CountsVar[iVar];
    AbsDiff.push_back(std::abs(def - var));
    RelDiff.push_back((def > 0) ? std::abs(1 - var / def) * 100.f : 0);
  }
  return;
}

TH1F* DreamSystematics::FillHisto(std::vector<float> Diff, const char* name) {
  TH1F* outHisto = new TH1F(TString::Format("%s%s", name, GetPairName().Data()),
                            TString::Format("%s%s", name, GetPairName().Data()),
                            Diff.size(), 0, Diff.size());
  for (unsigned int iBin = 1; iBin <= Diff.size(); ++iBin) {
    outHisto->GetXaxis()->SetBinLabel(iBin, GetVariation(iBin - 1).Data());
    outHisto->SetBinContent(iBin, Diff[iBin - 1]);
  }
  DreamPlot::SetStyleHisto(outHisto);
  outHisto->GetXaxis()->LabelsOption("v");
  return outHisto;
}

void DreamSystematics::EvalDifferenceInPairs() {
  if (fnPairsDefault.size() == 0 || fnPairsVar.size() == 0) {
    Error("DreamSystematics",
          "DreamSystematics::EvalDifferenceInPairs() : no variations set");
  } else {
    EvalDifference(fnPairsDefault, fnPairsVar, fnPairsAbsDiff, fnPairsRelDiff);
    fHistPairsAbsDiff = FillHisto(fnPairsAbsDiff, "AbsDiffPair");
    fHistPairsAbsDiff->GetYaxis()->SetTitle("Abs. variation");
    fHistPairsRelDiff = FillHisto(fnPairsRelDiff, "RelDiffPair");
    fHistPairsRelDiff->GetYaxis()->SetTitle("Rel. variation (%)");
  }
}

void DreamSystematics::EvalDifferenceInParticles() {
  if (fNPartOneDefault.size() == 0 || fNPartOneVariations.size() == 0) {
    Error(
        "DreamSystematics",
        "DreamSystematics::EvalDifferenceInParticles() : default or var not set for part one");
  } else {
    EvalDifference(fNPartOneDefault, fNPartOneVariations, fnPartOneAbsDiff,
                   fnPartOneRelDiff);
    fHistPartOneAbsDiff = FillHisto(fnPartOneAbsDiff, "AbsDiffPartOne");
    fHistPartOneAbsDiff->GetYaxis()->SetTitle("Abs. variation");
    fHistPartOneRelDiff = FillHisto(fnPartOneRelDiff, "RelDiffPartOne");
    fHistPartOneRelDiff->GetYaxis()->SetTitle("Rel. variation (%)");
  }
  if (fNPartTwoDefault.size() == 0 || fNPartTwoVariations.size() == 0) {
    Error(
        "DreamSystematics",
        "DreamSystematics::EvalDifferenceInParticles() : default or var not set for part two");
  } else {
    EvalDifference(fNPartTwoDefault, fNPartTwoVariations, fnPartTwoAbsDiff,
                   fnPartTwoRelDiff);
    fHistPartTwoAbsDiff = FillHisto(fnPartTwoAbsDiff, "AbsDiffPartTwo");
    fHistPartTwoAbsDiff->GetYaxis()->SetTitle("Abs. variation");
    fHistPartTwoRelDiff = FillHisto(fnPartTwoRelDiff, "RelDiffPartTwo");
    fHistPartTwoRelDiff->GetYaxis()->SetTitle("Rel. variation (%)");
  }
}

void DreamSystematics::EvalDifferenceInPurity() {
  if (fPurityOneDefault.size() == 0 || fPurityOneVariations.size() == 0) {
    Error(
        "DreamSystematics",
        "DreamSystematics::EvalDifferenceInPurity() : default or var not set for part one");
  } else {
    EvalDifference(fPurityOneDefault, fPurityOneVariations, fPurityOneAbsDiff,
                   fPurityOneRelDiff);
    fHistPurityOneAbsDiff = FillHisto(fPurityOneAbsDiff, "AbsDiffPartOne");
    fHistPurityOneAbsDiff->GetYaxis()->SetTitle("Abs. variation");
    fHistPurityOneRelDiff = FillHisto(fPurityOneRelDiff, "RelDiffPartOne");
    fHistPurityOneRelDiff->GetYaxis()->SetTitle("Rel. variation (%)");
  }
  if (fPurityTwoDefault.size() == 0 || fPurityTwoVariations.size() == 0) {
    Error(
        "DreamSystematics",
        "DreamSystematics::EvalDifferenceInPurity() : default or var not set for part two");
  } else {
    EvalDifference(fPurityTwoDefault, fPurityTwoVariations, fPurityTwoAbsDiff,
                   fPurityTwoRelDiff);
    fHistPurityTwoAbsDiff = FillHisto(fPurityTwoAbsDiff, "AbsDiffPartTwo");
    fHistPurityOneAbsDiff->GetYaxis()->SetTitle("Abs. variation");
    fHistPurityTwoRelDiff = FillHisto(fPurityTwoRelDiff, "RelDiffPartTwo");
    fHistPurityTwoRelDiff->GetYaxis()->SetTitle("Rel. variation (%)");
  }
}

void DreamSystematics::ComputeUncertainty() {

  fHistSystErrAbs = (TH1F*) fHistDefault->Clone(
      Form("%s_SystErrAbs", fHistDefault->GetName()));
  fHistSystErrAbs->GetYaxis()->SetTitle("Syst. error");
  fHistSystErrAbs->Reset("ICMS");

  fHistSystErrRel = (TH1F*) fHistDefault->Clone(
      Form("%s_SystErrRel", fHistDefault->GetName()));
  fHistSystErrRel->GetYaxis()->SetTitle("Rel. syst. error");
  fHistSystErrRel->Reset("ICMS");
  const int nBins = fHistDefault->GetXaxis()->FindBin(fSystematicFitRangeUp);
  for (int ikstar = 1; ikstar <= nBins; ++ikstar) {
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
      case Pair::pL:
        EvalProtonLambda(ikstar);
        break;
      default:
        Error("DreamSystematics", "Non implemented Particle mode");
        break;
    }
  }

  fRatio = new TF1("SystError", "pol0(0)+expo(1)", fSystematicFitRangeLow,
                   fSystematicFitRangeUp);

  fRatio->SetParameter(
      0,
      fHistSystErrRel->GetBinContent(
          fHistSystErrRel->FindBin(fSystematicFitRangeUp)));
  float startExp = fHistSystErrRel->GetBinContent(1)
      - fHistSystErrRel->GetBinContent(
          fHistSystErrRel->FindBin(fSystematicFitRangeUp));

  fRatio->SetParameter(1, TMath::Log(startExp));
  std::cout
      << "Start Parameter 0: "
      << fHistSystErrRel->GetBinContent(
          fHistSystErrRel->FindBin(fSystematicFitRangeUp))
      << std::endl;
  std::cout << "difference: " << startExp << std::endl;
  std::cout << "Start Parameter 1: " << TMath::Log(startExp) << std::endl;
  fRatio->SetLineColor(kGreen + 2);
  fHistSystErrRel->Fit(fRatio, "WWRQ");

  fGrFinalError = new TGraphErrors();

  const float xerrLL = fHistDefault->GetBinWidth(1) / 2.;
  for (int kstarBin = 0; kstarBin < nBins; ++kstarBin) {
    const float x = fHistDefault->GetBinCenter(kstarBin + 1);
    const float y = fHistDefault->GetBinContent(kstarBin + 1);
    fGrFinalError->SetPoint(kstarBin, x, y);
    fGrFinalError->SetPointError(kstarBin, xerrLL, y * fRatio->Eval(x));
  }
  fGrFinalError->SetFillColor(kGray + 1);
  fGrFinalError->SetLineColor(kGray + 1);
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

void DreamSystematics::EvalProtonLambda(const int kstar) {
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
  // v0 pT
  addSyst.push_back(
      (fHistAbsErr[9]->GetBinContent(kstar)
          + fHistAbsErr[10]->GetBinContent(kstar)) / 2.);
  // CPA up
  addSyst.push_back(fHistAbsErr[11]->GetBinContent(kstar));
  // nSigma Up
  addSyst.push_back(fHistAbsErr[12]->GetBinContent(kstar));
  // Track TPC Cluster Up
  addSyst.push_back(fHistAbsErr[13]->GetBinContent(kstar));
  // Eta
  addSyst.push_back(
      (fHistAbsErr[14]->GetBinContent(kstar)
          + fHistAbsErr[15]->GetBinContent(kstar)) / 2.);
  // dca daughter
  addSyst.push_back(fHistAbsErr[16]->GetBinContent(kstar));
  // dca daughters to pv
  addSyst.push_back(fHistAbsErr[17]->GetBinContent(kstar));

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
  WriteOutput(file, fHistVar, "Raw");
  WriteOutput(file, fHistAbsErr, "AbsErr");
  WriteOutput(file, fHistErrBudget, "ErrBudget");
  WriteOutput(file, fHistBarlow, "Barlow");
  file->cd();
  fHistDefault->Write("histDefault");
  fHistSystErrAbs->Write("SystErrAbs");
  fHistSystErrRel->Write("SystErrRel");
  fRatio->Write("SystError");
  fGrFinalError->Write("grError");

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
  if (fPurityOneDefault.size() > 0 || fPurityTwoDefault.size() > 0) {
    file->mkdir("Purity");
    file->cd("Purity");
    if (fHistPurityOneAbsDiff)
      fHistPurityOneAbsDiff->Write();
    if (fHistPurityOneRelDiff)
      fHistPurityOneRelDiff->Write();
    if (fHistPurityTwoAbsDiff)
      fHistPurityTwoAbsDiff->Write();
    if (fHistPurityTwoRelDiff)
      fHistPurityTwoRelDiff->Write();
  }
  file->Close();
}

void DreamSystematics::WriteOutput(TFile* file, std::vector<TH1F*>& histvec,
                                   const TString name) {
  file->cd();
  file->mkdir(name.Data());
  file->cd(name.Data());

  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextSize(1.3 * gStyle->GetTextSize());

  auto canvasName = TString::Format("%s_%s", name.Data(), GetPairName().Data());
  auto c = new TCanvas(canvasName.Data(), canvasName.Data(), 1400, 1000);
  switch (fParticlePairMode) {
    case Pair::pp:
      c->Divide(3, 4);
      break;
    case Pair::pSigma0:
      c->Divide(7, 5);
      break;
    case Pair::pXi:
      c->Divide(7, 5);
      break;
    case Pair::pL:
      c->Divide(3, 6);
      break;
  }

  int iVar = 0;
  for (auto &it : histvec) {
    it->Write(Form("histVar_%i", iVar++));
    c->cd(iVar);
    it->Draw("pe");
    if (name == TString("ErrBudget")) {
      auto fit = it->GetFunction(Form("%s_fit", it->GetName()));
      if (fit) {
        text.DrawLatex(
            0.35,
            0.8,
            TString::Format("err. budget = %.2f %%",
                            fit->GetParameter(0) * 100.f));
        text.DrawLatex(
            0.35,
            0.73,
            TString::Format("low #it{k}* = %.2f %%",
                            it->GetBinContent(1) * 100.f));
      }
    } else if (name == TString("Barlow")) {
      it->GetYaxis()->SetRangeUser(0,it->GetMaximum()*2.0);
      c->cd(iVar);
      text.DrawLatex(
          gPad->GetUxmax()-0.8,gPad->GetUymax()-0.3,fBarlowLabel.at(iVar-1));
    }
  }
  c->Write();
  c->Print(TString::Format("CF_%s_%s.pdf", name.Data(), GetPairName().Data()));
  delete c;

  // For the raw CF write also the Default and all variations in one plot
  if (name == TString("Raw")) {
    auto c1 = new TCanvas("CF_var", "CF_var");
    fHistDefault->Draw();
    int iCount = 0;
    for (auto &it : fHistVar) {
      DreamPlot::SetStyleHisto(it, 20 + iCount, ++iCount);
      it->Draw("same");
    }
    c1->Write();
    c1->Print(TString::Format("CF_var_%s.pdf", GetPairName().Data()));
    delete c1;
  }

}
