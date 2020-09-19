/*
 * DreamPair.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamDist.h"

#include <iostream>
DreamDist::DreamDist()
    : fSE(nullptr),
      fSEMult(nullptr),
      fME(nullptr),
      fMEMult(nullptr),
      fCF(nullptr),
      fGrCF(nullptr) {
}
DreamDist::DreamDist(DreamDist* pair, const char* name)
    : fSE(nullptr),
      fSEMult(nullptr),
      fME(nullptr),
      fMEMult(nullptr),
      fCF(nullptr),
      fGrCF(nullptr) {
  this->SetSEDist(pair->GetSEDist(), name);
  if (pair->GetSEMultDist())
    this->SetSEMultDist(pair->GetSEMultDist(), name);
  this->SetMEDist(pair->GetMEDist(), name);
  if (pair->GetMEMultDist())
    this->SetMEMultDist(pair->GetMEMultDist(), name);
}

DreamDist::~DreamDist() {
  if (fSE)
    delete fSE;
  if (fSEMult)
    delete fSEMult;
  if (fME)
    delete fME;
  if (fMEMult)
    delete fMEMult;
  if (fCF)
    delete fCF;
  if (fGrCF)
    delete fGrCF;
}

unsigned int DreamDist::GetFemtoPairs(float kMin, float kMax) {
  unsigned int iPairs = 0;
  if (fSE) {
    if (fSE->GetXaxis()->GetXmin() <= kMin
        && kMax <= fSE->GetXaxis()->GetXmax()) {
      iPairs = fSE->Integral(fSE->FindBin(kMin), fSE->FindBin(kMax));
    } else {
      std::cout
          << "DreamDist::GetFemtoPairs: Histogram and kMin/Max range do not add up: \n"
          << "Histogram range (xmin - xmax): " << fSE->GetXaxis()->GetXmin()
          << " - " << fSE->GetXaxis()->GetXmax() << '\n'
          << "Integration range (kMin - kMax): " << kMin << " - " << kMax
          << std::endl;
    }
  } else {
    std::cout << "DreamDist::GetFemtoPairs SE does not exist \n";
  }
  return iPairs;
}

void DreamDist::Calculate_CF(float normleft, float normright,
                             TH1F* hSEnorebin) {
  if (!fCF) {
    TString CFname = fSE->GetName();
    CFname.Replace(CFname.Index("SE"), 2, "CF");
    fCF = (TH1F*) fSE->Clone(CFname.Data());
    fGrCF = new TGraphAsymmErrors(fCF);
    fGrCF->Set(0);

    Double_t norm_relK = 0;
    double IntegralSE = fSE->Integral(fSE->FindBin(normleft),
                                      fSE->FindBin(normright));
    double IntegralME = fME->Integral(fME->FindBin(normleft),
                                      fME->FindBin(normright));
    if (IntegralME != 0) {
      norm_relK = IntegralSE / IntegralME;
      fCF->Divide(fSE, fME, 1, norm_relK);
      // Shift the center of the bin in x according to the same event distribution
      // In case of very large bins, this can have a sizeable effect!
      int counter = 0;
      float xVal, xErrRight, xErrLeft, centrVal;
      for (int i = 1; i <= fCF->GetNbinsX(); ++i) {
        // In case we have a SE distribution available we use the mean in that k* bin
        if (hSEnorebin) {
          centrVal = fCF->GetBinCenter(i);
          double chepsilyon = fCF->GetBinWidth(i)*1e-4;
          hSEnorebin->GetXaxis()->SetRangeUser(
              fCF->GetBinLowEdge(i) + chepsilyon,
              fCF->GetBinLowEdge(i) + fCF->GetBinWidth(i) - chepsilyon);
          xVal = hSEnorebin->GetMean();
          xErrLeft = xVal - centrVal + fCF->GetBinWidth(i) / 2.;
          xErrRight = centrVal - xVal + fCF->GetBinWidth(i) / 2.;
        } else {
          // else just the bin center
          xVal = fCF->GetBinCenter(i);
          xErrLeft = fCF->GetBinWidth(i) / 2.;
          xErrRight = fCF->GetBinWidth(i) / 2.;
        }
        fGrCF->SetPoint(counter, xVal, fCF->GetBinContent(i));
        fGrCF->SetPointError(counter++, xErrLeft, xErrRight,
                             fCF->GetBinError(i), fCF->GetBinError(i));
      }
    } else {
      std::cout << "DreamDist::Calculate_CF division by 0 \n";
      std::cout << "Normalization left: " << normleft << " and right: "
                << normright << std::endl;
    }

  } else {
    std::cout << fCF->GetName() << " was already set, skipping\n";
  }
  return;
}
