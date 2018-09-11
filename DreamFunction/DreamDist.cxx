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
      fCF(nullptr) {
}
DreamDist::DreamDist(DreamDist* pair, const char* name)
    : fSE(nullptr),
      fSEMult(nullptr),
      fME(nullptr),
      fMEMult(nullptr),
      fCF(nullptr) {
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
}

void DreamDist::Calculate_CF(float normleft, float normright) {
  if (!fCF) {
    TString CFname = fSE->GetName();
    CFname.Replace(CFname.Index("SE"), 2, "CF");
    fCF = (TH1F*) fSE->Clone(CFname.Data());
    Double_t norm_relK = 0;
    double IntegralSE = fSE->Integral(fSE->FindBin(normleft),
                                      fSE->FindBin(normright));
    double IntegralME = fME->Integral(fME->FindBin(normleft),
                                      fME->FindBin(normright));
    norm_relK = IntegralSE / IntegralME;

    fCF->Divide(fSE, fME, 1, norm_relK);
  } else {
    std::cout << fCF->GetName() << " was already set, skipping\n";
  }
  return;
}
