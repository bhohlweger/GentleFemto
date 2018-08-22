/*
 * DreamCF.cxx
 *
 *  Created on: Aug 22, 2018
 *      Author: hohlweger
 */

#include "DreamCF.h"

DreamCF::DreamCF()
    : fSE(),
      fME(),
      fCF(),
      fPartPair(nullptr),
      fAntiPartPair(nullptr) {

}

DreamCF::~DreamCF() {
}

void DreamCF::GetCorrelations(const char* outname) {
  if (fPartPair->GetNDists() == fAntiPartPair->GetNDists()) {
    unsigned int iIterShifted = 0;
    while (fPartPair->GetPairShiftedEmpty(iIterShifted)
        && fAntiPartPair->GetPairShiftedEmpty(iIterShifted)) {

      iIterShifted++;
    }
  } else {
    std::cout << "Part Pair with " << fPartPair->GetNDists()
              << "Distributions, Anti Part Pair with "
              << fAntiPartPair->GetNDists() << std::endl;
  }
  return;
}

void DreamCF::WriteOutput() {
  return;
}
