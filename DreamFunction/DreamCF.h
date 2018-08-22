/*
 * DreamCF.h
 *
 *  Created on: Aug 22, 2018
 *      Author: hohlweger
 */

#ifndef DREAMCF_H_
#define DREAMCF_H_
#include <vector>
#include "DreamPair.h"
#include "TH1F.h"

class DreamCF {
 public:
  DreamCF();
  virtual ~DreamCF();
  void SetPairs(DreamPair* pairPart, DreamPair* pairAntiPart) {
    fPartPair = pairPart;
    fAntiPartPair = pairAntiPart;
  }
  ;
  void GetCorrelations(const char* outName);
  void SetOutputName(const char *name) {
    fOutputNames.push_back(name);
  }
  ;
  void WriteOutput();
 private:
  std::vector<TH1F*> fSE;
  std::vector<TH1F*> fME;
  std::vector<TH1F*> fCF;
  DreamPair* fPartPair;
  DreamPair* fAntiPartPair;
  std::vector<const char*> fOutputNames;
};

#endif /* DREAMCF_H_ */
