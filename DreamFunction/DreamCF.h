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
  void GetCorrelations();
  void LoopCorrelations(std::vector<DreamDist*> partPair,
                        std::vector<DreamDist*> antipartPair, const char* name);
  void WriteOutput(const char* name);
  TH1F* AddCF(TH1F* CF1, TH1F* CF2, const char* name);
 private:
  std::vector<TH1F*> fCF;
  DreamPair* fPartPair;
  DreamPair* fAntiPartPair;
};

#endif /* DREAMCF_H_ */
