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
  void SetPairs(DreamPair* pairOne, DreamPair* pairTwo) {
    fPairOne = pairOne;
    fPairTwo = pairTwo;
  }
  ;
  void GetCorrelations(const char* pairName = "");
  void LoopCorrelations(std::vector<DreamDist*> PairOne,
                        std::vector<DreamDist*> PairTwo, const char* name);
  void LoopCorrelations(std::vector<DreamDist*> Pair,const char* name);
  void WriteOutput(const char* name);
  TH1F* AddCF(TH1F* CF1, TH1F* CF2, const char* name);
  TH1F* ConvertToOtherUnit(TH1F* HistCF, int Scale, const char* name);
  std::vector<TH1F*> GetCorrelationFunctions() {
    return fCF;
  }
  ;
 private:
  std::vector<TH1F*> fCF;
  DreamPair* fPairOne;
  DreamPair* fPairTwo;
};

#endif /* DREAMCF_H_ */
