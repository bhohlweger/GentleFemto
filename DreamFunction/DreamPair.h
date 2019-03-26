/*
 * DreamCF.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMPAIR_H_
#define DREAMFUNCTION_DREAMPAIR_H_
#include <vector>
#include "TFile.h"
#include "DreamDist.h"
class DreamPair {
 public:
  DreamPair(const char* name, float normleft=0, float normright=0);
  virtual ~DreamPair();
  void SetPair(DreamDist* Pair) {
    fPair = Pair;
    fPair->Calculate_CF(fNormLeft, fNormRight);
  }
  ;
  DreamDist* GetPair() {
    return fPair;
  }
  ;
  unsigned int GetFemtoPairs(float kMin, float kMax) {
    return fPair?fPair->GetFemtoPairs(kMin,kMax):0;
  }
  DreamDist* GetPairShiftedEmpty(unsigned int iIter) {
    return iIter < fPairShifted.size() ? fPairShifted.at(iIter) : nullptr;
  }
  ;
  std::vector<DreamDist*> GetShiftedEmpty() {
    return fPairShifted;
  }
  ;
  DreamDist* GetPairFixShifted(unsigned int iIter) {
    return iIter < fPairFixShifted.size() ? fPairFixShifted.at(iIter) : nullptr;
  }
  ;
  std::vector<DreamDist*> GetFixShifted() {
    return fPairFixShifted;
  }
  ;
  DreamDist* GetPairRebinned(unsigned int iIter) {
    return iIter < fPairRebinned.size() ? fPairRebinned.at(iIter) : nullptr;
  }
  ;
  std::vector<DreamDist*> GetRebinned() {
    return fPairRebinned;
  }
  ;
  DreamDist* GetPairReweighted(unsigned int iIter) {
    return iIter < fPairReweighted.size() ? fPairReweighted.at(iIter) : nullptr;
  }
  ;
  std::vector<DreamDist*> GetReweighted() {
    return fPairReweighted;
  }
  ;
  const char* GetName() {
    return fName;
  }
  ;
  float GetFirstBin() {
    return fFirstBin;
  }
  ;
  int GetNDists();
  void ShiftForEmpty(DreamDist* pair);
  void FixShift(DreamDist* pair, DreamDist* otherDist, float kMin,
                const bool fixedShift = false);
  void FixShift(DreamDist* pair, DreamDist* otherPair1, DreamDist* otherPair2,
                float kMin1, float kMin2);
  void Rebin(DreamDist* pair, int rebin);
  void ReweightMixedEvent(DreamDist* pair, float kSMin, float kSMax);
  void WriteOutput(TList *Outlist);
 private:
  DreamDist* fPair;
  std::vector<DreamDist*> fPairShifted;
  std::vector<DreamDist*> fPairFixShifted;
  std::vector<DreamDist*> fPairRebinned;
  std::vector<DreamDist*> fPairReweighted;
  float fFirstBin;
  float fNormLeft;
  float fNormRight;
  const char* fName;
};

#endif /* DREAMFUNCTION_DREAMPAIR_H_ */
