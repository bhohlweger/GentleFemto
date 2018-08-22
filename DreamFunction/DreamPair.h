/*
 * DreamCF.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMPAIR_H_
#define DREAMFUNCTION_DREAMPAIR_H_
#include <vector>

#include "DreamDist.h"
class DreamPair {
 public:
  DreamPair(const char* name, float normleft, float normright);
  virtual ~DreamPair();
  void SetPair(DreamDist* Pair) {
    fPair = Pair;
    fPair->Calculate_CF(fNormLeft,fNormRight);
  }
  ;
  DreamDist* GetPair() {
    return fPair;
  }
  ;
  DreamDist* GetPairShiftedEmpty(unsigned int iIter) {
    return iIter < fPairShifted.size() ? fPairShifted.at(iIter) : nullptr;
  }
  ;
  DreamDist* GetPairFixShifted(unsigned int iIter) {
    return iIter < fPairFixShifted.size() ? fPairFixShifted.at(iIter) : nullptr;
  }
  ;
  DreamDist* GetPairRebinned(unsigned int iIter) {
    return iIter < fPairRebinned.size() ? fPairRebinned.at(iIter) : nullptr;
  }
  ;
  DreamDist* GetPairReweighted(unsigned int iIter) {
    return iIter < fPairReweighted.size() ? fPairReweighted.at(iIter) : nullptr;
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
  void FixShift(DreamDist* pair, DreamDist* otherDist, int kMin);
  void Rebin(DreamDist* pair, int rebin);
  void ReweightMixedEvent(DreamDist* pair, float kSMin, float kSMax);
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
