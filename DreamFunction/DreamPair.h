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
  DreamPair();
  void SetPair(DreamDist* Pair) {fPair=Pair;};

  DreamDist* GetPair() {return fPair;};
  DreamDist* GetPairShifted(unsigned int iIter) {
    return iIter<fPairShifted.size()?fPairShifted.at(iIter):nullptr;};
  DreamDist* GetPairRebinned(unsigned int iIter) {
    return iIter<fPairRebinned.size()?fPairRebinned.at(iIter):nullptr;};
  DreamDist* GetPairReweighted(unsigned int iIter) {
    return iIter<fPairReweighted.size()?fPairReweighted.at(iIter):nullptr;};
  virtual ~DreamPair();
  void ShiftForEmpty(DreamDist* pair);
  void Rebin(DreamDist* pair, int rebin);
  void ReweightMixedEvent(DreamDist* pair,float kSMin,float kSMax);
  DreamDist*                fPair;
  std::vector<DreamDist*>   fPairShifted;
  std::vector<DreamDist*>   fPairRebinned;
  std::vector<DreamDist*>   fPairReweighted;
  int                       fFirstBin;
};

#endif /* DREAMFUNCTION_DREAMPAIR_H_ */
