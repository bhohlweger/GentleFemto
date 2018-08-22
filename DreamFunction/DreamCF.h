/*
 * DreamCF.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMCF_H_
#define DREAMFUNCTION_DREAMCF_H_
#include "DreamPair.h"

class DreamCF {
 public:
  DreamCF();
  void SetPair(DreamPair* Pair) {fPair=Pair;};
  virtual ~DreamCF();
  void ShiftForEmpty(DreamPair* pair);
  DreamPair*  fPair;
  DreamPair*  fPairShifted;
  DreamPair*  fPairLast;
};

#endif /* DREAMFUNCTION_DREAMCF_H_ */
