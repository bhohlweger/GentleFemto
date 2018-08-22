/*
 * DreamCF.h
 *
 *  Created on: Aug 22, 2018
 *      Author: hohlweger
 */

#ifndef DREAMCF_H_
#define DREAMCF_H_
#include "DreamPair.h"
class DreamCF {
 public:
  DreamCF();
  virtual ~DreamCF();
 private:
  DreamPair* fParticlePair;
  DreamPair* fAntiParticlePair;
};

#endif /* DREAMCF_H_ */
