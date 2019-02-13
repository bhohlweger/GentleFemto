/*
 * DecayQA.h
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_DECAYQA_H_
#define FORGIVINGQA_DECAYQA_H_
#include "TList.h"
#include "ForgivingReader.h"
#include "MakeHistosGreat.h"

class DecayQA {
 public:
  DecayQA();
  virtual ~DecayQA();
  void SetDecayCuts(TList* trkCuts) {
    fDecayCuts = trkCuts;
  }
  ;
  void SetAntiDecayCuts(TList* trkCuts) {
    fAntiDecayCuts = trkCuts;
  }
  ;
 private:
  ForgivingReader* fReader;
  MakeHistosGreat* fHairyPlotter;
  TList *fDecayCuts;
  TList *fAntiDecayCuts;
};

#endif /* FORGIVINGQA_DECAYQA_H_ */
