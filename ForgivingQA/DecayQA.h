/*
 * DecayQA.h
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_DECAYQA_H_
#define FORGIVINGQA_DECAYQA_H_
#include "TList.h"
#include "TH2F.h"
#include "ForgivingReader.h"
#include "ForgivingFitter.h"
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
  void SetRangesFitting(float signalMin, float signalMax, float bkgMin,
                        float bkgMax);
  void InvariantMassLambda(float CutMin, float CutMax);
  void KaonRejectionWindow(TH1F* invMassKaon);
  void SetCanvasDivisions(unsigned int divX, unsigned int divY) {
    fDivCanX = divX;
    fDivCanY = divY;
  }
  ;
 private:
  void FitInvariantMass(TH2F* invMasspT, float CutMin, float CutMax,
                        const char* outname);
  ForgivingReader* fReader;
  MakeHistosGreat* fHairyPlotter;
  ForgivingFitter* fFitter;
  TList *fDecayCuts;
  TList *fAntiDecayCuts;
  unsigned int fDivCanX;
  unsigned int fDivCanY;
};

#endif /* FORGIVINGQA_DECAYQA_H_ */
