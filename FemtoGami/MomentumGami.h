/*
 * MomentumGami.h
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#ifndef FEMTOGAMI_MOMENTUMGAMI_H_
#define FEMTOGAMI_MOMENTUMGAMI_H_
#include "TH2F.h"
#include "TH1F.h"

class MomentumGami {
 public:
  MomentumGami(float maxkStar);
  double Eval(double *x, double *p);
  virtual ~MomentumGami();
  void SetResolution(TH2F* resoMatrix, float UnitConversion = 1);
  void Unfold(TH1F* InputDist);
  TH1F* Fold(TH1F* InputDist);
  Double_t func(int ii, Double_t *par);
  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  TH1F* UnfoldviaTSVD(TH1F* InputDist, TList* QA);
 private:
  TH2F* fResolution;
  std::vector<TH1F*> fResProjection;
  TH1F* fToUnfold;
  float fMaxkStar;
};

#endif /* FEMTOGAMI_MOMENTUMGAMI_H_ */
