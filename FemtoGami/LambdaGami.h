/*
 * LambdaGami.h
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#ifndef FEMTOGAMI_LAMBDAGAMI_H_
#define FEMTOGAMI_LAMBDAGAMI_H_
#include <vector>
#include "TH1F.h"

class LambdaGami {
 public:
  LambdaGami();
  virtual ~LambdaGami();
  TH1F* UnfoldResidual(TH1F* cf, TH1F* res, double lamRes);
  TH1F* UnfoldGenuine(TH1F* cf, double lamGen);
  void SetLambdaPar(double LamPars) {
    fLamPars.push_back(LamPars);
  }
  double GetLamdaPar(size_t idx) {
    return idx < fLamPars.size()?fLamPars.at(idx):999;
  }
 private:
  std::vector<double> fLamPars;
};

#endif /* FEMTOGAMI_LAMBDAGAMI_H_ */
