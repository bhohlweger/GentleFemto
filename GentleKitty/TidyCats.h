/*
 * TidyCats.h
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#ifndef GENTLEKITTY_TIDYCATS_H_
#define GENTLEKITTY_TIDYCATS_H_

class TidyCats {
 public:
  enum Sources {
    sGaussian,
    sResonance,
    sLevy
  };
  enum pXimPot {
    pHALQCD,
    pHaidenbauer,
    pRikken
  };
  TidyCats();
  virtual ~TidyCats();
  void GetCatsProtonProton(CATS* AB_pp, int momBins, double kMin, double kMax,
                           TidyCats::Sources source);
  void GetCatsProtonLambda(CATS* AB_pL, int momBins, double kMin, double kMax,
                           TidyCats::Sources source);
  void GetCatsProtonXiMinus(CATS* AB_pXim, int momBins, double kMin,
                            double kMax, TidyCats::Sources source,
                            TidyCats::pXimPot pot, double QCDTime);
  void GetCatsProtonXiMinusCutOff(CATS* AB_pXim, int momBins, double kMin,
                                  double kMax, bool StrongOn, double QCDTime,
                                  double cutOff);
  void GetCatsProtonXiMinus1530(CATS* AB_pXim1530, int momBins, double kMin,
                                double kMax);
 private:
};

#endif /* GENTLEKITTY_TIDYCATS_H_ */
