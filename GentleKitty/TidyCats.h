/*
 * TidyCats.h
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "TString.h"
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
    pCoulomb,
    pGamow,
    pHALQCD,
    pHALQCDGamow,
    pHaidenbauer,
    pRikkenWF,
    pRikkenPot
  };
  enum pLPot {
    pUsmani,
    pNLOLed,
    pNLOWF,
    pLOLed,
    pLOWF
  };
  enum pSigma0Pot {
    pSigma0Haidenbauer,
    pSigma0ESC16
  };
  TidyCats();
  virtual ~TidyCats();
  void GetCatsProtonProton(CATS* AB_pp, int momBins, double kMin, double kMax,
                           TidyCats::Sources source);
  DLM_CleverMcLevyReso* GetSourceProtonProton() {
    return fppCleverMcLevy;
  }
  void GetCatsProtonLambda(CATS* AB_pL, int momBins, double kMin, double kMax,
                           TidyCats::Sources source, TidyCats::pLPot pot);
  void GetCatsProtonXiMinus(CATS* AB_pXim, int momBins, double kMin,
                            double kMax, TidyCats::Sources source,
                            TidyCats::pXimPot pot, double QCDTime);
  void GetCatsProtonXiMinusCutOff(CATS* AB_pXim, int momBins, double kMin,
                                  double kMax, bool StrongOn, double QCDTime,
                                  double cutOff);
  void GetCatsProtonXiMinus1530(CATS* AB_pXim1530, int momBins, double kMin,
                                double kMax, TidyCats::Sources source);
  void GetCatsProtonSigma0(CATS* AB_pSigma0, int momBins, double kMin,
                           double kMax, TidyCats::Sources source,
                           TidyCats::pSigma0Pot pot);
  static double ESC16_pXim_EXAMPLE(double* Parameters);
 private:
  DLM_CleverLevy* fppCleverLevy;
  DLM_CleverMcLevyReso* fppCleverMcLevy;
  DLM_CleverLevy* fpLCleverLevy;
  DLM_CleverMcLevyReso* fpLCleverMcLevy;
  DLM_CleverLevy* fpXimCleverLevy;
  DLM_CleverMcLevyReso* fpXimCleverMcLevy;
  DLM_CleverLevy* fpXim1530CleverLevy;
  DLM_CleverMcLevyReso* fpSigma0CleverMcLevy;
};

#endif /* GENTLEKITTY_TIDYCATS_H_ */
