/*
 * DreamHEP.h
 *
 *  Created on: Mar 12, 2019
 *      Author: schmollweger
 */

#ifndef DREAMFUNCTION_DREAMHEP_H_
#define DREAMFUNCTION_DREAMHEP_H_
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

class DreamHEP {
 public:
  DreamHEP();
  virtual ~DreamHEP();
  void SetRootS(float energy) { fRootS = energy; }
  void SetMaxkStar(float maxk) { fkStarMax = maxk; }

  void printTH1HEPdata(const TH1* hist, const TGraphErrors* syst, const char* outname);
  void printTGAsymmHEPdata(const TGraphAsymmErrors* gr, const TGraphErrors* syst, const char* outname);
  TGraphErrors* GetSystErrHist(TH1* hist, TF1* syst);
  TGraphErrors* GetSystErrHist(TGraphAsymmErrors* gr, TF1* syst);

 private:
  float fRootS;
  float fkStarMax;
};

#endif /* DREAMFUNCTION_DREAMHEP_H_ */
