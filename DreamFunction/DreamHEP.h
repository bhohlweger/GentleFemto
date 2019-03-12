/*
 * DreamHEP.h
 *
 *  Created on: Mar 12, 2019
 *      Author: schmollweger
 */

#ifndef DREAMFUNCTION_DREAMHEP_H_
#define DREAMFUNCTION_DREAMHEP_H_
#include "TH1.h"
#include "TGraphErrors.h"

class DreamHEP {
 public:
  DreamHEP();
  virtual ~DreamHEP();
  void printTH1HEPdata(const TH1* hist, const TGraphErrors* syst, const char* outname);
};

#endif /* DREAMFUNCTION_DREAMHEP_H_ */
