/*
 * DreamPair.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMDIST_H_
#define DREAMFUNCTION_DREAMDIST_H_
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

class DreamDist {
 public:
  DreamDist();
  DreamDist(DreamDist* pair, const char* name);
  virtual ~DreamDist();
  void SetSEDist(TH1F* SE, const char* name) {
    fSE = (TH1F*) SE->Clone(Form("%s%s", SE->GetName(), name));
  }
  TH1F* GetSEDist() {
    return fSE;
  }
  ;
  void SetSEMultDist(TH2F* SEMult, const char* name) {
    fSEMult = (TH2F*) SEMult->Clone(Form("%s%s", SEMult->GetName(), name));
  }
  TH2F* GetSEMultDist() {
    return fSEMult;
  }
  ;
  void SetMEDist(TH1F* ME, const char* name) {
    fME = (TH1F*) ME->Clone(Form("%s%s", ME->GetName(), name));
  }
  TH1F* GetMEDist() {
    return fME;
  }
  ;
  void SetMEMultDist(TH2F* MEMult, const char* name) {
    fMEMult = (TH2F*) MEMult->Clone(Form("%s%s", MEMult->GetName(), name));
  }
  TH2F* GetMEMultDist() {
    return fMEMult;
  }
  ;
  TH1F* GetCF() {
    return fCF;
  }
  ;
  void Calculate_CF(float normleft,float normright);
 private:
  TH1F* fSE;
  TH2F* fSEMult;
  TH1F* fME;
  TH2F* fMEMult;
  TH1F* fCF;
};

#endif /* DREAMFUNCTION_DREAMDIST_H_ */
