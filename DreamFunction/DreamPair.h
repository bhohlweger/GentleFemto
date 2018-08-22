/*
 * DreamPair.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMPAIR_H_
#define DREAMFUNCTION_DREAMPAIR_H_
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

class DreamPair {
 public:
  DreamPair(const char* name);
  DreamPair(DreamPair* pair,const char* name);
  virtual ~DreamPair();
  void SetSEDist(TH1F* SE) {fSE=(TH1F*)SE->Clone(Form("%s%s",SE->GetName(),fName));}
  TH1F* GetSEDist() {return fSE;};
  void SetSEMultDist(TH2F* SEMult) {fSEMult=(TH2F*)SEMult->Clone(Form("%s%s",SEMult->GetName(),fName));}
  TH2F* GetSEMultDist() {return fSEMult;};
  void SetMEDist(TH1F* ME) {fME=(TH1F*)ME->Clone(Form("%s%s",ME->GetName(),fName));}
  TH1F* GetMEDist() {return fME;};
  void SetMEMultDist(TH2F* MEMult) {fMEMult=(TH2F*)MEMult->Clone(Form("%s%s",MEMult->GetName(),fName));}
  TH2F* GetMEMultDist() {return fMEMult;};
  const char*   fName;
  TH1F*         fSE;
  TH2F*         fSEMult;
  TH1F*         fME;
  TH2F*         fMEMult;
};

#endif /* DREAMFUNCTION_DREAMPAIR_H_ */
