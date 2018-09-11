/*
 * DreamKayTee.h
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_DREAMKAYTEE_H_
#define DREAMFUNCTION_DREAMKAYTEE_H_
#include <vector>
#include <iterator>
#include "TH2F.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "DreamPair.h"
#include "DreamCF.h"
class DreamKayTee {
 public:
  DreamKayTee();
  virtual ~DreamKayTee();
  void SetKayTeeBins(std::vector<float> x) {
    fKayTeeBins=x;
    fNKayTeeBins=(int)x.size();
  };
  void SetSEkTDist(int iPart,TH2F* SEkT) {
    fSEkT[iPart]=SEkT;
  };
  void SetMEkTDist(int iPart,TH2F* MEkT) {
    fMEkT[iPart]=MEkT;
  };
  void SetNormalization(float left, float right) {
    fNormleft=left;
    fNormright=right;
  };
  void ObtainTheCorrelationFunction();
  void AveragekT();
 private:
  std::vector<float> fKayTeeBins;
  int fNKayTeeBins;
  TH2F* fSEkT[2];
  TH2F* fMEkT[2];
  TGraphErrors* fAveragekT;
  DreamPair*** fCFPart;
  DreamCF** fSum;
  float fNormleft;
  float fNormright;
};

#endif /* DREAMFUNCTION_DREAMKAYTEE_H_ */
