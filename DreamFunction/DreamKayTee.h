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
    fKayTeeBins = x;
    fNKayTeeBins = (int) x.size();
  }
  ;
  void SetSEkTDist(int iPart, TH2F* SEkT) {
    fSEkT[iPart] = SEkT;
  }
  ;
  void SetMEkTDist(int iPart, TH2F* MEkT) {
    fMEkT[iPart] = MEkT;
  }
  ;
  void SetSEmTDist(int iPart, TH2F* SEmT) {
    fSEkT[iPart] = SEmT;
    fIskT = false;
  }
  void SetMEmTDist(int iPart, TH2F* MEmT) {
    fMEkT[iPart] = MEmT;
    fIskT = false;
  }
  void SetNormalization(float left, float right) {
    fNormleft = left;
    fNormright = right;
  }
  ;
  void ObtainTheCorrelationFunction(const char* outFolder, const char* prefix =
                                        "MB",
                                    const char* pair = "pp");
  void AveragekT(const char *pair);
  void SetSEMEReweightingRatio(const char* pathToFile, TString pair = "pp");
 private:
  bool fIskT;
  std::vector<float> fKayTeeBins;
  int fNKayTeeBins;
  TH2F* fSEkT[2];
  TH2F* fMEkT[2];
  TGraphErrors* fAveragekT;
  DreamPair*** fCFPart;
  DreamCF** fSum;
  float fNormleft;
  float fNormright;
  TH1F* fSEMEReweighting;
  TH1F* fSEMEReweightingMeV;
};

#endif /* DREAMFUNCTION_DREAMKAYTEE_H_ */
