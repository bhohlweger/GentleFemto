/*
 * AnalyseProXi.h
 *
 *  Created on: 29 Oct 2019
 *      Author: bernhardhohlweger
 */

#ifndef FEMTOGAMI_ANALYSEPROXI_H_
#define FEMTOGAMI_ANALYSEPROXI_H_

#include "TH1F.h"
#include "TSystem.h"
#include "LambdaGami.h"
#include "TFile.h"

class AnalyseProXi {
 public:
  AnalyseProXi(double cutoff);
  virtual ~AnalyseProXi();
  double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
                         double Xi1530Var);
  TH1F* BaseLine(TH1F* dataCF);
  TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF, unsigned int varSideNorm);
  TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF);
  void StoreModels(TH1F* unfoldedGenuine, TFile* QAOutput);
  TH1F* GetVariation(int varnumber);
  void SetAnalysisFile(const char* Path, const char* Prefix) {
    fFilename = Path;
    fPrefix = Prefix;
  }
 private:
  const char* fFilename;
  const char* fPrefix;
  double fcutOff;  // at this value the calculation and doing of the cf stops
  TFile* fQAOutput;

};

#endif /* FEMTOGAMI_ANALYSEPROXI_H_ */
