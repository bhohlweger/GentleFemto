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
#include "TGraphErrors.h"
#include "LambdaGami.h"
#include "MomentumGami.h"
#include "TFile.h"

class AnalyseProXi {
 public:
  AnalyseProXi(double cutoff, double smearing);
  virtual ~AnalyseProXi();
  double SetupLambdaPars(LambdaGami* XiGami, double ProVar, double OmegaVar,
                         double Xi1530Var);
  TH1F* BaseLine(TH1F* dataCF);
  TH1F* XimSideband(LambdaGami* XiGami, TH1F* dataCF);
  TH1F* Xim1530FeedDown(LambdaGami* XiGami, TH1F* dataCF);
  TGraphErrors* GetCoulomb(TH1F* unfoldedGenuine);
  TGraphErrors* GetHalQCD(TH1F* unfoldedGenuine);
  TGraphErrors* GetESC16(TH1F* unfoldedGenuine);
  TH1F* GetVariation(int varnumber, bool getModels = false);
  void SetAnalysisFile(const char* Path, const char* Prefix) {
    fFilename = Path;
    fPrefix = Prefix;
  }
  void Default() {
    fNormVar = 1;
    fSideNormVar = 1;
    fBaselineVar = 1;
    fLamVarProton = 1;
    fLamVarOmega = 1;
    fLamVarXim1530 = 1;
    fRadVarXim1530 = 1;

  }
  ;
  void SetNormVar(int iVar) {
    fNormVar = iVar;
  }
  ;
  void SetSideNormVar(int iVar) {
    fSideNormVar = iVar;
  }
  ;
  void SetBaselineVar(int iVar) {
    fBaselineVar = iVar;
  }
  ;
  void SetLambdaVar(int iVarProton, int iVarOmega, int iVarXim) {
    fLamVarProton = iVarProton;
    fLamVarOmega = iVarOmega;
    fLamVarXim1530 = iVarXim;
  }
  ;
  void SetRadXim1530Var(int iVar) {
    fRadVarXim1530 = iVar;
  }
  ;
 private:
  double fcutOff;  // at this value the calculation and doing of the cf stops
  const char* fFilename;
  const char* fPrefix;
  TFile* fQAOutput;
  LambdaGami* fXiGami;
  MomentumGami* fMomGami;
  int fNormVar;
  std::vector<float> fNormMin;
  std::vector<float> fNormMax;
  int fSideNormVar;
  std::vector<float> fSidebandNormMin;
  std::vector<float> fSidebandNormMax;
  int fBaselineVar;
  std::vector<TString> fBLfunct;
  int fLamVarProton;
  int fLamVarOmega;
  int fLamVarXim1530;
  std::vector<float> fLamVars;
  int fRadVarXim1530;
  std::vector<float> fXim1530Rad;
};

#endif /* FEMTOGAMI_ANALYSEPROXI_H_ */
