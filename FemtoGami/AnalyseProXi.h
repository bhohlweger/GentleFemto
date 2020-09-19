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
#include "DreamDist.h"
#include "DreamCF.h"

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
  DreamCF* ObtainCorrFunction(const char* name, DreamDist* partDist,
                              DreamDist* APartDist);
  void LimitCFRange(double val) {
    fLimitCFRange = true;
    fCFLimit = val;
  }
  void ResetLimits(DreamDist* dist);
  TH1F* LimitRange(TH1F* hist, double limit, const char* name);
  TH2F* LimitRange(TH2F* hist, double limit, const char* name);
  void SetAnalysisFile(const char* Path, const char* Prefix,
                       const char* Suffix = "0") {
    fFilename = Path;
    fPrefix = Prefix;
    fSuffix = Suffix;
  }
  void Default() {
    fNormVar = 2;
    fSideNormVar = 1;
    fBaselineVar = 0;
    fLamVarProton = 1;
    fLamVarOmega = 1;
    fLamVarXim1530 = 1;
    fRadVarXim1530 = 1;
    fMomGami->SetResponseVariation(1);
    fMomGami->SetIterVariation(5);
    fMomGami->SetUnfoldingMethod(MomentumGami::kBayes);
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
  void SetMomentumResolutionVar(int iVar) {
    fMomGami->SetResponseVariation(iVar);
  }
  void SetMomentumResolutionIter(int iIter) {
    fMomGami->SetIterVariation(iIter);
  }
  void SetMomentumUnfoldMethod(MomentumGami::Unfolding meth) {
    fMomGami->SetUnfoldingMethod(meth);
  }
 private:
  double fcutOff;  // at this value the calculation and doing of the cf stops
  bool fLimitCFRange;
  double fCFLimit;
  const char* fFilename;
  const char* fPrefix;
  const char* fSuffix;
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
  float fnormalizationLeft;
  float fnormalizationRight;
};

#endif /* FEMTOGAMI_ANALYSEPROXI_H_ */
