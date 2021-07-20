/*
 * MomentumGami.h
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#ifndef FEMTOGAMI_MOMENTUMGAMI_H_
#define FEMTOGAMI_MOMENTUMGAMI_H_
#include "TH2F.h"
#include "TH1F.h"

#include "RooUnfoldResponse.h"

class MomentumGami {
 public:
  enum Unfolding {
    kBayes = 0,
    kB2B = 1,
    kIDS = 2
  };
  MomentumGami(float maxkStar);
  double Eval(double *x, double *p);
  virtual ~MomentumGami();
  void SetResolution(TH2F* resoMatrix, float UnitConversion = 1);
  void UnfoldGuessing(TH1F* InputDist);
  TH1F* Fold(TH1F* InputDist);
  Double_t func(int ii, Double_t *par);
  void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  TH1F* UnfoldviaTSVD(TH1F* InputDist, TList* QA);
  TH1F* UnfoldviaRooResp(TH1F* InputDist, double Rescaling = 1);
  void TrainRooResponse(TH2F* resp, RooUnfoldResponse* roo_resp);
  void SetResponseVariation(int iVar) {
    fResp = iVar;
  }
  void SetIterVariation(int iVar) {
    fIter = iVar;
  }
  void SetUnfoldingMethod(Unfolding meth) {
    fHowToUnfold = meth;
  }
  void SetDoMomentumResolutionParametrization(bool doIt) {
    fDoMomResParametrization = doIt;
  }
  TList* GetQAList() {return fQAList; }
 private:
  bool fResponseMatrixSetup; 
  TList* fQAList;
  TH2F* fResolution;
  std::vector<TH1F*> fResProjection;
  TH1F* fToUnfold;
  TH2F* fTrainQA;
  int fResp;
  int fIter;
  RooUnfoldResponse *fResponseDefault;
  RooUnfoldResponse *fResponseLower;
  RooUnfoldResponse *fResponseUpper;
  float fMaxkStar;
  float fUnitConversion;
  Unfolding fHowToUnfold;
  bool fDoMomResParametrization;
};

#endif /* FEMTOGAMI_MOMENTUMGAMI_H_ */
