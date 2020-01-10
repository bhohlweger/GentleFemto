/*
 * DecayQA.h
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_DECAYQA_H_
#define FORGIVINGQA_DECAYQA_H_
#include "TList.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "ForgivingReader.h"
#include "ForgivingFitter.h"
#include "MakeHistosGreat.h"

class DecayQA {
 public:
  DecayQA(const char* partLatex, const char* latexProducts);
  DecayQA(const char* partLatex, const char* latexProducts, const char* outname);
  virtual ~DecayQA();
  void SetDecayCuts(TList* DecayCuts) {
    fDecayCuts = DecayCuts;
  }
  ;
  void SetAntiDecayCuts(TList* DecayCuts) {
    fAntiDecayCuts = DecayCuts;
  }
  ;
  void InvariantMassLambda(float CutMin, float CutMax, bool minBook = false,
                           float KaonCutMin = 0.48, float KaonCutMax = 0.515);
  void InvariantMassLambda(float CutMin, float CutMax, bool minBook = false);
  void InvariantMassPartLambda(float CutMin, float CutMax, bool minBook = false,
                           float KaonCutMin = 0.48, float KaonCutMax = 0.515);
  void InvariantMassAntiPartLambda(float CutMin, float CutMax, bool minBook = false,
                           float KaonCutMin = 0.48, float KaonCutMax = 0.515);
  void InvariantMassSigma0(float massCuts, const char* name = "Sigma0",
                           bool isSum = true);
  void SetStyler(DrawStyle styler) { fStyler = styler; }
  void GetPeriodQA(float CutMin, float CutMax,
                   std::vector<const char*> pathToList, const char* histname);
  void GetPeriodQASigma(float CutMin, float CutMax, const char* period);
  void GetPeriodQASigma0(float massCuts, const char* period);
  void InvariantMassXi(float CutMin, float CutMax);
  void InvariantMassXiMinBooking(float CutMin, float CutMax);
  void IvariantMassXiLambda();
  void PlotKaonRejection(TH1F* invMassKaon, const char* outname, float KaonCutMin, float KaonCutMax);
  void SetCanvasDivisions(unsigned int divX, unsigned int divY) {
    fDivCanX = divX;
    fDivCanY = divY;
  }
  ;
  void SetRangesFitting(float signalMin, float signalMax, float bkgMin,
                        float bkgMax);
  void SetInvMasspTStartBin(unsigned int start) {
    fInvMassPtStartBin = start;
  }
  ;
  void PlotQATopologyLambda();
  void PlotQATopologyPartLambda();
  void PlotQATopologyAntiPartLambda();
  void PlotQATopologyLambda(TList *v0Cuts, const char* outname);
  void PlotQATopologySigma0Daughter(TList *v0Cuts, const char* outname);
  void PlotQATopologySigma0(TList *v0Cuts, const char* outname);
  void PlotPIDSigma0Daughter(TList *v0Cuts, const char* outname);
  void PlotPIDLambda();
  void PlotPIDPartLambda();
  void PlotPIDAntiPartLambda();

  void PlotPIDLambda(TList *v0Cuts, const char* outname);
  void SetIMHistoScale(float scaleMaximum, float TexOffX, float TexOffY) {
    fScaleMax = scaleMaximum;
    fTexOffX = TexOffX;
    fTexOffY = TexOffY;
  }
  ;
  double GetSignalCounts() const {
    return fFitter->GetSignalCounts();
  }
  double GetSignalCountsErr() const {
    return fFitter->GetSignalCountsErr();
  }
  double GetBackgroundCounts() const {
    return fFitter->GetBackgroundCounts();
  }
  double GetBackgroundCountsErr() const {
    return fFitter->GetBackgroundCountsErr();
  }
  double GetPurity() const {
    return fFitter->GetPurity();
  }
  double GetPurityErr() const {
    return fFitter->GetPurityErr();
  }
  TGraphErrors* GetPurityGraph() const {
    return fPurity;
  }
  double GetPurity(double pt) const {
    return fPurity->Eval(pt);
  }
  float GetIntegratedPurity(size_t i) {
    return i<fIntegratedPurities.size()?fIntegratedPurities.at(i):-1;
  }
 private:
  void FitInvariantMass(TH2F* invMasspT, float CutMin, float CutMax,
                        const char* outname);
  void FitInvariantMassSigma0(TH2F* invMasspT, float massCuts,
                              const char* outname);
  ForgivingReader* fReader;
  MakeHistosGreat* fHairyPlotter;
  ForgivingFitter* fFitter;
  TList *fDecayCuts;
  TList *fAntiDecayCuts;
  unsigned int fDivCanX;
  unsigned int fDivCanY;
  unsigned int fInvMassPtStartBin;
  float fScaleMax;
  float fTexOffX;
  float fTexOffY;
  const char* fPartLatex;
  const char* fDecChannel;
  TGraphErrors *fPurity;
  std::vector<float> fIntegratedPurities;
  DrawStyle fStyler;
};

#endif /* FORGIVINGQA_DECAYQA_H_ */
