/*
 * ForgivingFitter.h
 *
 *  Created on: 13 Feb 2019
 *      Author: bernhardhohlweger
 */

#ifndef FORGIVINGQA_FORGIVINGFITTER_H_
#define FORGIVINGQA_FORGIVINGFITTER_H_
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"

class ForgivingFitter {
 public:
  ForgivingFitter();
  virtual ~ForgivingFitter();
  void FitInvariantMass(TH1F* histo, float massCutMin, float massCutMax, int signalColor, int backgroundColor);
  void SetRanges(float SigMin, float SigMax, float BkgRangeMin,
                 float BkgRangeMax);
  void SetRangesSigma(float SigMin, float SigMax, float BkgRangeMin,
                      float BkgRangeMax);
  int GetSignalCounts() const {
    return fSignalCounts;
  }
  ;
  int GetSignalCountsErr() const {
    return fSignalCountsErr;
  }
  int GetBackgroundCounts() const {
    return fBackgroundCounts;
  }
  ;
  int GetBackgroundCountsErr() const {
    return fBackgroundCountsErr;
  }
  double GetPurity() const;
  double GetPurityErr() const;
  float GetMeanMass() const {
    return fMeanMass;
  }
  float GetMeanWidth() const {
    return fMeanWidth;
  }
  void ShittyInvariantMass(TH1F* histo, TPad* c1, float pTMin, float pTMax,
                           const char* part);
  void FitInvariantMassSigma(TH1F* histo, float massCuts, int signalColor, int backgroundColor);
  TF1* GetBackgroundFunction() const { return fContinousBackGround; }
  TF1* GetSingleGaussian() const { return fSingleGaussian; }
  TF1* GetDoubleGaussian() const { return fDoubleGaussian; }
  TF1* GetFullFitFunction() const { return fFullFitFnct; }
 private:
  void CreateBackgroundFunction();
  void CreateContinousBackgroundFunction();
  void CreateSignalFunctions();
  void CreateFullFitFunction(TH1F* targetHisto);
  void SetStartParsDoubleGaussian(TH1F* targetHisto);
  float weightedMeanError(float weightA, float A, float weightB, float B,
                          float weightAErr, float AErr, float weightBErr,
                          float BErr);
  float weightedMean(float weightA, float A, float weightB, float B);
  TH1F *getSignalHisto(TF1 *function, TH1F *histo, float rangeLow,
                       float rangeHigh, const char *name);
  void CalculateBackgorund(TH1F* targetHisto, float massCutMin,
                           float massCutMax, int backgroundColor);
  TF1* fBackGround;
  TF1* fContinousBackGround;
  TF1* fSingleGaussian;
  TF1* fDoubleGaussian;
  TF1* fFullFitFnct;
  bool fRangesSet;
  float fBkgRangeMin;
  float fBkgRangeMax;
  float fSigRangeMin;
  float fSigRangeMax;
  int fSignalCounts;
  int fSignalCountsErr;
  int fWeightA;
  int fWeightB;
  int fBackgroundCounts;
  int fBackgroundCountsErr;
  double fMeanMass;
  double fMeanWidth;
};

inline
double ForgivingFitter::GetPurity() const {
  const double signal = static_cast<double>(fSignalCounts);
  const double bck = static_cast<double>(fBackgroundCounts);
  if (bck < 1E-6) {
    return 0;
  } else {
    return signal / (signal + bck);
  }
}

inline
double ForgivingFitter::GetPurityErr() const {
  const double signal = static_cast<double>(fSignalCounts);
  const double bck = static_cast<double>(fBackgroundCounts);
  const double signalErr = static_cast<double>(fSignalCountsErr);
  const double bckErr = static_cast<double>(fBackgroundCountsErr);
  const double signalBck = signal + bck;
  const double signalBckForth = signalBck * signalBck * signalBck * signalBck;
  if (bck < 1E-6) {
    return 0;
  } else {
    return std::sqrt(
        signalErr * signalErr * bck * bck / signalBckForth
            + bckErr * bckErr * signal * signal / signalBckForth);
  }
}
#endif /* FORGIVINGQA_FORGIVINGFITTER_H_ */
