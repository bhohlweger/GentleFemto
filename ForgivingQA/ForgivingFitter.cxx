/*
 * ForgivingFitter.cxx
 *
 *  Created on: 13 Feb 2019
 *      Author: bernhardhohlweger
 */

#include "ForgivingFitter.h"

ForgivingFitter::ForgivingFitter()
    : fBackGround(nullptr),
      fContinousBackGround(nullptr),
      fSingleGaussian(nullptr),
      fDoubleGaussian(nullptr),
      fFullFitFnct(nullptr),
      fRangesSet(false),
      fBkgRangeMin(0),
      fBkgRangeMax(0),
      fSigRangeMin(0),
      fSigRangeMax(0),
      fSignalCounts(0),
      fBackgroundCounts(0) {

}

ForgivingFitter::~ForgivingFitter() {
  // TODO Auto-generated destructor stub
}

TH1F *getSignalHisto(TF1 *function, TH1F *histo, float rangeLow,
                     float rangeHigh, const char *name) {
  const int firstBin = histo->FindBin(rangeLow);
  const int lastBin = histo->FindBin(rangeHigh);
  TH1F *result = new TH1F(
      Form("result_%.2f_%.2f_%s", rangeLow, rangeHigh, name), "",
      histo->GetNbinsX(), histo->GetXaxis()->GetXmin(),
      histo->GetXaxis()->GetXmax());
  for (int i = firstBin; i < lastBin; ++i) {
    float weight = histo->GetBinContent(i)
        - function->Eval(histo->GetBinCenter(i));
    result->Fill(histo->GetBinCenter(i), weight);
    result->SetBinError(i, histo->GetBinError(i));
  }
  result->SetFillColor(kGray + 1);
  result->SetLineColor(kGray + 1);
  return result;
}

// second order polynomial + double gaus to Lambda peak
void ForgivingFitter::FitInvariantMass(TH1F* histo, float massCutMin,
                                float massCutMax) {
  // Fit Background with second order polynomial, excluding Mlambda +/- 10 MeV
  if (!fRangesSet) {
    std::cout
        << "no BackGround Function defined via SetBackGroundRange! Exiting \n";
    return;
  }
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SRQ0", "",
                                         fBkgRangeMin * 1.01,
                                         fBkgRangeMax * 0.99);
  CreateContinousBackgroundFunction();
  // remove background from signal
  TH1F *signalOnly = getSignalHisto(fContinousBackGround, histo,
                                    fSigRangeMin * 0.95, fSigRangeMax * 1.05,
                                    Form("%s_signal_only", histo->GetName()));
  //  signalOnly->DrawCopy();
  signalOnly->Fit("fSignalSingleGauss");

  TFitResultPtr r = signalOnly->Fit("fSignalDoubleGauss", "SRQ0", "",
                                    fSigRangeMin, fSigRangeMax);
  // Extract signal as integral
  fSignalCounts = fDoubleGaussian->Integral(massCutMin, massCutMax)
      / double(histo->GetBinWidth(1));
//  fSignalErr = fDoubleGaussian->IntegralError(
//      lowerBound, upperBound, r->GetParams(),
//      r->GetCovarianceMatrix().GetMatrixArray())
//      / double(histo->GetBinWidth(1));
  CreateFullFitFunction(histo);
  histo->Fit("fLambda", "SRQ", "", fBkgRangeMin * 1.01, fBkgRangeMax * 0.99);
  CalculateBackgorund(histo, massCutMin, massCutMax);
}

void ForgivingFitter::SetRanges(float SigMin, float SigMax, float BkgRangeMin,
                                float BkgRangeMax) {
  fBkgRangeMin = BkgRangeMin;
  fBkgRangeMax = BkgRangeMax;
  fSigRangeMin = SigMin;
  fSigRangeMax = SigMax;
  CreateBackgroundFunction();
  CreateSignalFunctions();
}

void ForgivingFitter::CreateBackgroundFunction() {
  if (fBackGround) {
    delete fBackGround;
  }
  fBackGround = new TF1("fBackground", [&](double *x, double *p) {
    if (x[0] > fSigRangeMin && x[0] < fSigRangeMax) {TF1::RejectPoint(); return
      (double)0;}return p[0] + p[1]*x[0] + p[2]*x[0]*x[0];},
                        fBkgRangeMin, fBkgRangeMax, 3);
}

void ForgivingFitter::CreateContinousBackgroundFunction() {
  if (fContinousBackGround) {
    delete fContinousBackGround;
  }
  fContinousBackGround = new TF1("fBackground2", "pol2", fBkgRangeMin,
                                 fBkgRangeMax);
  fContinousBackGround->SetParameter(0, fBackGround->GetParameter(0));
  fContinousBackGround->SetParameter(1, fBackGround->GetParameter(1));
  fContinousBackGround->SetParameter(2, fBackGround->GetParameter(2));
}

void ForgivingFitter::CreateSignalFunctions() {
  fSingleGaussian = new TF1("fSignalSingleGauss", "gaus(0)", fSigRangeMin,
                            fSigRangeMax);
  fDoubleGaussian = new TF1("fSignalDoubleGauss", "gaus(0) + gaus(3)",
                            fSigRangeMin, fSigRangeMax);
}

void ForgivingFitter::SetStartParsDoubleGaussian(TH1F* targetHisto) {
  fDoubleGaussian->SetParameter(0, 0.75 * targetHisto->GetMaximum());
  fDoubleGaussian->SetParameter(1, fSingleGaussian->GetParameter(1));
  fDoubleGaussian->SetParameter(2, 2.f * fSingleGaussian->GetParameter(2));
  fDoubleGaussian->SetParLimits(2, 0.5 * fSingleGaussian->GetParameter(2),
                                1e2 * 2.f * fSingleGaussian->GetParameter(2));
  fDoubleGaussian->SetParameter(3, 0.2 * targetHisto->GetMaximum());
  fDoubleGaussian->SetParameter(4, fSingleGaussian->GetParameter(1));
  fDoubleGaussian->SetParLimits(
      4, fSingleGaussian->GetParameter(1) - fSingleGaussian->GetParameter(2),
      fSingleGaussian->GetParameter(1) + fSingleGaussian->GetParameter(2));
  fDoubleGaussian->SetParameter(5, 0.5 * fSingleGaussian->GetParameter(2));
  fDoubleGaussian->SetParLimits(5, 0.5 * fSingleGaussian->GetParameter(2),
                                1e2 * 2.f * fSingleGaussian->GetParameter(2));
}

void ForgivingFitter::CreateFullFitFunction(TH1F* targetHisto) {
  if (fFullFitFnct) {
    delete fFullFitFnct;
  }
  TF1 *fFullFitFnct = new TF1("fLambda", "fBackground2 + fSignalDoubleGauss",
                              fBkgRangeMin, fBkgRangeMax);
  fFullFitFnct->SetNpx(1000);
  fFullFitFnct->SetParameter(3, 0.75 * targetHisto->GetMaximum());
  fFullFitFnct->SetParameter(4, fDoubleGaussian->GetParameter((1)));
  fFullFitFnct->SetParameter(5, fDoubleGaussian->GetParameter((2)));
  fFullFitFnct->SetParameter(6, 0.2 * targetHisto->GetMaximum());
  fFullFitFnct->SetParameter(7, fDoubleGaussian->GetParameter((4)));
  fFullFitFnct->SetParameter(8, fDoubleGaussian->GetParameter((5)));
  fFullFitFnct->SetLineColor(kBlue);
}
void ForgivingFitter::CalculateBackgorund(TH1F* targetHisto, float massCutMin,
                                          float massCutMax) {
  TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)",
                                    fBkgRangeMin, fBkgRangeMax);
  fLambda_background->SetParameter(0, fFullFitFnct->GetParameter(0));
  fLambda_background->SetParameter(1, fFullFitFnct->GetParameter(1));
  fLambda_background->SetParameter(2, fFullFitFnct->GetParameter(2));
  fLambda_background->SetLineStyle(3);
  fLambda_background->SetLineColor(kBlue);

//  backgroundErr = fLambda_background->IntegralError(
//      lowerBound, upperBound, backgroundR->GetParams(),
//      backgroundR->GetCovarianceMatrix().GetMatrixArray())
//      / double(histo->GetBinWidth(1));
  fBackgroundCounts = fLambda_background->Integral(massCutMin, massCutMax)
      / double(targetHisto->GetBinWidth(1));
  targetHisto->GetListOfFunctions()->Add(fLambda_background);
}
