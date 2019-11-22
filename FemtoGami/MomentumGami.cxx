/*
 * MomentumGami.cxx
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#include "MomentumGami.h"
#include "TError.h"
#include "TF1.h"
#include "TRandom.h"

#include "TSVDUnfold.h"
MomentumGami::MomentumGami(float maxkStar)
    : fResolution(nullptr),
      fResProjection(),
      fToUnfold(nullptr),
      fMaxkStar(maxkStar) {
  // TODO Auto-generated constructor stub
  gRandom->SetSeed(0);
}

MomentumGami::~MomentumGami() {
}

void MomentumGami::SetResolution(TH2F* resoMatrix, float UnitConversion) {
  //User needs to make sure that the target binning is set properly!
  TString hname;
  fResolution = resoMatrix;
  std::cout << "resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion): "
            << resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion)
            << std::endl;
  for (int i = 0;
      i < resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion); i++) {
    hname = Form("hprojY%i", i);
    fResProjection.push_back(
        (TH1F*) resoMatrix->ProjectionY(hname, i + 1, i + 1));
    //normalize
    fResProjection[i]->Scale(
        1. / fResProjection[i]->Integral(1, fResProjection[i]->GetNbinsX()));
  }
}

void MomentumGami::Unfold(TH1F* InputDist) {
  if (InputDist->FindBin(fMaxkStar) > InputDist->GetNbinsX()) {
    Error("MomentumGami::Unfold",
          "Distribution has smaller kStar range than unfolding");
    return;
  }

  fToUnfold = (TH1F*) InputDist->Clone("ToUnfold");
  TF1 * momSmearing = new TF1("momSmearing", this, &MomentumGami::Eval, 0,
                              fMaxkStar, fToUnfold->FindBin(fMaxkStar) - 5,
                              "momSmearing", "momSmearing");  // create TF1 class.
  momSmearing->SetParameter(0, gRandom->Uniform(0.98, 1.3));
  momSmearing->SetParLimits(0, 0., 2.);
  momSmearing->SetParameter(1, gRandom->Uniform(0.98, 1.1));
  momSmearing->SetParLimits(1, 0.5, 1.5);
  for (int iPar = 2; iPar < fToUnfold->FindBin(fMaxkStar) - 5; ++iPar) {
    momSmearing->SetParameter(iPar, gRandom->Uniform(0.95, 1.05));
    momSmearing->SetParLimits(iPar, 0.5, 1.5);
  }
  std::cout << "Result: " << fToUnfold->Fit("momSmearing", "R") << std::endl;
  for (int iBim = 1; iBim < InputDist->FindBin(fMaxkStar) - 4; ++iBim) {
    int ParNmb = iBim - 1;
    InputDist->SetBinContent(
        iBim,
        momSmearing->GetParameter(ParNmb) * InputDist->GetBinContent(iBim));
    //keep the same relative error
    InputDist->SetBinError(
        iBim, momSmearing->GetParameter(ParNmb) * InputDist->GetBinError(iBim));
  }
  delete momSmearing;
  delete fToUnfold;
  return;
}

TH1F* MomentumGami::Fold(TH1F* InputDist) {
  TH1F* uncorr = (TH1F*) InputDist->Clone(
      TString::Format("%s_Refolded", InputDist->GetName()).Data());
  uncorr->Reset();
  const int nbinsProj = fToUnfold->FindBin(fMaxkStar);
  for (int iTRUE = 1; iTRUE < nbinsProj; iTRUE++) {
    for (int iRECO = 1; iRECO <= nbinsProj - 5; iRECO++) {
      float binCont = uncorr->GetBinContent(iRECO);
      binCont += fResProjection[iTRUE - 1]->GetBinContent(iRECO)
          * InputDist->GetBinContent(iTRUE);
      uncorr->SetBinContent(iRECO, binCont);
    }
  }
  return uncorr;
}

double MomentumGami::Eval(double *x, double *p) {
  if (fToUnfold->FindBin(x[0]) == 0
      || fToUnfold->FindBin(x[0]) > fToUnfold->GetNbinsX()) {
    Error("MomentumGami::operator()",
          "Trying to call with x-Value out of range");
    return -999;
  }
  //corrected "imaginary histo"
  const int nbinsProj = fToUnfold->FindBin(fMaxkStar);
  std::vector<float> newGuess;
  for (int i = 0; i < nbinsProj - 3; i++) {
    newGuess.push_back(p[i] * fToUnfold->GetBinContent(i + 1));
  }

//now build uncorrected: take zz "imaginary corrected" and smear it:
  std::vector<float> PossibleUncorrected;
  for (int ibin = 0; ibin < nbinsProj; ++ibin) {
    PossibleUncorrected.push_back(0.);
  }
  for (int iproj = 0; iproj < nbinsProj; iproj++) {
    for (int ibin = 0; ibin <= nbinsProj - 5; ibin++) {
      PossibleUncorrected.at(ibin) += fResProjection[iproj]->GetBinContent(
          ibin + 1) * newGuess[iproj];
    }
  }

  return PossibleUncorrected[fToUnfold->FindBin(x[0]) - 1];
}

TH1F* MomentumGami::UnfoldviaTSVD(TH1F* InputDist, TList* QA) {
  //Left here not to delete work. However does not seem to do the job.
  TH1D *xini = fResolution->ProjectionX("MC truth");
  TH1D *bini = fResolution->ProjectionY("MC reco");
  TH1D *oneDClone = (TH1D*) InputDist->Clone(
      TString::Format("%s_DClone", InputDist->GetName()));
  TH2D *statcov = new TH2D("statcov", "covariance matrix",
                           oneDClone->GetNbinsX(),
                           oneDClone->GetXaxis()->GetXmin(),
                           oneDClone->GetXaxis()->GetXmax(),
                           oneDClone->GetNbinsX(),
                           oneDClone->GetXaxis()->GetXmin(),
                           oneDClone->GetXaxis()->GetXmax());
  for (int i = 1; i <= oneDClone->GetNbinsX(); i++) {
    statcov->SetBinContent(
        i, i, oneDClone->GetBinError(i) * oneDClone->GetBinError(i));
  }

  TSVDUnfold *tsvdunf = new TSVDUnfold(oneDClone, statcov, bini, xini,
                                       (TH2D*)fResolution);
  // It is possible to normalise unfolded spectrum to unit area
  tsvdunf->SetNormalize(kFALSE);  // no normalisation here

  // Perform the unfolding with regularisation parameter kreg = 13
  // - the larger kreg, the finer grained the unfolding, but the more fluctuations occur
  // - the smaller kreg, the stronger is the regularisation and the bias
  TH1D* unfres = tsvdunf->Unfold(13);

  // Get the distribution of the d to cross check the regularization
  // - choose kreg to be the point where |d_i| stop being statistically significantly >>1
  TH1D* ddist = tsvdunf->GetD();

  // Get the distriutaucovbution of the singular values
  TH1D* svdist = tsvdunf->GetSV();

  // Compute the error matrix for the unfolded spectrum using toy MC
  // using the measured covariance matrix as input to generate the toys
  // 100 toys should usually be enough
  // The same method can be used for different covariance matrices separately.
  TH2D* ustatcov = tsvdunf->GetUnfoldCovMatrix(statcov, 100);

  // Now compute the error matrix on the unfolded distribution originating
  // from the finite detector matrix statistics
  TH2D* uadetcov = tsvdunf->GetAdetCovMatrix(100);

  // Sum up the two (they are uncorrelated)
  ustatcov->Add(uadetcov);

  // Get the computed regularized covariance matrix
  // (always corresponding to total uncertainty passed in constructor)
  // and add uncertainties from finite MC statistics.
  TH2D* utaucov = tsvdunf->GetXtau();
  utaucov->Add(uadetcov);

  //Get the computed inverse of the covariance matrix
  TH2D* uinvcov = tsvdunf->GetXinv();

  for (int i = 1; i <= unfres->GetNbinsX(); i++) {
    unfres->SetBinError(i, TMath::Sqrt(utaucov->GetBinContent(i, i)));
  }
  TH1D* QAratio = (TH1D*) unfres->Clone(
      TString::Format("%sQARatio", InputDist->GetName()).Data());
  QAratio->Divide(InputDist);

  QA->Add(
      ddist->Clone(TString::Format("%s_dDist", InputDist->GetName()).Data()));
  QA->Add(
      utaucov->Clone(
          TString::Format("%s_utaucov", InputDist->GetName()).Data()));
  QA->Add(QAratio);
  return (TH1F*) unfres;
}
