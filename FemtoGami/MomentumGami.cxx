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

MomentumGami::MomentumGami(float maxkStar)
    : fResolution(nullptr),
      fResProjection(),
      fToUnfold(nullptr),
      fMaxkStar(maxkStar) {
  // TODO Auto-generated constructor stub

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

void MomentumGami::Unfold(TH1F* InputDist, TH1F* OutputDist) {
  if (InputDist->FindBin(fMaxkStar) > InputDist->GetNbinsX()) {
    Error("MomentumGami::Unfold",
          "Distribution has smaller kStar range than unfolding");
    return;
  }

  fToUnfold = (TH1F*) InputDist->Clone("ToUnfold");
  TF1 * momSmearing = new TF1("momSmearing", this, &MomentumGami::Eval, 0,
                              fMaxkStar, fToUnfold->FindBin(fMaxkStar) - 3,
                              "momSmearing", "momSmearing");  // create TF1 class.
  for (int iPar = 0; iPar < fToUnfold->FindBin(fMaxkStar) - 3; ++iPar) {
    momSmearing->SetParameter(iPar, gRandom->Uniform(0.8, 1.2));
  }
  fToUnfold->Fit("momSmearing", "R");

  for (int iBim = 1; iBim < InputDist->FindBin(fMaxkStar) - 2; ++iBim) {
    int ParNmb = iBim - 1;
    OutputDist->SetBinContent(
        iBim,
        momSmearing->GetParameter(ParNmb) * InputDist->GetBinContent(iBim));
    OutputDist->SetBinError(iBim, InputDist->GetBinError(iBim));
  }
  for (int iBim = InputDist->FindBin(fMaxkStar) - 2;
      iBim < InputDist->GetNbinsX() + 1; ++iBim) {
    OutputDist->SetBinContent(iBim, InputDist->GetBinContent(iBim));
    OutputDist->SetBinError(iBim, InputDist->GetBinError(iBim));
  }
  delete momSmearing;
  delete fToUnfold;
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
  float zz[nbinsProj];
  for (int i = 0; i < nbinsProj - 3; i++) {
    zz[i] = p[i] * fToUnfold->GetBinContent(i + 1);
  }

//now build uncorrected: take zz "imaginary corrected" and smear it:
  float uncorr[fToUnfold->GetNbinsX()];
  for (int i = 0; i < fToUnfold->GetNbinsX(); i++) {
    uncorr[i] = 0.;
  }

  for (int iproj = 0; iproj < nbinsProj; iproj++) {
    for (int ibin = 1; ibin <= nbinsProj - 3; ibin++) {
      uncorr[ibin - 1] += fResProjection[iproj]->GetBinContent(ibin)
          * zz[ibin - 1];
    }
  }
  return uncorr[fToUnfold->FindBin(x[0]) - 1];
}

