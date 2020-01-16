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
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#include "TSVDUnfold.h"
MomentumGami::MomentumGami(float maxkStar)
    : fQAList(new TList()),
      fResolution(nullptr),
      fResProjection(),
      fToUnfold(nullptr),
      fTrainQA(nullptr),
      fResp(1),
      fIter(1),
      fResponseDefault(nullptr),
      fResponseLower(nullptr),
      fResponseUpper(nullptr),
      fMaxkStar(maxkStar),
      fUnitConversion(1.) {
  // TODO Auto-generated constructor stub
  gRandom->SetSeed(0);
  fQAList->SetName("MomentumGamiQA");
  fQAList->SetOwner();
}

MomentumGami::~MomentumGami() {
}

void MomentumGami::SetResolution(TH2F* resoMatrix, float UnitConversion) {
  //User needs to make sure that the target binning is set properly!
  TString hname;
  TList* QAResolution = new TList();
  QAResolution->SetName("QAResolution");
  QAResolution->SetOwner(true);
  fQAList->Add(QAResolution);
  fResolution = resoMatrix;
  int nBins = resoMatrix->GetNbinsX();
  double xmin = resoMatrix->GetXaxis()->GetXmin();
  double xmax = resoMatrix->GetXaxis()->GetXmax();
  fResponseDefault = new RooUnfoldResponse(nBins, xmin, xmax);
  fResponseLower = new RooUnfoldResponse(nBins, xmin, xmax);
  fResponseUpper = new RooUnfoldResponse(nBins, xmin, xmax);
  TH1F* Mean = new TH1F("MeanFit", "MeanFit", nBins, xmin, xmax);
  fQAList->Add(Mean);
  TH1F* MeanErr = new TH1F("MeanErrFit", "MeanErrFit", nBins, xmin, xmax);
  fQAList->Add(MeanErr);
  TH1F* Width = new TH1F("WidthFit", "WidthFit", nBins, xmin, xmax);
  fQAList->Add(Width);
  TH1F* WidthErr = new TH1F("WidthErrFit", "WidthErrFit", nBins, xmin, xmax);
  fQAList->Add(WidthErr);
  TF1* HalfAGauss =
      new TF1(
          "HalfAGauss",
          "[2]*(TMath::Gaus(x,[0],[1],0)+TMath::Gaus(x,-[0],[1],0))/TMath::Sqrt(2*TMath::Pi()*[1]*[1])",
          0, 1);
  TF1* HalfAWideGauss =
      new TF1(
          "HalfAWideGauss",
          "[2]*(TMath::Gaus(x,[0],[1],0)+TMath::Gaus(x,-[0],[1],0))/TMath::Sqrt(2*TMath::Pi()*[1]*[1])",
          0, 1);
  TF1* HalfANarrowGauss =
      new TF1(
          "HalfANarrowGauss",
          "[2]*(TMath::Gaus(x,[0],[1],0)+TMath::Gaus(x,-[0],[1],0))/TMath::Sqrt(2*TMath::Pi()*[1]*[1])",
          0, 1);

  HalfAGauss->SetNpx(5000);
  for (int i = 0;
      i < resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion); i++) {
    hname = Form("hprojY%i", i);
    fResProjection.push_back(
        (TH1F*) resoMatrix->ProjectionY(hname, i + 1, i + 1));
    int nEntries = fResProjection[i]->Integral(1,
                                               fResProjection[i]->GetNbinsX());
    hname = Form("SamplingQAY%i", i);
    TH1F* samplingQA = new TH1F(hname, hname, 1000, 0, 1);
    //normalize
    fResProjection[i]->Scale(1. / (float) nEntries);
    //Fit
    HalfAGauss->SetParameter(0, fResProjection[i]->GetMean());
    HalfAGauss->SetParameter(1, fResProjection[i]->GetRMS());
    HalfAGauss->SetParameter(
        2,
        fResProjection[i]->GetBinContent(
            fResProjection[i]->FindBin(fResProjection[i]->GetMean())));
    fResProjection[i]->Fit(HalfAGauss, "R");

    HalfAWideGauss->SetParameter(0, HalfAGauss->GetParameter(0));
    HalfAWideGauss->SetParameter(1, HalfAGauss->GetParameter(1) * 1.03);
    HalfAWideGauss->SetParameter(2, HalfAGauss->GetParameter(2));

    HalfANarrowGauss->SetParameter(0, HalfAGauss->GetParameter(0));
    HalfANarrowGauss->SetParameter(1, HalfAGauss->GetParameter(1) * .97);
    HalfANarrowGauss->SetParameter(2, HalfAGauss->GetParameter(2));

    Mean->SetBinContent(
        i + 1, HalfAGauss->GetParameter(0) - Mean->GetBinCenter(i + 1));
    Mean->SetBinError(i + 1, HalfAGauss->GetParError(0));
    MeanErr->SetBinContent(
        i + 1, HalfAGauss->GetParError(0) / HalfAGauss->GetParameter(0));
    Width->SetBinContent(i + 1, HalfAGauss->GetParameter(1));
    Width->SetBinError(i + 1, HalfAGauss->GetParError(1));
    WidthErr->SetBinContent(
        i + 1, HalfAGauss->GetParError(1) / HalfAGauss->GetParameter(1));
    //Train Responses
    float kTrue = resoMatrix->GetXaxis()->GetBinCenter(i + 1);
    for (int nFills = 0; nFills < fResProjection[i]->GetEntries(); ++nFills) {
      float kMeas = HalfAGauss->GetRandom();
      fResponseDefault->Fill(kMeas, kTrue);
      samplingQA->Fill(kMeas);
      float kMeasLower = HalfANarrowGauss->GetRandom();
      fResponseLower->Fill(kMeasLower, kTrue);
      float kMeasUpper = HalfAWideGauss->GetRandom();
      fResponseUpper->Fill(kMeasUpper, kTrue);
    }
    samplingQA->Scale(1. / samplingQA->GetEntries());
    QAResolution->Add(fResProjection[i]);
    QAResolution->Add(samplingQA);
  }
  TH2D* responsematrix = fResponseDefault->HresponseNoOverflow();
  responsematrix->SetName(TString::Format("%sRespMatrix", "").Data());
  TH2D* responsematrixUp = fResponseUpper->HresponseNoOverflow();
  responsematrixUp->SetName(TString::Format("%sRespMatrixUp", "").Data());
  TH2D* responsematrixLow = fResponseLower->HresponseNoOverflow();
  responsematrixLow->SetName(TString::Format("%sRespMatrixLow", "").Data());
  fQAList->Add(responsematrix);
  fQAList->Add(responsematrixUp);
  fQAList->Add(responsematrixLow);
  fUnitConversion = UnitConversion;
}

void MomentumGami::UnfoldGuessing(TH1F* InputDist) {
  std::cout << "Unfold Guessing \n";
  if (InputDist->FindBin(0.3) > InputDist->GetNbinsX()) {
    Error("MomentumGami::Unfold",
          "Distribution has smaller kStar range than unfolding");
    return;
  }

  fToUnfold = (TH1F*) InputDist->Clone("ToUnfold");
  TF1 * momSmearing = new TF1("momSmearing", this, &MomentumGami::Eval, 0, 0.3,
                              fToUnfold->FindBin(0.3) - 5);  // create TF1 class.
  momSmearing->SetParameter(0, gRandom->Uniform(0.98, 1.3));
  momSmearing->SetParLimits(0, 0., 2.);
  momSmearing->SetParameter(1, gRandom->Uniform(0.98, 1.1));
  momSmearing->SetParLimits(1, 0.5, 1.5);
  for (int iPar = 2; iPar < fToUnfold->FindBin(0.3) - 5; ++iPar) {
    momSmearing->SetParameter(iPar, gRandom->Uniform(0.95, 1.05));
    momSmearing->SetParLimits(iPar, 0.5, 1.5);
  }
  std::cout << "Starting Fit \n";
  fToUnfold->Fit("momSmearing", "R");
  TH1F* parameters = new TH1F("QAGuessing", "QAGuessing",
                              fToUnfold->FindBin(0.3) - 5, 0, 0.3);
  fQAList->Add(parameters);
  for (int iBim = 1; iBim < InputDist->FindBin(0.3) - 10; ++iBim) {
    int ParNmb = iBim - 1;
    parameters->SetBinContent(iBim, momSmearing->GetParameter(ParNmb));
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
  const int nbinsProj = InputDist->FindBin(fMaxkStar);
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
  const int nbinsProj = fToUnfold->FindBin(0.3);
  std::vector<float> newGuess;
  if (newGuess.size() > 0) {
    newGuess.resize(0);
  }
  for (int i = 0; i < nbinsProj; i++) {
    if (i < nbinsProj - 5) {
      newGuess.push_back(p[i] * fToUnfold->GetBinContent(i + 1));
    } else {
      newGuess.push_back(fToUnfold->GetBinContent(i + 1));
    }
  }

//now build uncorrected: take zz "imaginary corrected" and smear it:
  std::vector<float> PossibleUncorrected;
  if (PossibleUncorrected.size() > 0) {
    PossibleUncorrected.resize(0);
  }
  for (int ibin = 0; ibin < nbinsProj; ++ibin) {
    PossibleUncorrected.push_back(0.);
  }
  for (int ibin = 0; ibin <= nbinsProj - 5; ibin++) {
    for (int iproj = 0; iproj < nbinsProj; iproj++) {
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
                                       (TH2D*) fResolution);
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

TH1F* MomentumGami::UnfoldviaRooResp(TH1F* InputDist, double Rescaling) {
//Its up to the user to ensure the same binning, this is just a simple check
  if (fResolution->GetXaxis()->GetBinWidth(1) / fUnitConversion
      != InputDist->GetXaxis()->GetBinWidth(1)) {
    std::cout
        << "MomentumGami::UnfoldviaRooResp: Bin width of Resolution matrix "
        << "different than bin widht of the input dist! \n"
        << "Histname of the input in question: " << InputDist->GetName()
        << std::endl;
    std::cout << "fResolution->GetXaxis()->GetBinWidth(1)/fUnitConversion: "
              << fResolution->GetXaxis()->GetBinWidth(1) / fUnitConversion
              << std::endl;
    std::cout << "InputDist->GetXaxis()->GetBinWidth(1): "
              << InputDist->GetXaxis()->GetBinWidth(1) << std::endl;
  }

//only unfold within the region the resolution matrix is defined.
  TString toUnfoldName = TString::Format("%s_ToUnfold", InputDist->GetName());
  RooUnfoldResponse* resp[3] = { fResponseLower, fResponseDefault,
      fResponseUpper };

  RooUnfoldBayes unfolderer(resp[fResp], InputDist, fIter);
//  RooUnfoldBinByBin unfolderer(resp[fResp], InputDist);
  TH1F* unfoldedHist = (TH1F*) unfolderer.Hreco(RooUnfold::kCovariance);  //RooUnfold::kCovToy
  unfolderer.Print();
  TH1F* Ratio = (TH1F*) unfoldedHist->Clone(
      TString::Format("%s_ratioUnfolded", InputDist->GetName()));
  Ratio->Divide(InputDist);

  TH1F* unfoldedQA = (TH1F*) unfoldedHist->Clone(
      TString::Format("%s_BayesUnfolded_nIter%u", InputDist->GetName(), fIter)
          .Data());
  TH1F* inputRelErrQA = (TH1F*) InputDist->Clone(
      TString::Format("%s_RelErrInput_nIter%u", InputDist->GetName(), fIter)
          .Data());
  TH1F* unfoldedRelErrQA = (TH1F*) unfoldedHist->Clone(
      TString::Format("%s_RelErrBayesUnfolded_nIter%u", InputDist->GetName(),
                      fIter).Data());
  for (int iBin = 1; iBin <= unfoldedRelErrQA->GetNbinsX(); ++iBin) {
    // The statistical error is set to the Sqrt(Entries) - if the mixed event
    // is already rescaled e.g. due to reweighting, it needs to be accounted
    // for when taking the error!
    unfoldedHist->SetBinError(
        iBin,
        TMath::Sqrt(unfoldedHist->GetBinContent(iBin))
            / TMath::Sqrt(Rescaling));

    if (TMath::Abs(inputRelErrQA->GetBinContent(iBin)) > 1e-5) {
      inputRelErrQA->SetBinContent(
          iBin,
          inputRelErrQA->GetBinError(iBin)
              / inputRelErrQA->GetBinContent(iBin));
    } else {
      inputRelErrQA->SetBinContent(iBin, 0);
    }
    inputRelErrQA->SetBinError(iBin, 0);
    if (TMath::Abs(unfoldedRelErrQA->GetBinContent(iBin)) > 1e-5) {
      unfoldedRelErrQA->SetBinContent(
          iBin,
          unfoldedRelErrQA->GetBinError(iBin)
              / unfoldedRelErrQA->GetBinContent(iBin));
    } else {
      unfoldedRelErrQA->SetBinContent(iBin, 0);
    }
    unfoldedRelErrQA->SetBinError(iBin, 0);
  }
  fQAList->Add(unfoldedQA);
  fQAList->Add(unfoldedRelErrQA);
  fQAList->Add(inputRelErrQA);
  fQAList->Add(Ratio);
  TString outName = TString::Format("%s_Unfolded", InputDist->GetName());
  return (TH1F*) unfoldedHist->Clone(outName.Data());
}

void MomentumGami::TrainRooResponse(TH2F* momMatrix,
                                    RooUnfoldResponse* roo_resp) {
  fTrainQA = new TH2F("QATrainRooResponse", "QATrainRooResponse",
                      momMatrix->GetNbinsX(),
                      momMatrix->GetXaxis()->GetXmin() / fUnitConversion,
                      momMatrix->GetXaxis()->GetXmax() / fUnitConversion,
                      momMatrix->GetNbinsY(),
                      momMatrix->GetYaxis()->GetXmin() / fUnitConversion,
                      momMatrix->GetXaxis()->GetXmax() / fUnitConversion);
  for (Int_t iBinx = 1; iBinx < momMatrix->GetNbinsX() + 1; iBinx++) {
    double ksTrue = momMatrix->GetXaxis()->GetBinCenter(iBinx)
        / fUnitConversion;
    for (Int_t iBiny = 1; iBiny < momMatrix->GetNbinsY() + 1; iBiny++) {
      double ksMeas = momMatrix->GetYaxis()->GetBinCenter(iBiny)
          / fUnitConversion;
      int nEntries = momMatrix->GetBinContent(iBinx, iBiny);
      fTrainQA->Fill(ksTrue,ksMeas);
      for (int iFill = 0; iFill < nEntries; iFill++) {
        roo_resp->Fill(ksTrue, ksMeas);
      }
    }
  }
  fQAList->Add(fTrainQA);
}
