#include "TidyCats.h"
#include "TRandom3.h"
#include "SidebandSigma.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "TDatabasePDG.h"
#include <iostream>
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TSystem.h"

float GetkStar(const TLorentzVector &Part1Momentum,
               const TLorentzVector &Part2Momentum) {
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  Part1Momentum.M());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  Part2Momentum.M());
  trackSum = SPtrack + TPProng;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  SPtrackCMS = SPtrack;
  TPProngCMS = TPProng;

  SPtrackCMS.Boost(-betax, -betay, -betaz);
  TPProngCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5 * trackRelK.P();
  return results * 1000.;  // in MeV
}

TH2F* ComputeSmearingMatrix(TString &InputDir, TString &trigger,
                            TString &suffix, int nParticles) {
  auto file = TFile::Open(Form("%s/AnalysisResults.root", InputDir.Data()));
  auto filename = TString::Format("%s/SmearSideband.root", InputDir.Data());

  if (TFile::Open(filename)) {
    std::cout << "Using smearing matrix from existing file\n";
    auto outfile = TFile::Open( filename);
    auto histSmear = (TH2F*) outfile->Get("histSmear");
    outfile->Close();
    return histSmear;
  } else {
    std::cout << "No smearing matrix found - recomputing \n";
    auto outfile = new TFile(filename, "RECREATE");

    TLorentzVector protonPart;
    TString protonName = trigger;
    protonName += "TrackCuts";
    protonName += suffix;
    TDirectory *protonDir = file->GetDirectory(protonName);
    auto protonList = (TList *) protonDir->Get(protonName);
    protonList = (TList *) protonList->FindObject("after");
    auto protonpT = (TH1F*) protonList->FindObject("pTDist_after");
    auto protoneta = (TH1F*) protonList->FindObject("EtaDist_after");
    auto protonphi = (TH1F*) protonList->FindObject("phiDist_after");
    auto protonMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    outfile->cd();
    protonpT->Write("protonpT");
    protoneta->Write("protonEta");
    protonphi->Write("protonPhi");

    TLorentzVector lambdaPart;
    TString lambdaName = trigger;
    lambdaName += "v0Cuts";
    lambdaName += suffix;
    TDirectory *lambdaDir = file->GetDirectory(lambdaName);
    auto lambdaList = (TList *) lambdaDir->Get(lambdaName);
    lambdaList = (TList *) lambdaList->FindObject("v0Cuts");
    lambdaList = (TList *) lambdaList->FindObject("after");
    auto lambdapT = (TH1F*) lambdaList->FindObject("pTDist_after");
    auto lambdaeta = (TH1F*) lambdaList->FindObject("EtaDist_after");
    auto lambdaphi = (TH1F*) lambdaList->FindObject("PhiDist_after");
    auto lambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    outfile->cd();
    lambdapT->Write("lambdapT");
    lambdaeta->Write("lambdaEta");
    lambdaphi->Write("lambdaPhi");

    TLorentzVector photonPart;
    TString photonName = trigger;
    photonName += "PhotonCuts";
    photonName += suffix;
    TDirectory *photonDir = file->GetDirectory(photonName);
    auto photonList = (TList *) photonDir->Get(photonName);
    auto photonpT = (TH1F*) photonList->FindObject("fHistV0Pt");
    auto photonetaphi = (TH2F*) photonList->FindObject("fHistEtaPhi");
    auto photoneta = (TH1F*) photonetaphi->ProjectionX();
    auto photonphi = (TH1F*) photonetaphi->ProjectionY();
    outfile->cd();
    photonpT->Write("photonpT");
    photoneta->Write("photonEta");
    photonphi->Write("photonPhi");

    TString sigmaName = trigger;
    sigmaName += "Sigma0Cuts";
    sigmaName += suffix;
    TDirectory *sigmaDir = file->GetDirectory(sigmaName);
    auto sigmaList = (TList *) sigmaDir->Get(sigmaName);
    auto sigmaMassRec = (TH1F*) sigmaList->FindObject("fHistInvMass");
    sigmaMassRec->Write("sigmaMassRec");
    float entries = sigmaMassRec->Integral(sigmaMassRec->FindBin(1.2),
                                           sigmaMassRec->FindBin(1.25));
    auto sigmaMass = TDatabasePDG::Instance()->GetParticle(3212)->Mass();

    TLorentzVector sigmaPart;
    auto histSigma = new TH1F("histSigma", "", 1000, 1, 2);
    auto histSmear =
        new TH2F(
            "histSmear",
            "; #it{k}*_{p#minus#Sigma^{0}} (GeV/#it{c}); #it{k}*_{p#minus#Lambda} (GeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);
    auto relMom = new TNtuple("relMom", "relMom", "pSigma:pLambda:pPhoton");

    float counter = 0;
    int modulo = (nParticles > 10000) ? nParticles / 10000 : 1;
    int i = 0;
    float kstarpSigma, kstarpLambda;
    float nCount = nParticles;
    while (i < nParticles) {
      protonPart.SetPtEtaPhiM(protonpT->GetRandom(), protoneta->GetRandom(),
                              protonphi->GetRandom(), protonMass);
      lambdaPart.SetPtEtaPhiM(lambdapT->GetRandom(), lambdaeta->GetRandom(),
                              lambdaphi->GetRandom(), lambdaMass);
      photonPart.SetPtEtaPhiM(photonpT->GetRandom(), photoneta->GetRandom(),
                              photonphi->GetRandom(), 0);
      sigmaPart = photonPart + lambdaPart;
      histSigma->Fill(sigmaPart.M());
      if (std::abs(sigmaPart.M() - sigmaMass) > 0.005) {
        continue;
      }

      kstarpSigma = GetkStar(protonPart, sigmaPart);
      kstarpLambda = GetkStar(protonPart, lambdaPart);

      relMom->Fill(kstarpSigma, kstarpLambda, GetkStar(protonPart, photonPart));
      histSmear->Fill(kstarpSigma, kstarpLambda);
      ++i;
      if ((i % 100) == 0)
        std::cout << "\r " << 100. / nCount * i << "%";
    }
    histSigma->Scale(
        entries
            / histSigma->Integral(histSigma->FindBin(1.2),
                                  histSigma->FindBin(1.25)));
    histSigma->Write();
    histSmear->Write();
    relMom->Write();
    outfile->Close();
    return histSmear;
  }
}

TGraph *GetSmearedCF(TGraph* CF, TH2F* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF);
  }
  return smearedCF;
}

int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");

  TString InputDir = argv[1];
  TString trigger = argv[2];
  TString suffix = argv[3];
  int nParticles = atoi(argv[4]);

  auto filename = TString::Format("%s/SherlockSideband.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");

  TRandom3 rangen(0);
  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

  auto side = new SidebandSigma();
  side->SetRebin(10);
  side->SetSideBandFile(InputDir.Data(), trigger.Data(), suffix.Data());
  const double sidebandNormDown = 250;
  const double sidebandNormUp = 400;
  side->SetNormalizationRange(sidebandNormDown, sidebandNormUp);

  side->SideBandCFs();
  auto SBmerge = side->GetSideBands(5);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Get the smearing matrix
  auto histSmear = ComputeSmearingMatrix(InputDir, trigger, suffix, nParticles);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                       suffix.Data());
  CATSinput->ObtainCFs(10, 250, 400);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  int NumMomBins = 25;
  double kMin = -9.99;
  double kMax = 490.01;

  // pp radius systematic variations
  const double ppRadius = 1.154;

  DLM_Ck* Ck_pL = new DLM_Ck(1, 4, NumMomBins, kMin, kMax,
                             Lednicky_SingletTriplet);

  // NLO
  Ck_pL->SetPotPar(0, 2.91);
  Ck_pL->SetPotPar(1, 2.78);
  Ck_pL->SetPotPar(2, 1.54);
  Ck_pL->SetPotPar(3, 2.72);

  // LO
//  Ck_pL->SetPotPar(0, 1.91);
//  Ck_pL->SetPotPar(1, 1.4);
//  Ck_pL->SetPotPar(2, 1.23);
//  Ck_pL->SetPotPar(3, 2.13);

  Ck_pL->SetSourcePar(0, ppRadius);
  Ck_pL->Update();

  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  const double lambdaPurity = 0.95;
  const double lambdaPrimary = 0.619493;
  const Particle lambda(lambdaPurity, lambdaPrimary,
                        { { 1.f - lambdaPrimary } });

  const CATSLambdaParam lambdaParam(proton, lambda);

  // Lambda parameters and momentum resolution
  DLM_CkDecomposition CkDec_pL("pL", 1, *Ck_pL, CATSinput->GetSigmaFile(1));
  CkDec_pL.AddContribution(
      0, 1. - lambdaParam.GetLambdaParam(CATSLambdaParam::Primary),
      DLM_CkDecomposition::cFake);
  CkDec_pL.Update();

  auto grSideband = new TGraph();
  auto grSidebandRaw = new TGraph();
  for (unsigned int i = 0; i < Ck_pL->GetNbins(); ++i) {
    grSidebandRaw->SetPoint(i, Ck_pL->GetBinCenter(0, i),
                            Ck_pL->Eval(Ck_pL->GetBinCenter(0, i)));
    grSideband->SetPoint(i, Ck_pL->GetBinCenter(0, i),
                         CkDec_pL.EvalCk(Ck_pL->GetBinCenter(0, i)));

  }

  auto grSidebandSmeared = GetSmearedCF(grSideband, histSmear);

  // check with inverted axes
  auto histSmearInv = (TH2F*)histSmear->Clone("histSmearInv");
  for (int i = 0; i<histSmear->GetNbinsX(); ++i) {
    for (int j = 0; j<histSmear->GetNbinsY(); ++j) {
      histSmearInv->SetBinContent(j, i, histSmear->GetBinContent(i, j));
    }
  }
  auto grSidebandSmearedInv = GetSmearedCF(grSideband, histSmearInv);


//  auto c = new TCanvas();
//  SBmerge->Draw();
//  SBmerge->SetMaximum(1.7);
//  SBmerge->GetXaxis()->SetRangeUser(0, 350);
//  grSideband->Draw("L3 same");
//  c->Print("Sideband_fit.pdf");
//
//  auto d = new TCanvas();
//  histLambdaGamma->Draw("colz");
//  histLambdaGamma->GetXaxis()->SetRangeUser(0, 0.5);
//  histLambdaGamma->GetYaxis()->SetRangeUser(0, 0.5);
//  d->Print("momRes.pdf");
  outfile->cd();

  grSidebandRaw->Write("p-Lambda raw");
  grSideband->Write("p-Lambda smeared");
  grSidebandSmeared->Write("p-Lambda smeared with Photon");
  grSidebandSmearedInv->Write("p-Lambda ing");
  SBmerge->Write();
  CATSinput->GetSigmaFile(1)->Write();

  outfile->Close();
}
