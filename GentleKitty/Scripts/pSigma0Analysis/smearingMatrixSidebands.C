#include "TidyCats.h"
#include "TRandom3.h"
#include "SidebandSigma.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "TH1D.h"
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
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TStopwatch.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TStyle.h"

/// =====================================================================================
float GetkStar(const TLorentzVector &Part1Momentum,
               const TLorentzVector &Part2Momentum) {
  static float results = 0.;
  static TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
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

  static TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5 * trackRelK.P();
  return results * 1000.;  // in MeV
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  TString trigger = argv[2];
  TString suffix = argv[3];
  int nParticles = atoi(argv[4]);

  auto file = TFile::Open(Form("%s/AnalysisResults.root", InputDir.Data()));
  auto filename = TString::Format("%s/SmearSideband.root", InputDir.Data());

  if (TFile::Open(filename)) {
    std::cout << "Sure you want to overwrite?\n";
    return 0;
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
    auto protonpT = (TH1D*) protonList->FindObject("pTDist_after");
    auto protoneta = (TH1D*) protonList->FindObject("EtaDist_after");
    auto protonphi = (TH1D*) protonList->FindObject("phiDist_after");
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
    auto lambdapT = (TH1D*) lambdaList->FindObject("pTDist_after");
    auto lambdaeta = (TH1D*) lambdaList->FindObject("EtaDist_after");
    auto lambdaphi = (TH1D*) lambdaList->FindObject("PhiDist_after");
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
    auto photonpT = (TH1D*) photonList->FindObject("fHistV0Pt");
    auto photonetaphi = (TH2D*) photonList->FindObject("fHistEtaPhi");
    auto photoneta = (TH1D*) photonetaphi->ProjectionX();
    auto photonphi = (TH1D*) photonetaphi->ProjectionY();
    outfile->cd();
    photonpT->Write("photonpT");
    photoneta->Write("photonEta");
    photonphi->Write("photonPhi");

    TLorentzVector trueSigmaPart;
    TString sigmaName = trigger;
    sigmaName += "Sigma0Cuts";
    sigmaName += suffix;
    TDirectory *sigmaDir = file->GetDirectory(sigmaName);
    auto sigmaList = (TList *) sigmaDir->Get(sigmaName);
    auto sigmaMassRec = (TH1D*) sigmaList->FindObject("fHistInvMass");
    auto sigmapT = (TH1D*) sigmaList->FindObject("fHistMassCutPt");
    auto sigmaetaphi = (TH2D*) sigmaList->FindObject("fHistEtaPhi");
    auto sigmaeta = (TH1D*) sigmaetaphi->ProjectionX();
    auto sigmaphi = (TH1D*) sigmaetaphi->ProjectionY();
    sigmaMassRec->Write("sigmaMassRec");
    float entries = sigmaMassRec->Integral(sigmaMassRec->FindBin(1.2),
                                           sigmaMassRec->FindBin(1.25));
    auto sigmaMass = TDatabasePDG::Instance()->GetParticle(3212)->Mass();
    const float purity = 0.346;

    TLorentzVector sigmaPart;
    auto histSigma = new TH1D("histSigma", "", 1000, 1, 2);
    auto histSigmaSigma = new TH1D("histSigmaSigma", "", 1000, 1, 2);
    auto histSigmaPt = new TH1D("histSigmaPt", "", 1000, 0, 10);
    auto histSmear =
        new TH2D(
            "histSmear",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearPhoton =
        new TH2D(
            "histSmearPhoton",
            "; #it{k}*_{p#minus#gamma} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearSigma =
        new TH2D(
            "histSmearSigma",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearSideband =
        new TH2D(
            "histSmearSideband",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearFull100 =
        new TH2D(
            "histSmearFull100",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearFull200 =
        new TH2D(
            "histSmearFull200",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearFull =
        new TH2D(
            "histSmearFull",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto relMom = new TNtuple("relMom", "relMom", "pSigma:pLambda:pPhoton");

    double counter = 0;
    int modulo = (nParticles > 10000) ? nParticles / 10000 : 1;
    int i = 0;
    double kstarpSigma, kstarpLambda, kstarpPhoton, kstarpSigmaSigma,
        kstarpLambdaSigma, massCurrentSigma;
    double nCount = nParticles;
    TGenPhaseSpace eventSigma0;
    TLorentzVector lambdaSigmaDecay;
    double sigmaDecay[2] = { lambdaMass, 0. };
    TStopwatch watch;
    while (i < nParticles) {
      protonPart.SetPtEtaPhiM(protonpT->GetRandom(), protoneta->GetRandom(),
                              protonphi->GetRandom(), protonMass);
      lambdaPart.SetPtEtaPhiM(lambdapT->GetRandom(), lambdaeta->GetRandom(),
                              lambdaphi->GetRandom(), lambdaMass);
      photonPart.SetPtEtaPhiM(photonpT->GetRandom(), photoneta->GetRandom(),
                              photonphi->GetRandom(), 0);
      sigmaPart = photonPart + lambdaPart;
      massCurrentSigma = sigmaPart.M();
      histSigma->Fill(massCurrentSigma);

      if (sigmaPart.Pt() < 1)
        continue;
      histSigmaPt->Fill(sigmaPart.Pt());

      kstarpSigma = GetkStar(protonPart, sigmaPart);
      kstarpLambda = GetkStar(protonPart, lambdaPart);
      kstarpPhoton = GetkStar(protonPart, photonPart);
      histSmearFull->Fill(kstarpLambda, kstarpSigma);

      if (std::abs(massCurrentSigma - sigmaMass) < 0.1) {
        histSmearFull100->Fill(kstarpLambda, kstarpSigma);
      }
      if (std::abs(massCurrentSigma - sigmaMass) < 0.2) {
        histSmearFull200->Fill(kstarpLambda, kstarpSigma);
      }

      // realistic sideband
      if (std::abs(massCurrentSigma - sigmaMass) < 0.05
          && std::abs(massCurrentSigma - sigmaMass) > 0.005) {
        histSmearSideband->Fill(kstarpLambda, kstarpSigma);
      }

      if (std::abs(massCurrentSigma - sigmaMass) > 0.003) {
        continue;
      }

      // let's add a true sigma
      histSigmaSigma->Fill(massCurrentSigma);
      if (gRandom->Uniform() < purity) {
        trueSigmaPart.SetPtEtaPhiM(sigmapT->GetRandom(), sigmaeta->GetRandom(),
                                   sigmaphi->GetRandom(),
                                   gRandom->Gaus(sigmaMass, 0.0015));
        histSigmaSigma->Fill(trueSigmaPart.M());

        eventSigma0.SetDecay(trueSigmaPart, 2, sigmaDecay);
        eventSigma0.Generate();
        lambdaSigmaDecay = *(eventSigma0.GetDecay(0));
        kstarpSigmaSigma = GetkStar(protonPart, trueSigmaPart);
        kstarpLambdaSigma = GetkStar(protonPart, lambdaSigmaDecay);
        histSmearSigma->Fill(kstarpLambdaSigma, kstarpSigmaSigma);
      }

      if (i % 1000 == 0)
        relMom->Fill(kstarpSigma, kstarpLambda,
                     GetkStar(protonPart, photonPart));
      histSmear->Fill(kstarpLambda, kstarpSigma);
      histSmearPhoton->Fill(kstarpPhoton, kstarpSigma);
      histSmearSigma->Fill(kstarpLambda, kstarpSigma);
      ++i;
      if ((i % 100) == 0) {
        watch.Stop();
        std::cout
            << "\r "
            << 100. / nCount * i
            << " % - running "
            << watch.RealTime() / 60.
            << "min - E.T.A. "
            << (watch.RealTime() / (100. / nCount * i) * 100 - watch.RealTime())
                / 60
            << "min";
        watch.Continue();
      }
    }

    histSigma->Scale(
        entries
            / histSigma->Integral(histSigma->FindBin(1.2),
                                  histSigma->FindBin(1.25)));

    histSigma->Write();
    histSigmaSigma->Write();
    histSigmaPt->Write();
    histSmear->Write();
    histSmearPhoton->Write();
    histSmearSigma->Write();
    histSmearFull100->Write();
    histSmearFull200->Write();
    histSmearFull->Write();
    histSmearSideband->Write();
    relMom->Write();
    outfile->Close();
    return 0;
  }
}
