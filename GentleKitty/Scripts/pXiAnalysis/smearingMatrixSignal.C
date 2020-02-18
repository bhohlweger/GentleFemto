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

  auto file = TFile::Open(Form("%s/AnalysisResultsTMP.root", InputDir.Data()));
  if (!file) {
    std::cout << "no input file \n";
    return -1;
  }
  auto filename = TString::Format("%s/SmearSideband.root", gSystem->pwd());
  std::cout << "Filename: " << filename << std::endl;
  if (TFile::Open(filename)) {
    std::cout << "Using smearing matrix from existing file\n";
    auto outfile = TFile::Open(filename);
    auto histSmear = (TH2D*) outfile->Get("histSmearSideband");
    outfile->Close();
    return 0;
  } else {
    std::cout << "No smearing matrix found - recomputing \n";
    auto outfile = new TFile(filename, "RECREATE");

    TLorentzVector protonPart;
    TString protonName = trigger;
    protonName += "TrackCuts";
    protonName += suffix;
    TDirectory *protonDir = file->GetDirectory(protonName);
    if (!protonDir) {
      std::cout << "no proton dir \n";
      return -1;
    }
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

    TString CascName = trigger;
    CascName += "CascadeCuts";
    CascName += suffix;
    TDirectory *cascDir = file->GetDirectory(CascName);
    if (!cascDir) {
      std::cout << "no casc dir\n";
      return -1;
    }

    auto cascList = (TList *) cascDir->Get(CascName);
    cascList = (TList *) cascList->FindObject("Cascade");
    cascList = (TList *) cascList->FindObject("after");

    TLorentzVector lambdaPart;
    auto lambdapT = (TH1D*) cascList->FindObject("v0Pt_after");
    auto lambdaeta = (TH1D*) cascList->FindObject("v0Eta_after");
    auto lambdaphi = (TH1D*) cascList->FindObject("v0Phi_after");
    auto lambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    outfile->cd();
    lambdapT->Write("lambdapT");
    lambdaeta->Write("lambdaEta");
    lambdaphi->Write("lambdaPhi");

    auto piList = (TList *) cascDir->Get(CascName);
    piList = (TList *) piList->FindObject("BachelorCuts");
    piList = (TList *) piList->FindObject("after");

    TLorentzVector PiPart;
    TLorentzVector KaonPart;
    auto pipT = (TH1D*) piList->FindObject("pTDist_after");
    auto pieta = (TH1D*) piList->FindObject("phiDist_after");
    auto piphi = (TH1D*) piList->FindObject("EtaDist_after");
    auto piMass = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    auto kaMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();

    outfile->cd();
    pipT->Write("photonpT");
    pieta->Write("photonEta");
    piphi->Write("photonPhi");

    auto XiMass = TDatabasePDG::Instance()->GetParticle(3312)->Mass();

    TLorentzVector XiPart;
    TLorentzVector HypoXiPart;
    auto histSigma = new TH1D("histXi", "", 1000, 1, 2);
    auto histSigmaPt = new TH1D("histXiPt", "", 1000, 0, 10);
    auto histSmear =
        new TH2D(
            "histSmear",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Xi^{-}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearPion =
        new TH2D(
            "histSmearPi",
            "; #it{k}*_{p#minus#pi} (MeV/#it{c}); #it{k}*_{p#minus#Xi^{-}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearSideband =
        new TH2D(
            "histSmearSideband",
            "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Xi^{-}} (MeV/#it{c})",
            500, 0, 1000, 500, 0, 1000);

    auto histSmearPionSideband =
        new TH2D(
            "histSmearPiSideband",
            "; #it{k}*_{p#minus#pi} (MeV/#it{c}); #it{k}*_{p#minus#Xi^{-}} (MeV/#it{c})",
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
      PiPart.SetPtEtaPhiM(pipT->GetRandom(), pieta->GetRandom(),
                          piphi->GetRandom(), piMass);
      KaonPart.SetPtEtaPhiM(PiPart.Pt(), PiPart.Eta(), PiPart.Phi(), kaMass);
      HypoXiPart = KaonPart + lambdaPart;
      if (TMath::Abs(HypoXiPart.M() - 1.672) < 0.005) {
        continue;
      }
      XiPart = PiPart + lambdaPart;
      massCurrentSigma = XiPart.M();
      histSigma->Fill(massCurrentSigma);

      histSigmaPt->Fill(XiPart.Pt());

      kstarpSigma = GetkStar(protonPart, XiPart);
      kstarpLambda = GetkStar(protonPart, lambdaPart);
      kstarpPhoton = GetkStar(protonPart, PiPart);


      // realistic sideband
      if (std::abs(massCurrentSigma - 1.346) < 0.014
          || std::abs(massCurrentSigma - 1.298) < 0.014) {
        histSmearSideband->Fill(kstarpLambda, kstarpSigma);
        histSmearPionSideband->Fill(kstarpPhoton, kstarpSigma);
      }

      if (std::abs(massCurrentSigma - XiMass) > 0.005) {
        continue;
      }

      if (i % 1000 == 0)
        relMom->Fill(kstarpSigma, kstarpLambda, GetkStar(protonPart, PiPart));
      histSmear->Fill(kstarpLambda, kstarpSigma);
      histSmearPion->Fill(kstarpPhoton, kstarpSigma);
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

    histSigma->Write();
    histSigmaPt->Write();
    histSmear->Write();
    histSmearPion->Write();
    histSmearSideband->Write();
    histSmearPionSideband->Write();
    relMom->Write();
    outfile->Close();
    return 0;
  }
}
