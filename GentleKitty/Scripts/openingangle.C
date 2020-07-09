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


float GetOpeningAngle(const TLorentzVector &Part1Momentum,
               const TLorentzVector &Part2Momentum, float masscurrentphi) {
  static float Theta = 0.;
  float c=TMath::C();
  float m1=Part1Momentum.M();
  float m2=Part2Momentum.M(); //in gev
//  float p1=std::sqrt(Part1Momentum.X()*Part1Momentum.X()+Part1Momentum.Y()*Part1Momentum.Y()+Part1Momentum.Z()*Part1Momentum.Z());
//  float p2=std::sqrt(Part2Momentum.X()*Part2Momentum.X()+Part2Momentum.Y()*Part2Momentum.Y()+Part2Momentum.Z()*Part2Momentum.Z());
  float p1=Part1Momentum.P();
  float p2=Part2Momentum.P();
  float E1=Part1Momentum.E();
  float E2=Part2Momentum.E();
  float T1=Part1Momentum.Theta();
  float T2=Part2Momentum.Theta();
  Theta=180/TMath::Pi()*TMath::ACos((masscurrentphi*masscurrentphi-m1*m1-m2*m2*-2*E1*E2)/(-2*p1*p2));
  return Theta;


}


/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  const char* trigger = argv[2];
 // TString suffix = argv[3];
  int suffix= atoi(argv[3]);
  int nParticles = atoi(argv[4]);
  TList *mainDir;
  auto file = TFile::Open(Form("%s/AnalysisResults.root", InputDir.Data()));
  auto filename = TString::Format("%s/Opening.root", InputDir.Data());

  if (TFile::Open(filename)) {
    std::cout << "Sure you want to overwrite?\n";
    return 0;
  } else {
    std::cout << "No smearing matrix found - recomputing \n";
    auto outfile = new TFile(filename, "RECREATE");

//    TDirectoryFile *dir=(TDirectoryFile*)(file->FindObjectAny(Form("%sResults%d",trigger,suffix)));
//    std::cout<< "lol333"<<std::endl;

//    dir->GetObject(Form("%sResults%d",trigger,suffix), mainDir);

//    TDirectory *protonDir = file->GetDirectory(protonName);
//    auto protonList = (TList *) protonDir->Get(protonName);

// std::cout<<Form("%sResults%d",trigger,suffix)<<std::endl;
    TDirectory *mainDir = file->GetDirectory(Form("%sResults%d",trigger,suffix));
    mainDir->ls();
    TList* LIST = (TList *) mainDir->Get(Form("%sResults%d",trigger,suffix));
    //LIST->ls();
    TLorentzVector protonPart;
    TList* protonList = (TList *)LIST->FindObject("Proton");

    protonList = (TList *) protonList->FindObject("after");

    auto protonpT = (TH1F*) protonList->FindObject("pTDist_after");

    auto protoneta = (TH1F*) protonList->FindObject("EtaDist_after");

    auto protonphi = (TH1F*) protonList->FindObject("phiDist_after");

    auto protonMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    outfile->cd();

    protonpT->Write("protonpT");
    protoneta->Write("protonEta");
    protonphi->Write("protonPhi");


//    TLorentzVector phiPart;
//    TString phiName = trigger;
//    phiName += "v0Cuts";
//    phiName += suffix;
//    TDirectory *phiDir = file->GetDirectory(phiName);

//    auto phiList = (TList *) phiDir->Get(phiName);

    TLorentzVector truePhiPart;
    TList* phiList = (TList *)LIST->FindObject("Phi");
    phiList = (TList *) phiList->FindObject("v0Cuts");
    auto phiMassRec = (TH1F*) phiList->FindObject("InvMasswithCuts");
    phiMassRec->Write("PhiMassRec");
    phiList = (TList *) phiList->FindObject("after");
    auto phipT = (TH1F*) phiList->FindObject("pTDist_after");
    auto phieta = (TH1F*) phiList->FindObject("EtaDist_after");
    auto phiphi = (TH1F*) phiList->FindObject("PhiDist_after");
    auto phiMass = TDatabasePDG::Instance()->GetParticle(333)->Mass();
    float entries = phiMassRec->Integral(phiMassRec->FindBin(phiMass-0.008),
                                           phiMassRec->FindBin(phiMass+0.008));
    outfile->cd();
    phipT->Write("phipT");
    phieta->Write("phiEta");
    phiphi->Write("phiPhi");
    const float purity = 0.66;

    TLorentzVector phiPart;



//    TLorentzVector trueSigmaPart;
//    TString sigmaName = trigger;
//    sigmaName += "Sigma0Cuts";
//    sigmaName += suffix;
//    TDirectory *sigmaDir = file->GetDirectory(sigmaName);
//    auto sigmaList = (TList *) sigmaDir->Get(sigmaName);
//    auto sigmaMassRec = (TH1F*) sigmaList->FindObject("fHistInvMass");
//    auto sigmapT = (TH1F*) sigmaList->FindObject("fHistMassCutPt");
//    auto sigmaetaphi = (TH2F*) sigmaList->FindObject("fHistEtaPhi");
//    auto sigmaeta = (TH1F*) sigmaetaphi->ProjectionX();
//    auto sigmaphi = (TH1F*) sigmaetaphi->ProjectionY();
//    sigmaMassRec->Write("sigmaMassRec");
//    float entries = sigmaMassRec->Integral(sigmaMassRec->FindBin(1.2),
//                                           sigmaMassRec->FindBin(1.25));
//    auto sigmaMass = TDatabasePDG::Instance()->GetParticle(3212)->Mass();
//    const float purity = 0.346;

//    TLorentzVector phiPart;

    TLorentzVector kpPart;
    TList* kpList = (TList *)LIST->FindObject("Particle1");
    kpList = (TList *) kpList->FindObject("after");
    auto kppT = (TH1F*) kpList->FindObject("pTDist_after");
    auto kpeta = (TH1F*) kpList->FindObject("EtaDist_after");
    auto kpphi = (TH1F*) kpList->FindObject("phiDist_after");
    auto kpMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    outfile->cd();
    kppT->Write("kppT");
    kpeta->Write("kpEta");
    kpphi->Write("kpPhi");


    TLorentzVector kmPart;
    TList* kmList = (TList *)LIST->FindObject("Particle2");
    kmList = (TList *) kmList->FindObject("after");
    auto kmpT = (TH1F*) kmList->FindObject("pTDist_after");
    auto kmeta = (TH1F*) kmList->FindObject("EtaDist_after");
    auto kmphi = (TH1F*) kmList->FindObject("phiDist_after");
    auto kmMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    outfile->cd();
    kmpT->Write("kmpT");
    kmeta->Write("kmEta");
    kmphi->Write("kmPhi");


    auto histPhi = new TH1F("histPhi", "", 1000, 1, 2);
    auto histPhiPhi = new TH1F("histPhiPhi", "", 1000, 1, 2);
    auto histPhiPt = new TH1F("histPhiPt", "", 1000, 0, 10);
    auto opening =
        new TH2F(
            "openingAngle",
            "; opening angle in deg; Inv Mass (GeV/#it{c^{2}})",
            360, 0, 360, 200, 0, 20);

//    auto histSmear2 =
//        new TH2F(
//            "histSmear2Km",
//            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
//            500, 0, 1000, 500, 0, 1000);

    double counter = 0;
    int modulo = (nParticles > 10000) ? nParticles / 10000 : 1;
    int i = 0;
    double kstarpPhi, kstarpKp, kstarpKm, massCurrentPhi, kstarpPhitrue , openingangle;
    double nCount = nParticles;
//    TGenPhaseSpace eventPhi;
//    TLorentzVector phiPhiDecay;
//    double PhiDecay[2] = { phiMass, 0. };
    TStopwatch watch;
    while (i < nParticles) {
      protonPart.SetPtEtaPhiM(protonpT->GetRandom(), protoneta->GetRandom(),
                              protonphi->GetRandom(), protonMass);
      truePhiPart.SetPtEtaPhiM(phipT->GetRandom(), phieta->GetRandom(),
                              phiphi->GetRandom(), phiMass);
      kpPart.SetPtEtaPhiM(kppT->GetRandom(), kpeta->GetRandom(),
                              kpphi->GetRandom(), kpMass);
      kmPart.SetPtEtaPhiM(kmpT->GetRandom(), kmeta->GetRandom(),
                              kmphi->GetRandom(), kmMass);
      phiPart = kpPart + kmPart;
      massCurrentPhi = phiPart.M();
      histPhi->Fill(massCurrentPhi);
      histPhiPt->Fill(phiPart.Pt());

      openingangle=GetOpeningAngle(kpPart,kmPart,massCurrentPhi);


      opening->Fill(openingangle,massCurrentPhi);
//      histSmearFull->Fill(kstarpKp, kstarpPhi);
//      histSmear2Full->Fill(kstarpKm, kstarpPhi);


//      if (std::abs(massCurrentPhi - phiMass) < 0.008) {
//        histSmearFullPhiPeak->Fill(kstarpKp, kstarpPhi);
//        histSmearFullPhiPeak2->Fill(kstarpKm, kstarpPhi);

//        histSmearRatioPeak->Fill(kstarpKp, kstarpPhi);
//        histSmearTOTPeak->Fill(kstarpKp, kstarpPhi);


//      }

//      if ((0.987<massCurrentPhi) && (massCurrentPhi<1.011)) {
//        histSmearFullSBleft->Fill(kstarpKp, kstarpPhi);
//        histSmear2FullSBleft->Fill(kstarpKm, kstarpPhi);

//        histSmearRatioSBleft->Fill(kstarpKp, kstarpPhi);
//        histSmearTOTSBleft->Fill(kstarpKp, kstarpPhi);

//      }

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



//    histPhi->Scale(
//        entries
//            / histPhi->Integral(histPhi->FindBin(phiMass-0.008),
//                                  histPhi->FindBin(phiMass+0.008)));


    opening->Write();
    outfile->Close();
    return 0;
  }
}
