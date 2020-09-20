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
  auto filename = TString::Format("%s/SmearSideband.root", InputDir.Data());

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


    auto histPhi = new TH1F("histPhi", "", 2000, 1, 2);
    auto histPhiPhi = new TH1F("histPhiPhi", "", 2000, 1, 2);
    auto histPhiPt = new TH1F("histPhiPt", "", 2000, 0, 10);
    auto histSmear =
        new TH2F(
            "histSmearKp",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2 =
        new TH2F(
            "histSmear2Km",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);


    auto histSmearFull =
        new TH2F(
            "histSmearFull",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full =
        new TH2F(
            "histSmear2Full",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFullPhiPeak =
        new TH2F(
            "histSmearFullPhiPeakKp",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFullPhiPeak2 =
        new TH2F(
            "histSmearFullPhiPeakKm",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFullSBleft =
        new TH2F(
            "histSmearFullSBleftKp",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2FullSBleft =
        new TH2F(
            "histSmearFullSBleftKm",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFullSBright =
        new TH2F(
            "histSmearFullSBrightKp",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2FullSBright =
        new TH2F(
            "histSmearFullSBrightKm",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull12 =
        new TH2F(
            "histSmearFullKp12",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full12 =
        new TH2F(
            "histSmearFullKm12",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull23 =
        new TH2F(
            "histSmearFullKp23",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full23 =
        new TH2F(
            "histSmearFullKm23",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull34 =
        new TH2F(
            "histSmearFullKp34",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full34 =
        new TH2F(
            "histSmearFullKm34",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull45 =
        new TH2F(
            "histSmearFullKp45",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full45 =
        new TH2F(
            "histSmearFullKm45",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull56 =
        new TH2F(
            "histSmearFullKp56",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full56 =
        new TH2F(
            "histSmearFullKm56",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull67 =
        new TH2F(
            "histSmearFullKp67",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full67 =
        new TH2F(
            "histSmearFullKm67",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull78 =
        new TH2F(
            "histSmearFullKp78",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full78 =
        new TH2F(
            "histSmearFullKm78",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull89 =
        new TH2F(
            "histSmearFullKp89",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full89 =
        new TH2F(
            "histSmearFullKm89",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearFull90 =
        new TH2F(
            "histSmearFullKp90",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2Full90 =
        new TH2F(
            "histSmearFullKm90",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioPeak =
        new TH2F(
            "histSmearRatioPeak",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTPeak =
        new TH2F(
            "histSmearTOTPeak",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSBleft =
        new TH2F(
            "histSmearRatioSBleft",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSBleft =
        new TH2F(
            "histSmearTOTSBleft",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSBright =
        new TH2F(
            "histSmearRatioSBright",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSBright=
        new TH2F(
            "histSmearTOTSBright",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB12 =
        new TH2F(
            "histSmearRatioSB12",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB12 =
        new TH2F(
            "histSmearTOTSB12",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB23 =
        new TH2F(
            "histSmearRatioSB23",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB23 =
        new TH2F(
            "histSmearTOTSB23",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB34 =
        new TH2F(
            "histSmearRatioSB34",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB34 =
        new TH2F(
            "histSmearTOTSB34",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB45 =
        new TH2F(
            "histSmearRatioSB45",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB45 =
        new TH2F(
            "histSmearTOTSB45",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB56 =
        new TH2F(
            "histSmearRatioSB56",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB56 =
        new TH2F(
            "histSmearTOTSB56",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB67 =
        new TH2F(
            "histSmearRatioSB67",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB67 =
        new TH2F(
            "histSmearTOTSB67",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB78 =
        new TH2F(
            "histSmearRatioSB78",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB78 =
        new TH2F(
            "histSmearTOTSB78",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB89 =
        new TH2F(
            "histSmearRatioSB89",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB89 =
        new TH2F(
            "histSmearTOTSB89",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatioSB90 =
        new TH2F(
            "histSmearRatioSB90",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOTSB90 =
        new TH2F(
            "histSmearTOTSB90",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);


    auto histSmearsmall2762 =
        new TH2F(
            "histSmearsmall2762",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small2762 =
        new TH2F(
            "histSmear2small2762",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio2762 =
        new TH2F(
            "histSmearRatio2762",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT2762 =
        new TH2F(
            "histSmearTOT2762",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);


    auto histSmearsmall6210 =
        new TH2F(
            "histSmearsmall6210",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small6210 =
        new TH2F(
            "histSmear2small6210",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio6210 =
        new TH2F(
            "histSmearRatio6210",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT6210 =
        new TH2F(
            "histSmearTOT6210",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall1015 =
        new TH2F(
            "histSmearsmall1015",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small1015 =
        new TH2F(
            "histSmear2small1015",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio1015 =
        new TH2F(
            "histSmearRatio1015",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT1015 =
        new TH2F(
            "histSmearTOT1015",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall1520 =
        new TH2F(
            "histSmearsmall1520",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small1520 =
        new TH2F(
            "histSmear2small1520",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio1520 =
        new TH2F(
            "histSmearRatio1520",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT1520 =
        new TH2F(
            "histSmearTOT1520",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall2025 =
        new TH2F(
            "histSmearsmall2025",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small2025 =
        new TH2F(
            "histSmear2small2025",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio2025 =
        new TH2F(
            "histSmearRatio2025",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT2025 =
        new TH2F(
            "histSmearTOT2025",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall2530 =
        new TH2F(
            "histSmearsmall2530",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small2530 =
        new TH2F(
            "histSmear2small2530",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio2530 =
        new TH2F(
            "histSmearRatio2530",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT2530 =
        new TH2F(
            "histSmearTOT2530",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall3035 =
        new TH2F(
            "histSmearsmall3035",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small3035 =
        new TH2F(
            "histSmear2small3035",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio3035 =
        new TH2F(
            "histSmearRatio3035",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT3035 =
        new TH2F(
            "histSmearTOT3035",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall3540 =
        new TH2F(
            "histSmearsmall3540",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small3540 =
        new TH2F(
            "histSmear2small3540",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio3540 =
        new TH2F(
            "histSmearRatio3540",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT3540 =
        new TH2F(
            "histSmearTOT3540",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall4045 =
        new TH2F(
            "histSmearsmall4045",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small4045 =
        new TH2F(
            "histSmear2small4045",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio4045 =
        new TH2F(
            "histSmearRatio4045",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT4045 =
        new TH2F(
            "histSmearTOT4045",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall4550 =
        new TH2F(
            "histSmearsmall4550",
            "; #it{k}*_{p#minusK^{+} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small4550 =
        new TH2F(
            "histSmear2small4550",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio4550 =
        new TH2F(
            "histSmearRatio4550",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT4550 =
        new TH2F(
            "histSmearTOT4550",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall5055 =
        new TH2F(
            "histSmearsmall5055",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small5055 =
        new TH2F(
            "histSmear2small5055",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio5055 =
        new TH2F(
            "histSmearRatio5055",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT5055 =
        new TH2F(
            "histSmearTOT5055",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearsmall5560 =
        new TH2F(
            "histSmearsmall5560",
            "; #it{k}*_{p#minusK^{+}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmear2small5560 =
        new TH2F(
            "histSmear2small5560",
            "; #it{k}*_{p#minusK^{-}} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearRatio5560 =
        new TH2F(
            "histSmearRatio5560",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT5560 =
        new TH2F(
            "histSmearTOT5560",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT2223 =
        new TH2F(
            "histSmearTOT2223",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT2425 =
        new TH2F(
            "histSmearTOT2425",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT2627 =
        new TH2F(
            "histSmearTOT2627",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT2829 =
        new TH2F(
            "histSmearTOT2829",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT3031 =
        new TH2F(
            "histSmearTOT3031",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT3334 =
        new TH2F(
            "histSmearTOT3334",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);
    auto histSmearTOT3738 =
        new TH2F(
            "histSmearTOT3738",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);

    auto histSmearTOT4041 =
        new TH2F(
            "histSmearTOT4041",
            "; #it{k}*_{p#minusK} (MeV/#it{c}); #it{k}*_{p#minus(K^{+}K^{-})} (MeV/#it{c})",
            1000, 0, 2000, 1000, 0, 2000);







    auto relMom = new TNtuple("relMom", "relMom", "pPhi:pkp:pkm");

    double counter = 0;
    int modulo = (nParticles > 10000) ? nParticles / 10000 : 1;
    int i = 0;
    double kstarpPhi, kstarpKp, kstarpKm, massCurrentPhi, kstarpPhitrue ;
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

      kstarpPhi = GetkStar(protonPart, phiPart);
      kstarpKp = GetkStar(protonPart, kpPart);
      kstarpKm = GetkStar(protonPart, kmPart);
      kstarpPhitrue = GetkStar(protonPart, truePhiPart);


      histSmearFull->Fill(kstarpKp, kstarpPhi);
      histSmear2Full->Fill(kstarpKm, kstarpPhi);


      if (std::abs(massCurrentPhi - phiMass) < 0.008) {
        histSmearFullPhiPeak->Fill(kstarpKp, kstarpPhi);
        histSmearFullPhiPeak2->Fill(kstarpKm, kstarpPhi);

        histSmearRatioPeak->Fill(kstarpKp, kstarpPhi);
        histSmearTOTPeak->Fill(kstarpKp, kstarpPhi);


      }

      if ((0.987<massCurrentPhi) && (massCurrentPhi<1.011)) {
        histSmearFullSBleft->Fill(kstarpKp, kstarpPhi);
        histSmear2FullSBleft->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSBleft->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSBleft->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.027<massCurrentPhi) && (massCurrentPhi< 1.1)) {
        histSmearFullSBright->Fill(kstarpKp, kstarpPhi);
        histSmear2FullSBright->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSBright->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSBright->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.1<massCurrentPhi) && (massCurrentPhi< 1.2)) {
        histSmearFull12->Fill(kstarpKp, kstarpPhi);
        histSmear2Full12->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB12->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB12->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.2<massCurrentPhi) && (massCurrentPhi< 1.3)) {
        histSmearFull23->Fill(kstarpKp, kstarpPhi);
        histSmear2Full23->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB23->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB23->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.3<massCurrentPhi) && (massCurrentPhi< 1.4)) {
        histSmearFull34->Fill(kstarpKp, kstarpPhi);
        histSmear2Full34->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB34->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB34->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.4<massCurrentPhi) && (massCurrentPhi< 1.5)) {
        histSmearFull45->Fill(kstarpKp, kstarpPhi);
        histSmear2Full45->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB45->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB45->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.5<massCurrentPhi) && (massCurrentPhi< 1.6)) {
        histSmearFull56->Fill(kstarpKp, kstarpPhi);
        histSmear2Full56->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB56->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB56->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.6<massCurrentPhi) && (massCurrentPhi< 1.7)) {
        histSmearFull67->Fill(kstarpKp, kstarpPhi);
        histSmear2Full67->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB67->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB67->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.7<massCurrentPhi) && (massCurrentPhi< 1.8)) {
        histSmearFull78->Fill(kstarpKp, kstarpPhi);
        histSmear2Full78->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB78->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB78->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.8<massCurrentPhi) && (massCurrentPhi< 1.9)) {
        histSmearFull89->Fill(kstarpKp, kstarpPhi);
        histSmear2Full89->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB89->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB89->Fill(kstarpKp, kstarpPhi);

      }

      if ((1.9<massCurrentPhi) && (massCurrentPhi< 2.0)) {
        histSmearFull90->Fill(kstarpKp, kstarpPhi);
        histSmear2Full90->Fill(kstarpKm, kstarpPhi);

        histSmearRatioSB90->Fill(kstarpKp, kstarpPhi);
        histSmearTOTSB90->Fill(kstarpKp, kstarpPhi);

      }



      if ((1.027<massCurrentPhi) && (massCurrentPhi< 1.062)) {

      histSmearsmall2762->Fill(kstarpKp, kstarpPhi);
      histSmear2small2762->Fill(kstarpKm, kstarpPhi);

      histSmearRatio2762->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2762->Fill(kstarpKp, kstarpPhi);

}

      if ((1.062<massCurrentPhi) && (massCurrentPhi< 1.1)) {

      histSmearsmall6210->Fill(kstarpKp, kstarpPhi);
      histSmear2small6210->Fill(kstarpKm, kstarpPhi);

      histSmearRatio6210->Fill(kstarpKp, kstarpPhi);
      histSmearTOT6210->Fill(kstarpKp, kstarpPhi);

}
      if ((1.1<massCurrentPhi) && (massCurrentPhi< 1.15)) {

      histSmearsmall1015->Fill(kstarpKp, kstarpPhi);
      histSmear2small1015->Fill(kstarpKm, kstarpPhi);

      histSmearRatio1015->Fill(kstarpKp, kstarpPhi);
      histSmearTOT1015->Fill(kstarpKp, kstarpPhi);

}
      if ((1.15<massCurrentPhi) && (massCurrentPhi< 1.2)) {

      histSmearsmall1520->Fill(kstarpKp, kstarpPhi);
      histSmear2small1520->Fill(kstarpKm, kstarpPhi);

      histSmearRatio1520->Fill(kstarpKp, kstarpPhi);
      histSmearTOT1520->Fill(kstarpKp, kstarpPhi);

}
      if ((1.2<massCurrentPhi) && (massCurrentPhi< 1.25)) {

      histSmearsmall2025->Fill(kstarpKp, kstarpPhi);
      histSmear2small2025->Fill(kstarpKm, kstarpPhi);

      histSmearRatio2025->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2025->Fill(kstarpKp, kstarpPhi);

}

      if ((1.25<massCurrentPhi) && (massCurrentPhi< 1.3)) {

      histSmearsmall2530->Fill(kstarpKp, kstarpPhi);
      histSmear2small2530->Fill(kstarpKm, kstarpPhi);

      histSmearRatio2530->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2530->Fill(kstarpKp, kstarpPhi);


}
      if ((1.3<massCurrentPhi) && (massCurrentPhi< 1.35)) {

      histSmearsmall3035->Fill(kstarpKp, kstarpPhi);
      histSmear2small3035->Fill(kstarpKm, kstarpPhi);

      histSmearRatio3035->Fill(kstarpKp, kstarpPhi);
      histSmearTOT3035->Fill(kstarpKp, kstarpPhi);

}
      if ((1.35<massCurrentPhi) && (massCurrentPhi< 1.4)) {

      histSmearsmall3540->Fill(kstarpKp, kstarpPhi);
      histSmear2small3540->Fill(kstarpKm, kstarpPhi);

      histSmearRatio3540->Fill(kstarpKp, kstarpPhi);
      histSmearTOT3540->Fill(kstarpKp, kstarpPhi);


}
      if ((1.4<massCurrentPhi) && (massCurrentPhi< 1.45)) {

      histSmearsmall4045->Fill(kstarpKp, kstarpPhi);
      histSmear2small4045->Fill(kstarpKm, kstarpPhi);

      histSmearRatio4045->Fill(kstarpKp, kstarpPhi);
      histSmearTOT4045->Fill(kstarpKp, kstarpPhi);


}
      if ((1.45<massCurrentPhi) && (massCurrentPhi< 1.5)) {

      histSmearsmall4550->Fill(kstarpKp, kstarpPhi);
      histSmear2small4550->Fill(kstarpKm, kstarpPhi);

      histSmearRatio4550->Fill(kstarpKp, kstarpPhi);
      histSmearTOT4550->Fill(kstarpKp, kstarpPhi);

}

      if ((1.5<massCurrentPhi) && (massCurrentPhi< 1.55)) {

      histSmearsmall5055->Fill(kstarpKp, kstarpPhi);
      histSmear2small5055->Fill(kstarpKm, kstarpPhi);

      histSmearRatio5055->Fill(kstarpKp, kstarpPhi);
      histSmearTOT5055->Fill(kstarpKp, kstarpPhi);

}
      if ((1.55<massCurrentPhi) && (massCurrentPhi< 1.6)) {

      histSmearsmall5560->Fill(kstarpKp, kstarpPhi);
      histSmear2small5560->Fill(kstarpKm, kstarpPhi);

      histSmearRatio5560->Fill(kstarpKp, kstarpPhi);
      histSmearTOT5560->Fill(kstarpKp, kstarpPhi);


}

      if ((2.2<massCurrentPhi) && (massCurrentPhi< 2.3)) {

      histSmearTOT2223->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2223->Fill(kstarpKm, kstarpPhi);

}
      if ((2.4<massCurrentPhi) && (massCurrentPhi< 2.5)) {

      histSmearTOT2425->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2425->Fill(kstarpKm, kstarpPhi);

}
      if ((2.6<massCurrentPhi) && (massCurrentPhi< 2.7)) {

      histSmearTOT2627->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2627->Fill(kstarpKm, kstarpPhi);

}

      if ((2.8<massCurrentPhi) && (massCurrentPhi< 2.9)) {

      histSmearTOT2829->Fill(kstarpKp, kstarpPhi);
      histSmearTOT2829->Fill(kstarpKm, kstarpPhi);

}

      if ((3.0<massCurrentPhi) && (massCurrentPhi< 3.1)) {

      histSmearTOT3031->Fill(kstarpKp, kstarpPhi);
      histSmearTOT3031->Fill(kstarpKm, kstarpPhi);

}

      if ((3.3<massCurrentPhi) && (massCurrentPhi< 3.4)) {

      histSmearTOT3334->Fill(kstarpKp, kstarpPhi);
      histSmearTOT3334->Fill(kstarpKm, kstarpPhi);

}

      if ((3.7<massCurrentPhi) && (massCurrentPhi< 3.8)) {

      histSmearTOT3738->Fill(kstarpKp, kstarpPhi);
      histSmearTOT3738->Fill(kstarpKm, kstarpPhi);

}


      if ((4.0<massCurrentPhi) && (massCurrentPhi< 4.1)) {

      histSmearTOT4041->Fill(kstarpKp, kstarpPhi);
      histSmearTOT4041->Fill(kstarpKm, kstarpPhi);

}


//      if (std::abs(massCurrentPhi - phiMass) > 0.008) {
//        continue;
//      }

//      // let's add a true Phi
//      histPhiPhi->Fill(massCurrentPhi);
//      if (gRandom->Uniform() < purity) {
//        truePhiPart.SetPtEtaPhiM(phipT->GetRandom(), phieta->GetRandom(),
//                                   phiphi->GetRandom(),
//                                   gRandom->Gaus(phiMass, 0.0015));
//        histPhiPhi->Fill(truePhiPart.M());

//        eventPhi.SetDecay(truePhiPart, 2, PhiDecay);
//        eventPhi.Generate();
//        phiPhiDecay = *(eventPhi.GetDecay(0));
//        kstarpPhiPhi = GetkStar(protonPart, truePhiPart);
//        kstarpphiPhi = GetkStar(protonPart, phiPhiDecay);
//        histSmearPhi->Fill(kstarpphiPhi, kstarpPhiPhi);
//      }

      if (i % 1000 == 0)
        relMom->Fill(kstarpPhi, kstarpKp, kstarpKm);
     // histSmear->Fill(kstarpKp, kstarpPhi);
     // histSmear2->Fill(kstarpKm, kstarpPhi);
     // histSmearPhi->Fill(kstarpphi, kstarpPhi);
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

    histSmearRatioPeak->Divide(histSmearFullPhiPeak2);
    histSmearTOTPeak->Add(histSmearFullPhiPeak2);

    histSmearRatioSBleft->Divide(histSmear2FullSBleft);
    histSmearTOTSBleft->Add(histSmear2FullSBleft);

    histSmearRatioSBright->Divide(histSmear2FullSBright);
    histSmearTOTSBright->Add(histSmear2FullSBright);

    histSmearRatioSB12->Divide(histSmear2Full12);
    histSmearTOTSB12->Add(histSmear2Full12);

    histSmearRatioSB23->Divide(histSmear2Full23);
    histSmearTOTSB23->Add(histSmear2Full23);

    histSmearRatioSB34->Divide(histSmear2Full34);
    histSmearTOTSB34->Add(histSmear2Full34);

    histSmearRatioSB45->Divide(histSmear2Full45);
    histSmearTOTSB45->Add(histSmear2Full45);

    histSmearRatioSB56->Divide(histSmear2Full56);
    histSmearTOTSB56->Add(histSmear2Full56);

    histSmearRatioSB67->Divide(histSmear2Full67);
    histSmearTOTSB67->Add(histSmear2Full67);

    histSmearRatioSB78->Divide(histSmear2Full78);
    histSmearTOTSB78->Add(histSmear2Full78);

    histSmearRatioSB89->Divide(histSmear2Full89);
    histSmearTOTSB89->Add(histSmear2Full89);

    histSmearRatioSB90->Divide(histSmear2Full90);
    histSmearTOTSB90->Add(histSmear2Full90);




    histSmearRatio2762->Divide(histSmear2small2762);
    histSmearTOT2762->Add(histSmear2small2762);

    histSmearRatio6210->Divide(histSmear2small6210);
    histSmearTOT6210->Add(histSmear2small6210);

    histSmearRatio1015->Divide(histSmear2small1015);
    histSmearTOT1015->Add(histSmear2small1015);

    histSmearRatio1520->Divide(histSmear2small1520);
    histSmearTOT1520->Add(histSmear2small1520);

    histSmearRatio2025->Divide(histSmear2small2025);
    histSmearTOT2025->Add(histSmear2small2025);

    histSmearRatio2530->Divide(histSmear2small2530);
    histSmearTOT2530->Add(histSmear2small2530);

    histSmearRatio3035->Divide(histSmear2small3035);
    histSmearTOT3035->Add(histSmear2small3035);

    histSmearRatio3540->Divide(histSmear2small3540);
    histSmearTOT3540->Add(histSmear2small3540);

    histSmearRatio4045->Divide(histSmear2small4045);
    histSmearTOT4045->Add(histSmear2small4045);

    histSmearRatio4550->Divide(histSmear2small4550);
    histSmearTOT4550->Add(histSmear2small4550);

    histSmearRatio5055->Divide(histSmear2small5055);
    histSmearTOT5055->Add(histSmear2small5055);

    histSmearRatio5560->Divide(histSmear2small5560);
    histSmearTOT5560->Add(histSmear2small5560);


    histPhi->Scale(
        entries
            / histPhi->Integral(histPhi->FindBin(phiMass-0.008),
                                  histPhi->FindBin(phiMass+0.008)));


    histPhi->Write();
//    histPhiPhi->Write();
    histPhiPt->Write();
//    histSmear->Write();
//    histSmear2->Write();
    histSmearFullPhiPeak->Write();
    histSmearFullPhiPeak2->Write();
    histSmearFullSBright->Write();
    histSmear2FullSBright->Write();
    histSmearFull->Write();
    histSmear2Full->Write();
    histSmearFullSBleft->Write();
    histSmear2FullSBleft->Write();
    histSmearFull12->Write();
    histSmear2Full12->Write();
    histSmearFull23->Write();
    histSmear2Full23->Write();
    histSmearFull34->Write();
    histSmear2Full34->Write();
    histSmearFull45->Write();
    histSmear2Full45->Write();
    histSmearFull56->Write();
    histSmear2Full56->Write();
    histSmearFull67->Write();
    histSmear2Full67->Write();
    histSmearFull78->Write();
    histSmear2Full78->Write();
    histSmearFull89->Write();
    histSmear2Full89->Write();
    histSmearFull90->Write();
    histSmear2Full90->Write();

    histSmearRatioPeak->Write();
    histSmearTOTPeak->Write();
    histSmearRatioSBleft->Write();
    histSmearTOTSBleft->Write();
    histSmearRatioSBright->Write();
    histSmearTOTSBright->Write();
    histSmearRatioSB12->Write();
    histSmearTOTSB12->Write();
    histSmearRatioSB23->Write();
    histSmearTOTSB23->Write();
    histSmearRatioSB34->Write();
    histSmearTOTSB34->Write();
    histSmearRatioSB45->Write();
    histSmearTOTSB45->Write();
    histSmearRatioSB56->Write();
    histSmearTOTSB56->Write();
    histSmearRatioSB67->Write();
    histSmearTOTSB67->Write();
    histSmearRatioSB78->Write();
    histSmearTOTSB78->Write();
    histSmearRatioSB89->Write();
    histSmearTOTSB89->Write();
    histSmearRatioSB90->Write();
    histSmearTOTSB90->Write();

    histSmearsmall2762->Write();
    histSmear2small2762->Write();
    histSmearsmall6210->Write();
    histSmear2small6210->Write();
    histSmearsmall1015->Write();
    histSmear2small1015->Write();
    histSmearsmall1520->Write();
    histSmear2small1520->Write();
    histSmearsmall2025->Write();
    histSmear2small2025->Write();
    histSmearsmall2530->Write();
    histSmear2small2530->Write();
    histSmearsmall3035->Write();
    histSmear2small3035->Write();
    histSmearsmall3540->Write();
    histSmear2small3540->Write();
    histSmearsmall4045->Write();
    histSmear2small4045->Write();
    histSmearsmall4550->Write();
    histSmear2small4550->Write();
    histSmearsmall5055->Write();
    histSmear2small5055->Write();
    histSmearsmall5560->Write();
    histSmear2small5560->Write();

    histSmearRatio2762->Write();
    histSmearTOT2762->Write();
    histSmearRatio6210->Write();
    histSmearTOT6210->Write();
    histSmearRatio1015->Write();
    histSmearTOT1015->Write();
    histSmearRatio1520->Write();
    histSmearTOT1520->Write();
    histSmearRatio2025->Write();
    histSmearTOT2025->Write();
    histSmearRatio2530->Write();
    histSmearTOT2530->Write();
    histSmearRatio3035->Write();
    histSmearTOT3035->Write();
    histSmearRatio3540->Write();
    histSmearTOT3540->Write();
    histSmearRatio4045->Write();
    histSmearTOT4045->Write();
    histSmearRatio4550->Write();
    histSmearTOT4550->Write();
    histSmearRatio5055->Write();
    histSmearTOT5055->Write();
    histSmearRatio5560->Write();
    histSmearTOT5560->Write();


    histSmearTOT2223->Write();
    histSmearTOT2425->Write();
    histSmearTOT2627->Write();
    histSmearTOT2829->Write();

    histSmearTOT3031->Write();
    histSmearTOT3334->Write();
    histSmearTOT3738->Write();

    histSmearTOT4041->Write();


    relMom->Write();
    outfile->Close();
    return 0;
  }
}
