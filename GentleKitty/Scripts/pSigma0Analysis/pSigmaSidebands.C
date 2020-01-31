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
void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

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
TH2D* ComputeSmearingMatrix(TString &InputDir, TString &trigger,
                            TString &suffix, int nParticles) {
  auto file = TFile::Open(Form("%s/AnalysisResults.root", InputDir.Data()));
  auto filename = TString::Format("%s/SmearSideband.root", InputDir.Data());

  if (TFile::Open(filename)) {
    std::cout << "Using smearing matrix from existing file\n";
    auto outfile = TFile::Open(filename);
    auto histSmear = (TH2D*) outfile->Get("histSmearSigma");
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
      histSigmaSigma->Fill(massCurrentSigma);

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

      if (std::abs(massCurrentSigma - sigmaMass) > 0.005) {
        continue;
      }

      // let's add a true sigma
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

    histSigmaSigma->Scale(
        entries
            / histSigmaSigma->Integral(histSigmaSigma->FindBin(1.2),
                                       histSigmaSigma->FindBin(1.25)));
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
    return histSmearSigma;
  }
}

/// =====================================================================================
TGraph *GetSmearedCF(TGraph* CF, TH2* matrix) {
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

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

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
  auto SBmerge = side->GetSideBandGraph(5);

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

  // pp radius systematic variations
  const double ppRadius = 1.154;

  CATS AB_pL;
  tidy->GetCatsProtonLambda(&AB_pL, 41, -9.99, 400, TidyCats::sGaussian,
                            TidyCats::pNLOWF);
  AB_pL.KillTheCat();
  DLM_Ck* Ck_pL = new DLM_Ck(1, 0, AB_pL);
  Ck_pL->SetSourcePar(0, ppRadius);
  Ck_pL->Update();
  auto grSidebandRaw = new TGraph();
  FillWaveGraph(AB_pL, grSidebandRaw);

  const double protonPurity = 0.9943;
  const double protonPrimary = 0.823;
  const double protonLambda = 0.125;
  const double protonSecondary = protonLambda / (1. - protonPrimary);
  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  const double lambdaPurity = 0.946;
  const double lambdaPrimary = 0.6;
  const Particle lambda(lambdaPurity, lambdaPrimary,
                        { { 1.f - lambdaPrimary } });

  const CATSLambdaParam lambdaParam(proton, lambda);
  const float lambdaplambda = lambdaParam.GetLambdaParam(
      CATSLambdaParam::Primary);

  std::cout << "considering a primary fraction for the p-L CF of "
            << lambdaplambda << "\n";

  // apply the smearing with the photon, and the momentum smearing
  auto grSideband = GetSmearedCF(GetSmearedCF(grSidebandRaw, histSmear),
                                 CATSinput->GetSigmaFile(1));
  auto grSidebandSmeared = new TGraph();

  double x, y;
  for (int i = 0; i < 1000; ++i) {
    grSideband->GetPoint(i, x, y);
    grSidebandSmeared->SetPoint(i, x, 1 + (y - 1) * lambdaplambda);
  }

  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            350);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  DreamPlot::SetStyleHisto(histSmear);
  auto c = new TCanvas("c", "c", 650, 550);
  histSmear->GetZaxis()->SetTitleFont(43);
  histSmear->GetZaxis()->SetTitleSize(28);
  histSmear->GetZaxis()->SetLabelFont(43);
  histSmear->GetZaxis()->SetLabelSize(28);
  histSmear->GetXaxis()->SetNdivisions(505);
  histSmear->GetYaxis()->SetNdivisions(505);
  histSmear->Draw("col");

  histSmear->SetTitle(
      "; #it{k}*_{p#minus#Lambda} (MeV/#it{c}); #it{k}*_{p#minus#Sigma^{0}} (MeV/#it{c})");
  c->Print("plambdaSmearingMatrix.pdf");

  int lineWidth = 2;
  DreamPlot::SetStyleGraph(grSidebandRaw, 20, kGreen + 3, 0.7);
  grSidebandRaw->SetLineWidth(lineWidth);
  grSidebandRaw->SetLineStyle(3);
  DreamPlot::SetStyleGraph(grSideband, 20, kGreen + 3, 0.7);
  grSideband->SetLineWidth(lineWidth);
  grSideband->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSidebandSmeared, 20, kGreen + 3, 0.7);
  grSidebandSmeared->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(SBmerge, kOpenCircle, kBlue + 3);

  auto gc3 = new TCanvas("CFpL_lambda", "CFpL_lambda", 0, 0, 650, 550);
  gc3->SetRightMargin(right);
  gc3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 3.5);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grSidebandRaw->Draw("L3same");
  grSideband->Draw("L3same");
  grSidebandSmeared->Draw("L3same");
  SBmerge->Draw("pezsame");
  auto gleg5 = new TLegend(0.4, 0.6, 0.4 + 0.45, 0.9);
  gleg5->SetBorderSize(0);
  gleg5->SetTextFont(42);
  gleg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  gleg5->AddEntry(grSidebandRaw,
                  "Genuine p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)", "l");
  gleg5->AddEntry(
      grSideband,
      "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0}",
      "l");
  gleg5->AddEntry(
      grSidebandSmeared,
      Form(
          "p#minus#kern[-0.8]{ }#Lambda #rightarrow p#minus#kern[-0.95]{ }#Sigma^{0} (#lambda = %.3f)",
          lambdaplambda),
      "l");
  gleg5->AddEntry(
      SBmerge,
      "p#minus#kern[-0.65]{ }(#Lambda#gamma) #oplus #bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma)",
      "pez");
  gleg5->Draw("same");
  gc3->Print("sidebandsmearing.pdf");

  outfile->cd();
  histSmear->Write();
  grSidebandRaw->Write("p-Lambda raw");
  grSideband->Write("p-Lambda smeared");
  grSidebandSmeared->Write("p-Lambda smeared with Photon");
  SBmerge->Write();
  CATSinput->GetSigmaFile(1)->Write();
  outfile->Close();
}
