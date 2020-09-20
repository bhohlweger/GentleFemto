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
TGraph *GetSmearedCF(TGraph* CF, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  //std::cout<<"nbinsor:"<<nbins_original<<std::endl;
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();
  //std::cout<<"nbins_transformed:"<<nbins_transformed<<std::endl;

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
     // std::cout<<"momentum_original"<<momentum_original<<std::endl;

      //std::cout<<"jj"<<CF->Eval(momentum_original)<<std::endl;

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

TGraph *GetSmearedCF2(TGraph* CF, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetYaxis()->GetNbins();
  //std::cout<<"nbinsor:"<<nbins_original<<std::endl;
  const Int_t nbins_transformed = matrix->GetXaxis()->GetNbins();
  //std::cout<<"nbins_transformed:"<<nbins_transformed<<std::endl;

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetYaxis()->GetBinCenter(momOri + 1);
     // std::cout<<"momentum_original"<<momentum_original<<std::endl;

      //std::cout<<"jj"<<CF->Eval(momentum_original)<<std::endl;

      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetXaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF);
  }
  return smearedCF;
}

//TGraph *GetSmearedCF(TGraph CF, TH2* matrix) {
//  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
//    std::cout<<"hoho"<<std::endl;

//  const int nbins_original = matrix->GetXaxis()->GetNbins();
//  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();
//  std::cout<<"aha"<<std::endl;

//  TGraph *smearedCF = new TGraph();

//  int countPoint = 0;
//  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
//    Double_t matrixvalues_sum = 0.;
//    Double_t weighted_matrixvalues_sum = 0.;

//    for (int momOri = 0; momOri < nbins_original; momOri++) {
//      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
//      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
//      std::cout<<"alal"<<std::endl;

//      weighted_matrixvalues_sum += CF.Eval(momentum_original)
//          * matrix->GetBinContent(momOri + 1, momTrans + 1);
//      std::cout<<"olol"<<std::endl;

//    }
//    std::cout<<"kkk"<<std::endl;

//    Double_t transformed_CF = 0.;
//    if (matrixvalues_sum != 0.)
//      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

//    smearedCF->SetPoint(countPoint++,
//                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
//                        transformed_CF);
//  }
//  return smearedCF;
//}

TH1F* CalcCF( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/FemtoPhiHM_ROT/240340/CFOutput_pPhi_%s_%s.root",prefix, addon);
  auto file = TFile::Open(name);
 // file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}

TH1F* CalcCFMC( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/PhiTRUTH_2/phitruth/CFOutput_pPhi_%s_%s.root",prefix, addon);
  auto file = TFile::Open(name);
  //file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}

/// =====================================================================================
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  const char* trigger = argv[2];
  int suffix= atoi(argv[3]);

  auto filename = TString::Format("%s/SherlockSideband.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");

//  auto file = TFile::Open(Form("%s/500-800/CFOutput_%s_%d", trigger,suffix));
//  TH1F* CF = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");

//  CF->Fit("pol3");


//  TRandom3 rangen(0);
//  TidyCats* tidy = new TidyCats();  // for some reason we need this for the thing to compile

//  auto side = new SidebandSigma();
//  side->SetRebin(10);
//  side->SetSideBandFile(InputDir.Data(), trigger.Data(), suffix.Data());
//  const double sidebandNormDown = 250;
//  const double sidebandNormUp = 400;
//  side->SetNormalizationRange(sidebandNormDown, sidebandNormUp);

//  side->SideBandCFs();
//  auto SBmerge = side->GetSideBandGraph(5);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Get the smearing matrix



  auto filenameSmear = TString::Format("%s/SmearSideband.root", InputDir.Data());
  auto infile = TFile::Open(filenameSmear);
  if (!infile) {
    std::cout
        << "No smearing matrix found - start the sideband computation task!\n";
    return 0;
  }

//  auto histSmear = (TH2D*) infile->Get("histSmear");
  auto histSmearFullSBrightKp = (TH2F*) infile->Get("histSmearFullSBrightKp");
  auto histSmearFullSBrightKm = (TH2F*) infile->Get("histSmearFullSBrightKm");
  auto histSmearFullSBleftKp = (TH2F*) infile->Get("histSmearFullSBleftKp");
  auto histSmearFullSBleftKm = (TH2F*) infile->Get("histSmearFullSBleftKm");
  auto histSmearFullPhiPeakKp = (TH2F*) infile->Get("histSmearFullPhiPeakKp");
  auto histSmearFullPhiPeakKm = (TH2F*) infile->Get("histSmearFullPhiPeakKm");

  auto histSmearFullKp12 = (TH2F*) infile->Get("histSmearFullKp12");
  auto histSmearFullKm12 = (TH2F*) infile->Get("histSmearFullKm12");
  auto histSmearFullKp23 = (TH2F*) infile->Get("histSmearFullKp23");
  auto histSmearFullKm23 = (TH2F*) infile->Get("histSmearFullKm23");
  auto histSmearFullKp34 = (TH2F*) infile->Get("histSmearFullKp34");
  auto histSmearFullKm34 = (TH2F*) infile->Get("histSmearFullKm34");
  auto histSmearFullKp45 = (TH2F*) infile->Get("histSmearFullKp45");
  auto histSmearFullKm45 = (TH2F*) infile->Get("histSmearFullKm45");
  auto histSmearFullKp56 = (TH2F*) infile->Get("histSmearFullKp56");
  auto histSmearFullKm56 = (TH2F*) infile->Get("histSmearFullKm56");
  auto histSmearFullKp67 = (TH2F*) infile->Get("histSmearFullKp67");
  auto histSmearFullKm67 = (TH2F*) infile->Get("histSmearFullKm67");
  auto histSmearFullKp78 = (TH2F*) infile->Get("histSmearFullKp78");
  auto histSmearFullKm78 = (TH2F*) infile->Get("histSmearFullKm78");
  auto histSmearFullKp89 = (TH2F*) infile->Get("histSmearFullKp89");
  auto histSmearFullKm89 = (TH2F*) infile->Get("histSmearFullKm89");
  auto histSmearFullKp90 = (TH2F*) infile->Get("histSmearFullKp90");
  auto histSmearFullKm90 = (TH2F*) infile->Get("histSmearFullKm90");

  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 5,
                            500);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  DreamPlot::SetStyleHisto(histSmearFullKp12);
  DreamPlot::SetStyleHisto(histSmearFullKm12);
  DreamPlot::SetStyleHisto(histSmearFullKp23);
  DreamPlot::SetStyleHisto(histSmearFullKm23);
  DreamPlot::SetStyleHisto(histSmearFullKp34);
  DreamPlot::SetStyleHisto(histSmearFullKm34);
  DreamPlot::SetStyleHisto(histSmearFullKp45);
  DreamPlot::SetStyleHisto(histSmearFullKm45);
  DreamPlot::SetStyleHisto(histSmearFullKp56);
  DreamPlot::SetStyleHisto(histSmearFullKm56);
  DreamPlot::SetStyleHisto(histSmearFullKp67);
  DreamPlot::SetStyleHisto(histSmearFullKm67);
  DreamPlot::SetStyleHisto(histSmearFullKp78);
  DreamPlot::SetStyleHisto(histSmearFullKm67);
  DreamPlot::SetStyleHisto(histSmearFullKp78);
  DreamPlot::SetStyleHisto(histSmearFullKm78);
  DreamPlot::SetStyleHisto(histSmearFullKp89);
  DreamPlot::SetStyleHisto(histSmearFullKm89);
  DreamPlot::SetStyleHisto(histSmearFullKp90);
  DreamPlot::SetStyleHisto(histSmearFullKm90);





  auto c = new TCanvas("c", "c", 650, 550);
  histSmearFullPhiPeakKp->Draw("col");
  histSmearFullPhiPeakKp->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullPhiPeakKp->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullPhiPeakKp->GetXaxis()->SetNdivisions(505);
  histSmearFullPhiPeakKp->GetYaxis()->SetNdivisions(505);
  histSmearFullPhiPeakKp->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
  c->Print("pPhiSmearingMatrixPeakKp.pdf");
  c->Print("pPhiSmearingMatrixPeakKp.png");
 delete c;

  auto d = new TCanvas("d", "d", 650, 550);
  histSmearFullPhiPeakKm->Draw("col");
  histSmearFullPhiPeakKm->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullPhiPeakKm->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullPhiPeakKm->GetXaxis()->SetNdivisions(505);
  histSmearFullPhiPeakKm->GetYaxis()->SetNdivisions(505);
  histSmearFullPhiPeakKm->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
  d->Print("pPhiSmearingMatrixPeakKm.pdf");
  d->Print("pPhiSmearingMatrixPeakKm.png");
delete d;

  auto e = new TCanvas("e", "e", 650, 550);
  histSmearFullSBleftKp->Draw("col");
  histSmearFullSBleftKp->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullSBleftKp->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullSBleftKp->GetXaxis()->SetNdivisions(505);
  histSmearFullSBleftKp->GetYaxis()->SetNdivisions(505);
  histSmearFullSBleftKp->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  e->Print("pPhiSmearingMatrixSBleftKp.pdf");
  e->Print("pPhiSmearingMatrixSBleftKp.png");
delete e;

  auto f = new TCanvas("f", "f", 650, 550);
  histSmearFullSBleftKm->Draw("col");
  histSmearFullSBleftKm->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullSBleftKm->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullSBleftKm->GetXaxis()->SetNdivisions(505);
  histSmearFullSBleftKm->GetYaxis()->SetNdivisions(505);
  histSmearFullSBleftKm->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  f->Print("pPhiSmearingMatrixSBleftKm.pdf");
  f->Print("pPhiSmearingMatrixSBleftKm.png");
delete f;

  auto g = new TCanvas("g", "g", 650, 550);
  histSmearFullSBrightKp->Draw("col");
  histSmearFullSBrightKp->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullSBrightKp->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullSBrightKp->GetXaxis()->SetNdivisions(505);
  histSmearFullSBrightKp->GetYaxis()->SetNdivisions(505);
  histSmearFullSBrightKp->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  g->Print("pPhiSmearingMatrixSBrightKp.pdf");
  g->Print("pPhiSmearingMatrixSBrightKp.png");

  delete g;

  auto h = new TCanvas("h", "h", 650, 550);
  histSmearFullSBrightKm->Draw("col");
  histSmearFullSBrightKm->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullSBrightKm->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullSBrightKm->GetXaxis()->SetNdivisions(505);
  histSmearFullSBrightKm->GetYaxis()->SetNdivisions(505);
  histSmearFullSBrightKm->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  h->Print("pPhiSmearingMatrixSBrightKm.pdf");
  h->Print("pPhiSmearingMatrixSBrightKm.png");
delete h;


  auto hp12 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp12->Draw("col");
  histSmearFullKp12->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp12->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp12->GetXaxis()->SetNdivisions(505);
  histSmearFullKp12->GetYaxis()->SetNdivisions(505);
  histSmearFullKp12->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp12->Print("pPhiSmearingMatrixSB12Kp.pdf");
  hp12->Print("pPhiSmearingMatrixSB12Kp.png");
delete hp12;

  auto hm12 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm12->Draw("col");
  histSmearFullKm12->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm12->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm12->GetXaxis()->SetNdivisions(505);
  histSmearFullKm12->GetYaxis()->SetNdivisions(505);
  histSmearFullKm12->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm12->Print("pPhiSmearingMatrixSB12Km.pdf");
  hm12->Print("pPhiSmearingMatrixSB12Km.png");
delete hm12;


  auto hp23 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp23->Draw("col");
  histSmearFullKp23->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp23->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp23->GetXaxis()->SetNdivisions(505);
  histSmearFullKp23->GetYaxis()->SetNdivisions(505);
  histSmearFullKp23->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp23->Print("pPhiSmearingMatrixSB23Kp.pdf");
  hp23->Print("pPhiSmearingMatrixSB23Kp.png");
delete hp23;

  auto hm23 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm23->Draw("col");
  histSmearFullKm23->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm23->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm23->GetXaxis()->SetNdivisions(505);
  histSmearFullKm23->GetYaxis()->SetNdivisions(505);
  histSmearFullKm23->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm23->Print("pPhiSmearingMatrixSB23Km.pdf");
  hm23->Print("pPhiSmearingMatrixSB23Km.png");
delete hm23;

  auto hp34 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp34->Draw("col");
  histSmearFullKp34->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp34->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp34->GetXaxis()->SetNdivisions(505);
  histSmearFullKp34->GetYaxis()->SetNdivisions(505);
  histSmearFullKp34->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp34->Print("pPhiSmearingMatrixSB34Kp.pdf");
  hp34->Print("pPhiSmearingMatrixSB34Kp.png");
delete hp34;

  auto hm34 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm34->Draw("col");
  histSmearFullKm34->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm34->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm34->GetXaxis()->SetNdivisions(505);
  histSmearFullKm34->GetYaxis()->SetNdivisions(505);
  histSmearFullKm34->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm34->Print("pPhiSmearingMatrixSB34Km.pdf");
  hm34->Print("pPhiSmearingMatrixSB34Km.png");
delete hm34;

  auto hp45 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp45->Draw("col");
  histSmearFullKp45->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp45->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp45->GetXaxis()->SetNdivisions(505);
  histSmearFullKp45->GetYaxis()->SetNdivisions(505);
  histSmearFullKp45->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp45->Print("pPhiSmearingMatrixSB45Kp.pdf");
  hp45->Print("pPhiSmearingMatrixSB45Kp.png");
delete hp45;

  auto hm45 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm45->Draw("col");
  histSmearFullKm45->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm45->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm45->GetXaxis()->SetNdivisions(505);
  histSmearFullKm45->GetYaxis()->SetNdivisions(505);
  histSmearFullKm45->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm45->Print("pPhiSmearingMatrixSB45Km.pdf");
  hm45->Print("pPhiSmearingMatrixSB45Km.png");
delete hm45;

  auto hp56 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp56->Draw("col");
  histSmearFullKp56->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp56->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp56->GetXaxis()->SetNdivisions(505);
  histSmearFullKp56->GetYaxis()->SetNdivisions(505);
  histSmearFullKp56->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp56->Print("pPhiSmearingMatrixSB56Kp.pdf");
  hp56->Print("pPhiSmearingMatrixSB56Kp.png");
delete hp56;

  auto hm56 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm56->Draw("col");
  histSmearFullKm56->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm56->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm56->GetXaxis()->SetNdivisions(505);
  histSmearFullKm56->GetYaxis()->SetNdivisions(505);
  histSmearFullKm56->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm56->Print("pPhiSmearingMatrixSB56Km.pdf");
  hm56->Print("pPhiSmearingMatrixSB56Km.png");
delete hm56;

  auto hp67 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp67->Draw("col");
  histSmearFullKp67->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp67->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp67->GetXaxis()->SetNdivisions(505);
  histSmearFullKp67->GetYaxis()->SetNdivisions(505);
  histSmearFullKp67->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp67->Print("pPhiSmearingMatrixSB67Kp.pdf");
  hp67->Print("pPhiSmearingMatrixSB67Kp.png");
delete hp67;

  auto hm67 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm67->Draw("col");
  histSmearFullKm67->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm67->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm67->GetXaxis()->SetNdivisions(505);
  histSmearFullKm67->GetYaxis()->SetNdivisions(505);
  histSmearFullKm67->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm67->Print("pPhiSmearingMatrixSB67Km.pdf");
  hm67->Print("pPhiSmearingMatrixSB67Km.png");
delete hm67;

  auto hp78 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp78->Draw("col");
  histSmearFullKp78->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp78->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp78->GetXaxis()->SetNdivisions(505);
  histSmearFullKp78->GetYaxis()->SetNdivisions(505);
  histSmearFullKp78->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp78->Print("pPhiSmearingMatrixSB78Kp.pdf");
  hp78->Print("pPhiSmearingMatrixSB78Kp.png");
delete hp78;

  auto hm78 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm78->Draw("col");
  histSmearFullKm78->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm78->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm78->GetXaxis()->SetNdivisions(505);
  histSmearFullKm78->GetYaxis()->SetNdivisions(505);
  histSmearFullKm78->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm78->Print("pPhiSmearingMatrixSB78Km.pdf");
  hm78->Print("pPhiSmearingMatrixSB78Km.png");
delete hm78;

  auto hp89 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp89->Draw("col");
  histSmearFullKp89->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp89->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp89->GetXaxis()->SetNdivisions(505);
  histSmearFullKp89->GetYaxis()->SetNdivisions(505);
  histSmearFullKp89->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp89->Print("pPhiSmearingMatrixSB89Kp.pdf");
  hp89->Print("pPhiSmearingMatrixSB89Kp.png");
delete hp89;

  auto hm89 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm89->Draw("col");
  histSmearFullKm89->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm89->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm89->GetXaxis()->SetNdivisions(505);
  histSmearFullKm89->GetYaxis()->SetNdivisions(505);
  histSmearFullKm89->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm89->Print("pPhiSmearingMatrixSB89Km.pdf");
  hm89->Print("pPhiSmearingMatrixSB89Km.png");
delete hm89;

  auto hp90 = new TCanvas("h", "h", 650, 550);
  histSmearFullKp90->Draw("col");
  histSmearFullKp90->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKp90->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKp90->GetXaxis()->SetNdivisions(505);
  histSmearFullKp90->GetYaxis()->SetNdivisions(505);
  histSmearFullKp90->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hp90->Print("pPhiSmearingMatrixSB90Kp.pdf");
  hp90->Print("pPhiSmearingMatrixSB90Kp.png");
delete hp90;

  auto hm90 = new TCanvas("h", "h", 650, 550);
  histSmearFullKm90->Draw("col");
  histSmearFullKm90->GetXaxis()->SetRangeUser(0, 500);
  histSmearFullKm90->GetYaxis()->SetRangeUser(0, 500);
  histSmearFullKm90->GetXaxis()->SetNdivisions(505);
  histSmearFullKm90->GetYaxis()->SetNdivisions(505);
  histSmearFullKm90->SetTitle(
      "; #it{k}*_{p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  hm90->Print("pPhiSmearingMatrixSB90Km.pdf");
  hm90->Print("pPhiSmearingMatrixSB90Km.png");
delete hm90;



  auto hr1 = (TH2F*) infile->Get("histSmearFullSBrightKp");
  auto hr2 = (TH2F*) infile->Get("histSmearFullSBrightKm");
  auto hl1 = (TH2F*) infile->Get("histSmearFullSBleftKp");
  auto hl2 = (TH2F*) infile->Get("histSmearFullSBleftKm");
  auto hp1 = (TH2F*) infile->Get("histSmearFullPhiPeakKp");
  auto hp2 = (TH2F*) infile->Get("histSmearFullPhiPeakKm");

  auto h12p = (TH2F*) infile->Get("histSmearFullKp12");
  auto h12m = (TH2F*) infile->Get("histSmearFullKm12");
  auto h23p = (TH2F*) infile->Get("histSmearFullKp23");
  auto h23m = (TH2F*) infile->Get("histSmearFullKm23");
  auto h34p = (TH2F*) infile->Get("histSmearFullKp34");
  auto h34m = (TH2F*) infile->Get("histSmearFullKm34");
  auto h45p = (TH2F*) infile->Get("histSmearFullKp45");
  auto h45m = (TH2F*) infile->Get("histSmearFullKm45");
  auto h56p = (TH2F*) infile->Get("histSmearFullKp56");
  auto h56m = (TH2F*) infile->Get("histSmearFullKm56");
  auto h67p = (TH2F*) infile->Get("histSmearFullKp67");
  auto h67m = (TH2F*) infile->Get("histSmearFullKm67");
  auto h78p = (TH2F*) infile->Get("histSmearFullKp78");
  auto h78m = (TH2F*) infile->Get("histSmearFullKm78");
  auto h89p = (TH2F*) infile->Get("histSmearFullKp89");
  auto h89m = (TH2F*) infile->Get("histSmearFullKm89");
  auto h90p = (TH2F*) infile->Get("histSmearFullKp90");
  auto h90m = (TH2F*) infile->Get("histSmearFullKm90");



  DreamPlot::SetStyleHisto(hr1);
  DreamPlot::SetStyleHisto(hr2);
  DreamPlot::SetStyleHisto(hl1);
  DreamPlot::SetStyleHisto(hl2);
  DreamPlot::SetStyleHisto(hp1);
  DreamPlot::SetStyleHisto(hp2);
  DreamPlot::SetStyleHisto(h12p);
  DreamPlot::SetStyleHisto(h12m);
  DreamPlot::SetStyleHisto(h23p);
  DreamPlot::SetStyleHisto(h23m);
  DreamPlot::SetStyleHisto(h34p);
  DreamPlot::SetStyleHisto(h34m);
  DreamPlot::SetStyleHisto(h45p);
  DreamPlot::SetStyleHisto(h45m);
  DreamPlot::SetStyleHisto(h56p);
  DreamPlot::SetStyleHisto(h56m);
  DreamPlot::SetStyleHisto(h67p);
  DreamPlot::SetStyleHisto(h67m);
  DreamPlot::SetStyleHisto(h78p);
  DreamPlot::SetStyleHisto(h78m);
  DreamPlot::SetStyleHisto(h89p);
  DreamPlot::SetStyleHisto(h89m);
  DreamPlot::SetStyleHisto(h90p);
  DreamPlot::SetStyleHisto(h90m);


  hp1->Divide(hp2);
  hl1->Divide(hl2);
  hr1->Divide(hr2);
  h12p->Divide(h12m);
  h23p->Divide(h23m);
  h34p->Divide(h34m);
  h45p->Divide(h45m);
  h56p->Divide(h56m);
  h67p->Divide(h67m);
  h78p->Divide(h78m);
  h89p->Divide(h89m);
  h90p->Divide(h90m);






    auto i = new TCanvas("i", "i", 650, 550);
    hp1->Draw("col");
    hp1->GetXaxis()->SetRangeUser(0, 500);
    hp1->GetYaxis()->SetRangeUser(0, 500);
    hp1->GetXaxis()->SetNdivisions(505);
    hp1->GetYaxis()->SetNdivisions(505);
    hp1->SetTitle(
        ";  Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus #phi} (MeV/#it{c})");
    i->Print("ratioPeak.pdf");
    i->Print("ratioPeak.png");

delete i;

    auto j = new TCanvas("j", "j", 650, 550);
    hl1->Draw("col");
    hl1->GetXaxis()->SetRangeUser(0, 500);
    hl1->GetYaxis()->SetRangeUser(0, 500);
    hl1->GetXaxis()->SetNdivisions(505);
    hl1->GetYaxis()->SetNdivisions(505);
    hl1->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    j->Print("ratioSBleft.pdf");
    j->Print("ratioSBleft.png");

delete j;
    auto k = new TCanvas("k", "k", 650, 550);
    hr1->Draw("col");
    hr1->GetXaxis()->SetRangeUser(0, 500);
    hr1->GetYaxis()->SetRangeUser(0, 500);
    hr1->GetXaxis()->SetNdivisions(505);
    hr1->GetYaxis()->SetNdivisions(505);
    hr1->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k->Print("ratioSBright.pdf");
    k->Print("ratioSBright.png");
delete k;


    auto k12 = new TCanvas("k", "k", 650, 550);
    h12p->Draw("col");
    h12p->GetXaxis()->SetRangeUser(0, 500);
    h12p->GetYaxis()->SetRangeUser(0, 500);
    h12p->GetXaxis()->SetNdivisions(505);
    h12p->GetYaxis()->SetNdivisions(505);
    h12p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k12->Print("ratio12.pdf");
    k12->Print("ratio12.png");
delete k12;

    auto k23 = new TCanvas("k", "k", 650, 550);
    h23p->Draw("col");
    h23p->GetXaxis()->SetRangeUser(0, 500);
    h23p->GetYaxis()->SetRangeUser(0, 500);
    h23p->GetXaxis()->SetNdivisions(505);
    h23p->GetYaxis()->SetNdivisions(505);
    h23p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k23->Print("ratio23.pdf");
    k23->Print("ratio23.png");
delete k23;
    auto k34 = new TCanvas("k", "k", 650, 550);
    h34p->Draw("col");
    h34p->GetXaxis()->SetRangeUser(0, 500);
    h34p->GetYaxis()->SetRangeUser(0, 500);
    h34p->GetXaxis()->SetNdivisions(505);
    h34p->GetYaxis()->SetNdivisions(505);
    h34p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k34->Print("ratio34.pdf");
    k34->Print("ratio34.png");
delete k34;
    auto k45 = new TCanvas("k", "k", 650, 550);
    h45p->Draw("col");
    h45p->GetXaxis()->SetRangeUser(0, 500);
    h45p->GetYaxis()->SetRangeUser(0, 500);
    h45p->GetXaxis()->SetNdivisions(505);
    h45p->GetYaxis()->SetNdivisions(505);
    h45p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k45->Print("ratio45.pdf");
    k45->Print("ratio45.png");
delete k45;
    auto k56 = new TCanvas("k", "k", 650, 550);
    h56p->Draw("col");
    h56p->GetXaxis()->SetRangeUser(0, 500);
    h56p->GetYaxis()->SetRangeUser(0, 500);
    h56p->GetXaxis()->SetNdivisions(505);
    h56p->GetYaxis()->SetNdivisions(505);
    h56p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k56->Print("ratio56.pdf");
    k56->Print("ratio56.png");
delete k56;
    auto k67 = new TCanvas("k", "k", 650, 550);
    h67p->Draw("col");
    h67p->GetXaxis()->SetRangeUser(0, 500);
    h67p->GetYaxis()->SetRangeUser(0, 500);
    h67p->GetXaxis()->SetNdivisions(505);
    h67p->GetYaxis()->SetNdivisions(505);
    h67p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k67->Print("ratio67.pdf");
    k67->Print("ratio67.png");
delete k67;
    auto k78 = new TCanvas("k", "k", 650, 550);
    h78p->Draw("col");
    h78p->GetXaxis()->SetRangeUser(0, 500);
    h78p->GetYaxis()->SetRangeUser(0, 500);
    h78p->GetXaxis()->SetNdivisions(505);
    h78p->GetYaxis()->SetNdivisions(505);
    h78p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k78->Print("ratio78.pdf");
    k78->Print("ratio78.png");
delete k78;
    auto k89 = new TCanvas("k", "k", 650, 550);
    h89p->Draw("col");
    h89p->GetXaxis()->SetRangeUser(0, 500);
    h89p->GetYaxis()->SetRangeUser(0, 500);
    h89p->GetXaxis()->SetNdivisions(505);
    h89p->GetYaxis()->SetNdivisions(505);
    h89p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k89->Print("ratio89.pdf");
    k89->Print("ratio89.png");
delete k89;
    auto k90 = new TCanvas("k", "k", 650, 550);
    h90p->Draw("col");
    h90p->GetXaxis()->SetRangeUser(0, 500);
    h90p->GetYaxis()->SetRangeUser(0, 500);
    h90p->GetXaxis()->SetNdivisions(505);
    h90p->GetYaxis()->SetNdivisions(505);
    h90p->SetTitle(
        "; Ratio_{p #minus K^{+}/p #minus K^{-}}; Ratio_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    k90->Print("ratio90.pdf");
    k90->Print("ratio90.png");
delete k90;



    outfile->cd();

    hp1->Write("ratioPeak");
    hl1->Write("ratioSBleft");
    hr1->Write("ratioSBright");
    h12p->Write("ratio12");
    h23p->Write("ratio23");
    h34p->Write("ratio34");
    h45p->Write("ratio45");
    h56p->Write("ratio56");
    h67p->Write("ratio67");
    h78p->Write("ratio78");
    h89p->Write("ratio89");
    h90p->Write("ratio90");




     // outfile->Close();



  auto hhr1 = (TH2F*) infile->Get("histSmearFullSBrightKp");
  auto hhr2 = (TH2F*) infile->Get("histSmearFullSBrightKm");
  auto hhl1 = (TH2F*) infile->Get("histSmearFullSBleftKp");
  auto hhl2 = (TH2F*) infile->Get("histSmearFullSBleftKm");
  auto hhp1 = (TH2F*) infile->Get("histSmearFullPhiPeakKp");
  auto hhp2 = (TH2F*) infile->Get("histSmearFullPhiPeakKm");

  auto hh12p = (TH2F*) infile->Get("histSmearFullKp12");
  auto hh12m = (TH2F*) infile->Get("histSmearFullKm12");
  auto hh23p = (TH2F*) infile->Get("histSmearFullKp23");
  auto hh23m = (TH2F*) infile->Get("histSmearFullKm23");
  auto hh34p = (TH2F*) infile->Get("histSmearFullKp34");
  auto hh34m = (TH2F*) infile->Get("histSmearFullKm34");
  auto hh45p = (TH2F*) infile->Get("histSmearFullKp45");
  auto hh45m = (TH2F*) infile->Get("histSmearFullKm45");
  auto hh56p = (TH2F*) infile->Get("histSmearFullKp56");
  auto hh56m = (TH2F*) infile->Get("histSmearFullKm56");
  auto hh67p = (TH2F*) infile->Get("histSmearFullKp67");
  auto hh67m = (TH2F*) infile->Get("histSmearFullKm67");
  auto hh78p = (TH2F*) infile->Get("histSmearFullKp78");
  auto hh78m = (TH2F*) infile->Get("histSmearFullKm78");
  auto hh89p = (TH2F*) infile->Get("histSmearFullKp89");
  auto hh89m = (TH2F*) infile->Get("histSmearFullKm89");
  auto hh90p = (TH2F*) infile->Get("histSmearFullKp90");
  auto hh90m = (TH2F*) infile->Get("histSmearFullKm90");

//  auto ratPeak = histSmearFullPhiPeakKp;
//  ratPeak->Divide(histSmearFullPhiPeakKm);
//  auto ratlSB=histSmearFullSBleftKp;
//  ratlSB->Divide(histSmearFullSBleftKm);
//  auto ratrSB=histSmearFullSBrightKp;
//  ratrSB->Divide(histSmearFullSBrightKm);

  DreamPlot::SetStyleHisto(hhr1);
  DreamPlot::SetStyleHisto(hhr2);
  DreamPlot::SetStyleHisto(hhl1);
  DreamPlot::SetStyleHisto(hhl2);
  DreamPlot::SetStyleHisto(hhp1);
  DreamPlot::SetStyleHisto(hhp2);

  DreamPlot::SetStyleHisto(hh12p);
  DreamPlot::SetStyleHisto(hh12m);
  DreamPlot::SetStyleHisto(hh23p);
  DreamPlot::SetStyleHisto(hh23m);
  DreamPlot::SetStyleHisto(hh34p);
  DreamPlot::SetStyleHisto(hh34m);
  DreamPlot::SetStyleHisto(hh45p);
  DreamPlot::SetStyleHisto(hh45m);
  DreamPlot::SetStyleHisto(hh56p);
  DreamPlot::SetStyleHisto(hh56m);
  DreamPlot::SetStyleHisto(hh67p);
  DreamPlot::SetStyleHisto(hh67m);
  DreamPlot::SetStyleHisto(hh78p);
  DreamPlot::SetStyleHisto(hh78m);
  DreamPlot::SetStyleHisto(hh89p);
  DreamPlot::SetStyleHisto(hh89m);
  DreamPlot::SetStyleHisto(hh90p);
  DreamPlot::SetStyleHisto(hh90m);

  hhp1->Add(hhp2);
  hhl1->Add(hhl2);
  hhr1->Add(hhr2);


  hh12p->Add(hh12m);
  hh23p->Add(hh23m);
  hh34p->Add(hh34m);
  hh45p->Add(hh45m);
  hh56p->Add(hh56m);
  hh67p->Add(hh67m);
  hh78p->Add(hh78m);
  hh89p->Add(hh89m);
  hh90p->Add(hh90m);




  outfile->cd();




    auto l = new TCanvas("l", "l", 650, 550);
    hhp1->Draw("col");
    hhp1->GetXaxis()->SetRangeUser(0, 500);
    hhp1->GetYaxis()->SetRangeUser(0, 500);
    hhp1->GetXaxis()->SetNdivisions(505);
    hhp1->GetYaxis()->SetNdivisions(505);
    hhp1->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus #phi} (MeV/#it{c})");
    l->Print("totPeak.pdf");
    l->Print("totPeak.png");

    delete l;

    auto m = new TCanvas("m", "m", 650, 550);
    hhl1->Draw("col");
    hhl1->GetXaxis()->SetRangeUser(0, 500);
    hhl1->GetYaxis()->SetRangeUser(0, 500);
    hhl1->GetXaxis()->SetNdivisions(505);
    hhl1->GetYaxis()->SetNdivisions(505);
    hhl1->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m->Print("totSBleft.pdf");
    m->Print("totSBleft.png");

    delete m;
    auto nn = new TCanvas("n", "n", 650, 550);
    hhr1->Draw("col");
    hhr1->GetXaxis()->SetRangeUser(0, 500);
    hhr1->GetYaxis()->SetRangeUser(0, 500);
    hhr1->GetXaxis()->SetNdivisions(505);
    hhr1->GetYaxis()->SetNdivisions(505);
    hhr1->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    nn->Print("totSBright.pdf");
    nn->Print("totSBright.png");

delete nn;

    auto m12 = new TCanvas("m", "m", 650, 550);
    hh12p->Draw("col");
    hh12p->GetXaxis()->SetRangeUser(0, 500);
    hh12p->GetYaxis()->SetRangeUser(0, 500);
    hh12p->GetXaxis()->SetNdivisions(505);
    hh12p->GetYaxis()->SetNdivisions(505);
    hh12p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m12->Print("tot12.pdf");
    m12->Print("tot12left.png");

    delete m12;
    auto m23 = new TCanvas("m", "m", 650, 550);
    hh23p->Draw("col");
    hh23p->GetXaxis()->SetRangeUser(0, 500);
    hh23p->GetYaxis()->SetRangeUser(0, 500);
    hh23p->GetXaxis()->SetNdivisions(505);
    hh23p->GetYaxis()->SetNdivisions(505);
    hh23p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m23->Print("tot23.pdf");
    m23->Print("tot23left.png");

    delete m23;
    auto m34 = new TCanvas("m", "m", 650, 550);
    hh34p->Draw("col");
    hh34p->GetXaxis()->SetRangeUser(0, 500);
    hh34p->GetYaxis()->SetRangeUser(0, 500);
    hh34p->GetXaxis()->SetNdivisions(505);
    hh34p->GetYaxis()->SetNdivisions(505);
    hh34p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m34->Print("tot34.pdf");
    m34->Print("tot34left.png");

    delete m34;
    auto m45 = new TCanvas("m", "m", 650, 550);
    hh45p->Draw("col");
    hh45p->GetXaxis()->SetRangeUser(0, 500);
    hh45p->GetYaxis()->SetRangeUser(0, 500);
    hh45p->GetXaxis()->SetNdivisions(505);
    hh45p->GetYaxis()->SetNdivisions(505);
    hh45p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m45->Print("tot45.pdf");
    m45->Print("tot45left.png");

    delete m45;
    auto m56 = new TCanvas("m", "m", 650, 550);
    hh56p->Draw("col");
    hh56p->GetXaxis()->SetRangeUser(0, 500);
    hh56p->GetYaxis()->SetRangeUser(0, 500);
    hh56p->GetXaxis()->SetNdivisions(505);
    hh56p->GetYaxis()->SetNdivisions(505);
    hh56p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m56->Print("tot56.pdf");
    m56->Print("tot56left.png");

    delete m56;
    auto m67 = new TCanvas("m", "m", 650, 550);
    hh67p->Draw("col");
    hh67p->GetXaxis()->SetRangeUser(0, 500);
    hh67p->GetYaxis()->SetRangeUser(0, 500);
    hh67p->GetXaxis()->SetNdivisions(505);
    hh67p->GetYaxis()->SetNdivisions(505);
    hh67p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m67->Print("tot67.pdf");
    m67->Print("tot67left.png");

    delete m67;
    auto m78 = new TCanvas("m", "m", 650, 550);
    hh78p->Draw("col");
    hh78p->GetXaxis()->SetRangeUser(0, 500);
    hh78p->GetYaxis()->SetRangeUser(0, 500);
    hh78p->GetXaxis()->SetNdivisions(505);
    hh78p->GetYaxis()->SetNdivisions(505);
    hh78p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m78->Print("tot78.pdf");
    m78->Print("tot78left.png");

    delete m78;
    auto m89 = new TCanvas("m", "m", 650, 550);
    hh89p->Draw("col");
    hh89p->GetXaxis()->SetRangeUser(0, 500);
    hh89p->GetYaxis()->SetRangeUser(0, 500);
    hh89p->GetXaxis()->SetNdivisions(505);
    hh89p->GetYaxis()->SetNdivisions(505);
    hh89p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m89->Print("tot89.pdf");
    m89->Print("tot89left.png");

    delete m89;
    auto m90 = new TCanvas("m", "m", 650, 550);
    hh90p->Draw("col");
    hh90p->GetXaxis()->SetRangeUser(0, 500);
    hh90p->GetYaxis()->SetRangeUser(0, 500);
    hh90p->GetXaxis()->SetNdivisions(505);
    hh90p->GetYaxis()->SetNdivisions(505);
    hh90p->SetTitle(
        "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
    m90->Print("tot90.pdf");
    m90->Print("tot90left.png");

    delete m90;




    hhp1->Write("totalPeak");
    hhl1->Write("totalSBleft");
    hhr1->Write("totalSBright");
    hh12p->Write("total12");
    hh23p->Write("total23");
    hh34p->Write("total34");
    hh45p->Write("total45");
    hh56p->Write("total56");
    hh67p->Write("total67");
    hh78p->Write("total78");
    hh89p->Write("total89");
    hh90p->Write("total90");




//    auto infile2 = TFile::Open("/home/emma/FemtoPhiHM_ROTMC/SmearingMCALL/SmearSideband.root");

//    auto hhr1b = (TH2F*) infile2->Get("histSmearFullSBrightKp");
//    auto hhr2b = (TH2F*) infile2->Get("histSmearFullSBrightKm");
//    auto hhl1b = (TH2F*) infile2->Get("histSmearFullSBleftKp");
//    auto hhl2b = (TH2F*) infile2->Get("histSmearFullSBleftKm");
//    auto hhp1b = (TH2F*) infile2->Get("histSmearFullPhiPeakKp");
//    auto hhp2b = (TH2F*) infile2->Get("histSmearFullPhiPeakKm");

//    auto hh12pb = (TH2F*) infile2->Get("histSmearFullKp12");
//    auto hh12mb = (TH2F*) infile2->Get("histSmearFullKm12");
//    auto hh23pb = (TH2F*) infile2->Get("histSmearFullKp23");
//    auto hh23mb = (TH2F*) infile2->Get("histSmearFullKm23");
//    auto hh34pb = (TH2F*) infile2->Get("histSmearFullKp34");
//    auto hh34mb = (TH2F*) infile2->Get("histSmearFullKm34");
//    auto hh45pb = (TH2F*) infile2->Get("histSmearFullKp45");
//    auto hh45mb = (TH2F*) infile2->Get("histSmearFullKm45");
//    auto hh56pb = (TH2F*) infile2->Get("histSmearFullKp56");
//    auto hh56mb = (TH2F*) infile2->Get("histSmearFullKm56");
//    auto hh67pb = (TH2F*) infile2->Get("histSmearFullKp67");
//    auto hh67mb = (TH2F*) infile2->Get("histSmearFullKm67");
//    auto hh78pb = (TH2F*) infile2->Get("histSmearFullKp78");
//    auto hh78mb = (TH2F*) infile2->Get("histSmearFullKm78");
//    auto hh89pb = (TH2F*) infile2->Get("histSmearFullKp89");
//    auto hh89mb = (TH2F*) infile2->Get("histSmearFullKm89");
//    auto hh90pb = (TH2F*) infile2->Get("histSmearFullKp90");
//    auto hh90mb = (TH2F*) infile2->Get("histSmearFullKm90");

//  //  auto ratPeak = histSmearFullPhiPeakKp;
//  //  ratPeak->Divide(histSmearFullPhiPeakKm);
//  //  auto ratlSB=histSmearFullSBleftKp;
//  //  ratlSB->Divide(histSmearFullSBleftKm);
//  //  auto ratrSB=histSmearFullSBrightKp;
//  //  ratrSB->Divide(histSmearFullSBrightKm);

//    DreamPlot::SetStyleHisto(hhr1b);
//    DreamPlot::SetStyleHisto(hhr2b);
//    DreamPlot::SetStyleHisto(hhl1b);
//    DreamPlot::SetStyleHisto(hhl2b);
//    DreamPlot::SetStyleHisto(hhp1b);
//    DreamPlot::SetStyleHisto(hhp2b);

//    DreamPlot::SetStyleHisto(hh12pb);
//    DreamPlot::SetStyleHisto(hh12mb);
//    DreamPlot::SetStyleHisto(hh23pb);
//    DreamPlot::SetStyleHisto(hh23mb);
//    DreamPlot::SetStyleHisto(hh34pb);
//    DreamPlot::SetStyleHisto(hh34mb);
//    DreamPlot::SetStyleHisto(hh45pb);
//    DreamPlot::SetStyleHisto(hh45mb);
//    DreamPlot::SetStyleHisto(hh56pb);
//    DreamPlot::SetStyleHisto(hh56mb);
//    DreamPlot::SetStyleHisto(hh67pb);
//    DreamPlot::SetStyleHisto(hh67mb);
//    DreamPlot::SetStyleHisto(hh78pb);
//    DreamPlot::SetStyleHisto(hh78mb);
//    DreamPlot::SetStyleHisto(hh89pb);
//    DreamPlot::SetStyleHisto(hh89mb);
//    DreamPlot::SetStyleHisto(hh90pb);
//    DreamPlot::SetStyleHisto(hh90mb);

//    hhp1b->Add(hhp2b);
//    hhl1b->Add(hhl2b);
//    hhr1b->Add(hhr2b);


//    hh12pb->Add(hh12mb);
//    hh23pb->Add(hh23mb);
//    hh34pb->Add(hh34mb);
//    hh45pb->Add(hh45mb);
//    hh56pb->Add(hh56mb);
//    hh67pb->Add(hh67mb);
//    hh78pb->Add(hh78mb);
//    hh89pb->Add(hh89mb);
//    hh90pb->Add(hh90mb);




//    outfile->cd();




//      auto lb = new TCanvas("l", "l", 650, 550);
//      hhp1b->Draw("col");
//      hhp1b->GetXaxis()->SetRangeUser(0, 500);
//      hhp1b->GetYaxis()->SetRangeUser(0, 500);
//      hhp1b->GetXaxis()->SetNdivisions(505);
//      hhp1b->GetYaxis()->SetNdivisions(505);
//      hhp1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus #phi} (MeV/#it{c})");
//      lb->Print("totPeakmc.pdf");
//      lb->Print("totPeakmc.png");

//      delete lb;

//      auto mb = new TCanvas("m", "m", 650, 550);
//      hhl1b->Draw("col");
//      hhl1b->GetXaxis()->SetRangeUser(0, 500);
//      hhl1b->GetYaxis()->SetRangeUser(0, 500);
//      hhl1b->GetXaxis()->SetNdivisions(505);
//      hhl1b->GetYaxis()->SetNdivisions(505);
//      hhl1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      mb->Print("totSBleftmc.pdf");
//      mb->Print("totSBleftmc.png");

//      delete mb;
//      auto nnb = new TCanvas("n", "n", 650, 550);
//      hhr1b->Draw("col");
//      hhr1b->GetXaxis()->SetRangeUser(0, 500);
//      hhr1b->GetYaxis()->SetRangeUser(0, 500);
//      hhr1b->GetXaxis()->SetNdivisions(505);
//      hhr1b->GetYaxis()->SetNdivisions(505);
//      hhr1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      nnb->Print("totSBrightmc.pdf");
//      nnb->Print("totSBrightmc.png");

//  delete nnb;

//      auto m12b = new TCanvas("m", "m", 650, 550);
//      hh12pb->Draw("col");
//      hh12pb->GetXaxis()->SetRangeUser(0, 500);
//      hh12pb->GetYaxis()->SetRangeUser(0, 500);
//      hh12pb->GetXaxis()->SetNdivisions(505);
//      hh12pb->GetYaxis()->SetNdivisions(505);
//      hh12pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m12b->Print("tot12mc.pdf");
//      m12b->Print("tot12mc.png");

//      delete m12b;
//      auto m23b = new TCanvas("m", "m", 650, 550);
//      hh23pb->Draw("col");
//      hh23pb->GetXaxis()->SetRangeUser(0, 500);
//      hh23pb->GetYaxis()->SetRangeUser(0, 500);
//      hh23pb->GetXaxis()->SetNdivisions(505);
//      hh23pb->GetYaxis()->SetNdivisions(505);
//      hh23pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m23b->Print("tot23mc.pdf");
//      m23b->Print("tot23mc.png");

//      delete m23b;
//      auto m34b = new TCanvas("m", "m", 650, 550);
//      hh34pb->Draw("col");
//      hh34pb->GetXaxis()->SetRangeUser(0, 500);
//      hh34pb->GetYaxis()->SetRangeUser(0, 500);
//      hh34pb->GetXaxis()->SetNdivisions(505);
//      hh34pb->GetYaxis()->SetNdivisions(505);
//      hh34pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m34b->Print("tot34mc.pdf");
//      m34b->Print("tot34mc.png");

//      delete m34b;
//      auto m45b = new TCanvas("m", "m", 650, 550);
//      hh45pb->Draw("col");
//      hh45pb->GetXaxis()->SetRangeUser(0, 500);
//      hh45pb->GetYaxis()->SetRangeUser(0, 500);
//      hh45pb->GetXaxis()->SetNdivisions(505);
//      hh45pb->GetYaxis()->SetNdivisions(505);
//      hh45pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m45b->Print("tot45mc.pdf");
//      m45b->Print("tot45mc.png");

//      delete m45b;
//      auto m56b = new TCanvas("m", "m", 650, 550);
//      hh56pb->Draw("col");
//      hh56pb->GetXaxis()->SetRangeUser(0, 500);
//      hh56pb->GetYaxis()->SetRangeUser(0, 500);
//      hh56pb->GetXaxis()->SetNdivisions(505);
//      hh56pb->GetYaxis()->SetNdivisions(505);
//      hh56pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m56b->Print("tot56mc.pdf");
//      m56b->Print("tot56mc.png");

//      delete m56b;
//      auto m67b = new TCanvas("m", "m", 650, 550);
//      hh67pb->Draw("col");
//      hh67pb->GetXaxis()->SetRangeUser(0, 500);
//      hh67pb->GetYaxis()->SetRangeUser(0, 500);
//      hh67pb->GetXaxis()->SetNdivisions(505);
//      hh67pb->GetYaxis()->SetNdivisions(505);
//      hh67pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m67b->Print("tot67mc.pdf");
//      m67b->Print("tot67mc.png");

//      delete m67b;
//      auto m78b = new TCanvas("m", "m", 650, 550);
//      hh78pb->Draw("col");
//      hh78pb->GetXaxis()->SetRangeUser(0, 500);
//      hh78pb->GetYaxis()->SetRangeUser(0, 500);
//      hh78pb->GetXaxis()->SetNdivisions(505);
//      hh78pb->GetYaxis()->SetNdivisions(505);
//      hh78pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m78b->Print("tot78mc.pdf");
//      m78b->Print("tot78mc.png");

//      delete m78b;
//      auto m89b = new TCanvas("m", "m", 650, 550);
//      hh89pb->Draw("col");
//      hh89pb->GetXaxis()->SetRangeUser(0, 500);
//      hh89pb->GetYaxis()->SetRangeUser(0, 500);
//      hh89pb->GetXaxis()->SetNdivisions(505);
//      hh89pb->GetYaxis()->SetNdivisions(505);
//      hh89pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m89b->Print("tot89mc.pdf");
//      m89b->Print("tot89mc.png");

//      delete m89b;
//      auto m90b = new TCanvas("m", "m", 650, 550);
//      hh90pb->Draw("col");
//      hh90pb->GetXaxis()->SetRangeUser(0, 500);
//      hh90pb->GetYaxis()->SetRangeUser(0, 500);
//      hh90pb->GetXaxis()->SetNdivisions(505);
//      hh90pb->GetYaxis()->SetNdivisions(505);
//      hh90pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m90b->Print("tot90mc.pdf");
//      m90b->Print("tot90mc.png");

//      delete m90b;




//      hhp1b->Write("totalPeakmc");
//      hhl1b->Write("totalSBleftmc");
//      hhr1b->Write("totalSBrightmc");
//      hh12pb->Write("total12mc");
//      hh23pb->Write("total23mc");
//      hh34pb->Write("total34mc");
//      hh45pb->Write("total45mc");
//      hh56pb->Write("total56mc");
//      hh67pb->Write("total67mc");
//      hh78pb->Write("total78mc");
//      hh89pb->Write("total89mc");
//      hh90pb->Write("total90mc");


//      hhp1b->Divide(hhp1);
//      hhl1b->Divide(hhl1);
//      hhr1b->Divide(hhr1);
//      hh12pb->Divide(hh12p);
//      hh23pb->Divide(hh23p);
//      hh34pb->Divide(hh34p);
//      hh45pb->Divide(hh45p);
//      hh56pb->Divide(hh56p);
//      hh67pb->Divide(hh67p);
//      hh78pb->Divide(hh78p);
//      hh89pb->Divide(hh89p);
//      hh90pb->Divide(hh90p);


//      auto lb2 = new TCanvas("l", "l", 650, 550);
//      hhp1b->Draw("col");
//      hhp1b->GetXaxis()->SetRangeUser(0, 500);
//      hhp1b->GetYaxis()->SetRangeUser(0, 500);
//      hhp1b->GetXaxis()->SetNdivisions(505);
//      hhp1b->GetYaxis()->SetNdivisions(505);
//      hhp1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus #phi} (MeV/#it{c})");
//      lb2->Print("RATpeak.pdf");
//      lb2->Print("RATpeak.png");

//      delete lb2;

//      auto mb2 = new TCanvas("m", "m", 650, 550);
//      hhl1b->Draw("col");
//      hhl1b->GetXaxis()->SetRangeUser(0, 500);
//      hhl1b->GetYaxis()->SetRangeUser(0, 500);
//      hhl1b->GetXaxis()->SetNdivisions(505);
//      hhl1b->GetYaxis()->SetNdivisions(505);
//      hhl1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      mb2->Print("RATSBleftmc.pdf");
//      mb2->Print("RATSBleftmc.png");

//      delete mb2;
//      auto nnb2 = new TCanvas("n", "n", 650, 550);
//      hhr1b->Draw("col");
//      hhr1b->GetXaxis()->SetRangeUser(0, 500);
//      hhr1b->GetYaxis()->SetRangeUser(0, 500);
//      hhr1b->GetXaxis()->SetNdivisions(505);
//      hhr1b->GetYaxis()->SetNdivisions(505);
//      hhr1b->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      nnb2->Print("RATSBrightmc.pdf");
//      nnb2->Print("RATSBrightmc.png");

//  delete nnb2;

//      auto m12b2 = new TCanvas("m", "m", 650, 550);
//      hh12pb->Draw("col");
//      hh12pb->GetXaxis()->SetRangeUser(0, 500);
//      hh12pb->GetYaxis()->SetRangeUser(0, 500);
//      hh12pb->GetXaxis()->SetNdivisions(505);
//      hh12pb->GetYaxis()->SetNdivisions(505);
//      hh12pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m12b2->Print("RAT12mc.pdf");
//      m12b2->Print("RAT12mc.png");

//      delete m12b2;
//      auto m23b2 = new TCanvas("m", "m", 650, 550);
//      hh23pb->Draw("col");
//      hh23pb->GetXaxis()->SetRangeUser(0, 500);
//      hh23pb->GetYaxis()->SetRangeUser(0, 500);
//      hh23pb->GetXaxis()->SetNdivisions(505);
//      hh23pb->GetYaxis()->SetNdivisions(505);
//      hh23pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m23b2->Print("RAT23mc.pdf");
//      m23b2->Print("RAT23mc.png");

//      delete m23b2;
//      auto m34b2 = new TCanvas("m", "m", 650, 550);
//      hh34pb->Draw("col");
//      hh34pb->GetXaxis()->SetRangeUser(0, 500);
//      hh34pb->GetYaxis()->SetRangeUser(0, 500);
//      hh34pb->GetXaxis()->SetNdivisions(505);
//      hh34pb->GetYaxis()->SetNdivisions(505);
//      hh34pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m34b2->Print("RAT34mc.pdf");
//      m34b2->Print("RAT34mc.png");

//      delete m34b2;
//      auto m45b2 = new TCanvas("m", "m", 650, 550);
//      hh45pb->Draw("col");
//      hh45pb->GetXaxis()->SetRangeUser(0, 500);
//      hh45pb->GetYaxis()->SetRangeUser(0, 500);
//      hh45pb->GetXaxis()->SetNdivisions(505);
//      hh45pb->GetYaxis()->SetNdivisions(505);
//      hh45pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m45b2->Print("RAT45mc.pdf");
//      m45b2->Print("RAT45mc.png");

//      delete m45b2;
//      auto m56b2 = new TCanvas("m", "m", 650, 550);
//      hh56pb->Draw("col");
//      hh56pb->GetXaxis()->SetRangeUser(0, 500);
//      hh56pb->GetYaxis()->SetRangeUser(0, 500);
//      hh56pb->GetXaxis()->SetNdivisions(505);
//      hh56pb->GetYaxis()->SetNdivisions(505);
//      hh56pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m56b2->Print("RAT56mc.pdf");
//      m56b2->Print("RAT56mc.png");

//      delete m56b2;
//      auto m67b2 = new TCanvas("m", "m", 650, 550);
//      hh67pb->Draw("col");
//      hh67pb->GetXaxis()->SetRangeUser(0, 500);
//      hh67pb->GetYaxis()->SetRangeUser(0, 500);
//      hh67pb->GetXaxis()->SetNdivisions(505);
//      hh67pb->GetYaxis()->SetNdivisions(505);
//      hh67pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m67b2->Print("RAT67mc.pdf");
//      m67b2->Print("RAT67mc.png");

//      delete m67b2;
//      auto m78b2 = new TCanvas("m", "m", 650, 550);
//      hh78pb->Draw("col");
//      hh78pb->GetXaxis()->SetRangeUser(0, 500);
//      hh78pb->GetYaxis()->SetRangeUser(0, 500);
//      hh78pb->GetXaxis()->SetNdivisions(505);
//      hh78pb->GetYaxis()->SetNdivisions(505);
//      hh78pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m78b2->Print("RAT78mc.pdf");
//      m78b2->Print("RAT78mc.png");

//      delete m78b2;
//      auto m89b2 = new TCanvas("m", "m", 650, 550);
//      hh89pb->Draw("col");
//      hh89pb->GetXaxis()->SetRangeUser(0, 500);
//      hh89pb->GetYaxis()->SetRangeUser(0, 500);
//      hh89pb->GetXaxis()->SetNdivisions(505);
//      hh89pb->GetYaxis()->SetNdivisions(505);
//      hh89pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m89b2->Print("RAT89mc.pdf");
//      m89b2->Print("RAT89mc.png");

//      delete m89b2;
//      auto m90b2 = new TCanvas("m", "m", 650, 550);
//      hh90pb->Draw("col");
//      hh90pb->GetXaxis()->SetRangeUser(0, 500);
//      hh90pb->GetYaxis()->SetRangeUser(0, 500);
//      hh90pb->GetXaxis()->SetNdivisions(505);
//      hh90pb->GetYaxis()->SetNdivisions(505);
//      hh90pb->SetTitle(
//          "; #it{k}*_{p #minus K^{+}+p #minus K^{-}} (MeV/#it{c}); #it{k}*_{p #minus (K^{+}K^{-})} (MeV/#it{c})");
//      m90b2->Print("RAT90mc.pdf");
//      m90b2->Print("RAT90mc.png");

//      delete m90b2;















 TH1F* phitruth=CalcCFMC("HMPhi",0);

 TH1F* SBleft=CalcCF("HMPhi",1);
 TH1F* SBright=CalcCF("HMPhi",2);
 TH1F* Peak=CalcCF("HMPhi",0);


  phitruth->SetMarkerStyle(20);
  phitruth->SetMarkerSize(1);
  phitruth->SetMarkerColor(kBlack);
  phitruth->SetLineColor(kBlack); 
  phitruth->SetLineWidth(2);

  SBleft->SetMarkerStyle(20);
  SBleft->SetMarkerSize(1);
  SBleft->SetMarkerColor(kBlack);
  SBleft->SetLineColor(kBlack);
  SBleft->SetLineWidth(2);

  SBright->SetMarkerStyle(20);
  SBright->SetMarkerSize(1);
  SBright->SetMarkerColor(kBlack);
  SBright->SetLineColor(kBlack);
  SBright->SetLineWidth(2);

  Peak->SetMarkerStyle(20);
  Peak->SetMarkerSize(1);
  Peak->SetMarkerColor(kBlack);
  Peak->SetLineColor(kBlack);
  Peak->SetLineWidth(2);


 int mPairs = 4;
 auto fit1 = new TF1("Fit", [&](double *x, double *p) {

          return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

         }, 0, 500, mPairs);
 auto fit2 = new TF1("Fit", [&](double *x, double *p) {

          return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

         }, 0, 500, mPairs);
 auto fit3 = new TF1("Fit", [&](double *x, double *p) {

          return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

         }, 0, 500, mPairs);


// phitruth->Fit(fit1,"F", "", 5, 400);

// SBleft->Fit(fit2,"F", "", 5, 360);
// SBright->Fit(fit3,"F", "", 5, 360);



 phitruth->Fit(fit1,"F", "", 5, 400);

 SBleft->Fit(fit2,"F", "", 5, 500);
 SBright->Fit(fit3,"F", "", 5, 500);

 const int binwidth = phitruth->GetBinWidth(1);
 const int NumMomBins = int(1000 / binwidth);
 float a=fit1->GetParameter(0);
 float b=fit1->GetParameter(1);
 float cc=fit1->GetParameter(2);
 float dd=fit1->GetParameter(3);


// TGraph fitfun;

// std::cout<<"lol222"<<std::endl;

// fitfun.SetName("GenuineCF");
// std::cout<<"lol3333"<<std::endl;

// for (int i = 0; i < NumMomBins; ++i) {

//   float f1=(a+b*i+cc*i*i+dd*i*i*i);

//   std::cout<<i<<" f1: "<<f1<<std::endl;
//   std::cout<<"a: "<<a<<"b: "<<b<<"c: "<<cc<<std::endl;

//   fitfun.SetPoint(i, i*binwidth, f1);  // genuine p-Phi CF with the parameters obtained in the fit

//}

 Int_t n =NumMomBins;
 Double_t x[n], y[n];
 for (Int_t i=0; i<n; i++) {
    x[i] = i*binwidth;
    y[i] = a+b*x[i]+cc*x[i]*x[i]+dd*x[i]*x[i]*x[i];
 }
 TGraph *fitfun = new TGraph (n, x, y);

 const int binwidth2 = SBleft->GetBinWidth(1);
 const int NumMomBins2 = int(1000 / binwidth2);
 float a2=fit2->GetParameter(0);
 float b2=fit2->GetParameter(1);
 float cc2=fit2->GetParameter(2);
 float dd2=fit2->GetParameter(3);

 Int_t n2 =NumMomBins2;
 Double_t x2[n2], y2[n2];
 for (Int_t i2=0; i2<n; i2++) {
    x2[i2] = i2*binwidth;
    y2[i2] = a2+b2*x2[i2]+cc2*x2[i2]*x2[i2]+dd2*x2[i2]*x2[i2]*x2[i2];
 }
 TGraph *fitfun2 = new TGraph (n2, x2, y2);

 const int binwidth3 = SBright->GetBinWidth(1);
 const int NumMomBins3 = int(1000 / binwidth3);
 float a3=fit3->GetParameter(0);
 float b3=fit3->GetParameter(1);
 float cc3=fit3->GetParameter(2);
 float dd3=fit3->GetParameter(3);

 Int_t n3 =NumMomBins3;
 Double_t x3[n3], y3[n3];
 for (Int_t i3=0; i3<n; i3++) {
    x3[i3] = i3*binwidth;
    y3[i3] = a3+b3*x3[i3]+cc3*x3[i3]*x3[i3]+dd3*x3[i3]*x3[i3]*x3[i3];
 }
 TGraph *fitfun3 = new TGraph (n3, x3, y3);

 //TGraph* fitfun;


  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// CATS input
//  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
//  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
//  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
//  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
//  CATSinput->ReadResFile();
//  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
//  CATSinput->ReadSigmaFile();
//  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
//                                       suffix.Data());
//  CATSinput->ObtainCFs(10, 250, 400);

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  /// Set up the CATS ranges, lambda parameters, etc.

  // pp radius
//  const double ppRadius = 1.249;

//  CATS AB_pXim;
//  tidy->GetCatsProtonXiMinus(&AB_pXim, 41, -4.9, 400, TidyCats::sGaussian,
//                             TidyCats::pHALQCD, 12);
//  AB_pXim.SetAnaSource(0, ppRadius);
//  AB_pXim.KillTheCat();
//  DLM_Ck* Ck_pXim = new DLM_Ck(1, 0, AB_pXim);
//  Ck_pXim->SetSourcePar(0, ppRadius);
//  Ck_pXim->Update();

//  CATS AB_pXim1530;
//  tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, 41, -4.9, 400,
//                                 TidyCats::sGaussian);
//  AB_pXim1530.SetAnaSource(0, ppRadius);
//  AB_pXim1530.KillTheCat();
//  DLM_Ck* Ck_pXim1530 = new DLM_Ck(1, 0, AB_pXim1530);
//  Ck_pXim1530->SetSourcePar(0, ppRadius);
//  Ck_pXim1530->Update();

//  CATS AB_pSigma0;
//  tidy->GetCatsProtonSigma0(&AB_pSigma0, 21, -9.99, 400, TidyCats::sGaussian,
//                            TidyCats::pSigma0ESC16);
//  AB_pSigma0.KillTheCat();
//  DLM_Ck* Ck_pSigma0;
//  Ck_pSigma0 = new DLM_Ck(1, 0, AB_pSigma0);
//  Ck_pSigma0->SetSourcePar(0, ppRadius);
//  Ck_pSigma0->Update();


//  CATS AB_pL;
//  tidy->GetCatsProtonLambda(&AB_pL, 41, -9.99, 400, TidyCats::sGaussian,
//                            TidyCats::pNLOWF);
//  AB_pL.KillTheCat();
//  DLM_Ck* Ck_pL = new DLM_Ck(1, 0, AB_pL);
//  Ck_pL->SetSourcePar(0, ppRadius);
//  Ck_pL->Update();

//  auto grLambdaGenuine = new TGraph();
//  FillWaveGraph(AB_pL, grLambdaGenuine);


//  // Proton lambda parameters
//  const double protonPurity = 0.9943;
//  const double protonPrimary = 0.823;
//  const double protonLambda = 0.125;
//  const double protonSecondary = protonLambda / (1. - protonPrimary);
//  const Particle proton(
//      protonPurity,
//      protonPrimary,
//      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
//          * (1 - protonSecondary) } });

//  // Lambda lambda parameters
//  // some parameters, stolen from PLOnly.C
//  const double lambdaPurity = 0.946;
//  const double lambdaPrimaryAndSigma0 = 0.76; // includes primary Lambda and feed-down from Sigma0
//  double PrimLambda = 3. / 4. * lambdaPrimaryAndSigma0; // isospin 3 : 1
//  double SecSigLambda = 1. / 4. * lambdaPrimaryAndSigma0;  // decay probability = 100%!
//  double SecXiLambda = (1. - lambdaPrimaryAndSigma0) / 2.; // same for Xi0 and Xim
//  const Particle lambda(lambdaPurity, PrimLambda, { SecSigLambda, SecXiLambda,
//                            SecXiLambda });

//  // Xi lambda parameters
//  // some parameters, stolen from PLOnly.C
//  const double xiPurity = 0.915;
//  const double Xi01530XimProdFraction = 1 / 2.;
//  const double Xim1530XimProdFraction = 1 / 2.;
//  const double Xi01530Xim_BR = 2 / 3.;
//  const double Xim1530Xim_BR = 1 / 3.;
//  const double OmegamXimProdFraction = 1 / 10.;
//  const double OmegamXim_BR = 0.086;  // Value given by PDG, 8.6 pm 0.4 %
//  double XiNormalization = 1 + OmegamXimProdFraction * OmegamXim_BR
//      + Xi01530XimProdFraction * Xi01530Xim_BR
//      + Xim1530XimProdFraction * Xim1530Xim_BR;
//  double SecOmegaXim = OmegamXimProdFraction * OmegamXim_BR
//      / (double) XiNormalization;
//  double SecXi01530Xim = Xi01530XimProdFraction * Xi01530Xim_BR
//      / (double) XiNormalization;
//  double SecXim1530Xim = Xim1530XimProdFraction * Xim1530Xim_BR
//      / (double) XiNormalization;
//  double PrimXim = 1. / (double) XiNormalization;
//  const Particle xi(xiPurity, PrimXim, { SecOmegaXim, SecXi01530Xim,
//                        SecXim1530Xim });

//  // Sigma0 lambda parameters
//  const Particle sigma0(0.274, 1., { { 0 } });


//  const CATSLambdaParam lambdaParamPL(proton, lambda);
//  const CATSLambdaParam lambdaParamPXi(proton, xi);
//  const CATSLambdaParam lambdaParamPSigma0(proton, sigma0);

//  const double lam_pL = lambdaParamPL.GetLambdaParam(CATSLambdaParam::Primary);
//  const double lam_pL_fake = lambdaParamPL.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::Fake)
//      + lambdaParamPL.GetLambdaParam(CATSLambdaParam::Fake,
//                                   CATSLambdaParam::Primary)
//      + lambdaParamPL.GetLambdaParam(CATSLambdaParam::Fake,
//                                   CATSLambdaParam::Fake);
//  const double lam_pL_pS0 = lambdaParamPL.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 0);
//  const double lam_pL_pXm = lambdaParamPL.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 1);

//  const double lam_pXim = lambdaParamPXi.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::Primary);
//  const double lam_pXim_pXim1530 = lambdaParamPXi.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::FeedDown, 0, 2);
//  const double lam_pXim_fake = lambdaParamPXi.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::Fake)
//      + lambdaParamPXi.GetLambdaParam(CATSLambdaParam::Fake,
//                                      CATSLambdaParam::Primary)
//      + lambdaParamPXi.GetLambdaParam(CATSLambdaParam::Fake,
//                                      CATSLambdaParam::Fake);

//  float lam_pS0_SB = lambdaParamPSigma0.GetLambdaParam(
//      CATSLambdaParam::Primary, CATSLambdaParam::Fake, 0, 0);
//  lam_pS0_SB += lambdaParamPSigma0.GetLambdaParam(
//      CATSLambdaParam::FeedDown, CATSLambdaParam::Fake, 0, 0);
//  lam_pS0_SB += lambdaParamPSigma0.GetLambdaParam(
//      CATSLambdaParam::FeedDown, CATSLambdaParam::Fake, 1, 0);
//  lam_pS0_SB += lambdaParamPSigma0.GetLambdaParam(
//      CATSLambdaParam::Fake, CATSLambdaParam::Fake);

//  DLM_CkDecomposition CkDec_pL("pLambda", 4, *Ck_pL, nullptr);  // we apply the momentum smearing later on, doesn't make a difference
//  DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0, nullptr);
//  DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim, nullptr);
//  DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530, nullptr);

//  CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
//                             DLM_CkDecomposition::cFeedDown, &CkDec_pXim1530,
//                             CATSinput->GetResFile(3));  //from Xi-(1530)
//  CkDec_pXim.AddContribution(1,
//                             1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
//                             DLM_CkDecomposition::cFeedDown);  //other feed-down (flat)
//  CkDec_pXim.AddContribution(2, lam_pXim_fake, DLM_CkDecomposition::cFake);
//  CkDec_pXim.Update();

//  CkDec_pL.AddContribution(0, lam_pL_pS0, DLM_CkDecomposition::cFeedDown,
//                           &CkDec_pSigma0, CATSinput->GetResFile(1));
//  CkDec_pL.AddContribution(1, lam_pL_pXm, DLM_CkDecomposition::cFeedDown,
//                           &CkDec_pXim, CATSinput->GetResFile(2));
//  CkDec_pL.AddContribution(2,
//                           1. - lam_pL - lam_pL_pS0 - lam_pL_pXm - lam_pL_fake,
//                           DLM_CkDecomposition::cFeedDown);
//  CkDec_pL.AddContribution(3, lam_pL_fake, DLM_CkDecomposition::cFake);
//  CkDec_pL.Update();

//  auto grLambdaDecomposition = new TGraph();
//  for (int i = 0; i < Ck_pL->GetNbins(); ++i) {
//    const double mom = Ck_pL->GetBinCenter(0, i);
//    grLambdaDecomposition->SetPoint(i, mom, CkDec_pL.EvalCk(mom));
//  }

//  std::cout << "Primary fraction for the p-L CF of " << lam_pL << "\n";
//  std::cout << "Sideband contribution to p-Sigma0 " << lam_pS0_SB << "\n";

  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // apply the smearing with the photon, and the momentum smearing

//  auto fitfunSmearedSigma = GetSmearedCF(GetSmearedCF(grLambdaGenuine, histSmear),
//                                 CATSinput->GetSigmaFile(1));
//  auto grDecompSmearedSigma = GetSmearedCF(GetSmearedCF(grLambdaDecomposition, histSmear),
//                                       CATSinput->GetSigmaFile(1));

//  auto fitfunSmearedSideband = GetSmearedCF(
//      GetSmearedCF(grLambdaGenuine, histSmearSideband),
//      CATSinput->GetSigmaFile(1));
//  auto grDecompSmearedSideband = GetSmearedCF(
//      GetSmearedCF(grLambdaDecomposition, histSmearSideband),
//      CATSinput->GetSigmaFile(1));

//  auto grSmearedSigmaForSigma = new TGraph();
//  auto grSmearedSidebandForSigma= new TGraph();

//  double x, y;
//  for (int i = 0; i < 1000; ++i) {
//    grDecompSmearedSigma->GetPoint(i, x, y);
//    grSmearedSigmaForSigma->SetPoint(i, x, 1 + (y - 1) * lam_pS0_SB);
//  }

//  for (int i = 0; i < 1000; ++i) {
//    grDecompSmearedSideband->GetPoint(i, x, y);
//    grSmearedSidebandForSigma->SetPoint(i, x, 1 + (y - 1) * lam_pS0_SB);
//  }

//  auto grDeviationSigma = new TGraph();
//  double x1, x2, y1, y2;
//  for (int i = 0; i < 1000; ++i) {
//    grSmearedSigmaForSigma->GetPoint(i, x1, y1);
//    grSmearedSidebandForSigma->GetPoint(i, x2, y2);
//    grDeviationSigma->SetPoint(i, x1, (y1-y2)/y1);
//  }


 auto smearingp=GetSmearedCF2(fitfun, histSmearFullPhiPeakKp);
 auto smearingm=GetSmearedCF2(fitfun, histSmearFullPhiPeakKm);
 auto smearingleftp=GetSmearedCF2(fitfun2, histSmearFullSBleftKp);
 auto smearingleftm=GetSmearedCF2(fitfun2, histSmearFullSBleftKm);
 auto smearingrightp=GetSmearedCF2(fitfun3, histSmearFullSBrightKp);
 auto smearingrightm=GetSmearedCF2(fitfun3, histSmearFullSBrightKm);


 auto smearingtot=GetSmearedCF2(fitfun, hhp1);
 auto smearinglefttot=GetSmearedCF2(fitfun2, hhl1);
 auto smearingrighttot=GetSmearedCF2(fitfun3, hhr1);

 auto smearingtotleft=GetSmearedCF2(fitfun, hhl1);
 auto smearingtotright=GetSmearedCF2(fitfun, hhr1);

 auto smearingsblp=GetSmearedCF2(fitfun2, hhp1);
 auto smearingsbrp=GetSmearedCF2(fitfun3, hhp1);

 auto smearingsblr=GetSmearedCF2(fitfun2, hhr1);
 auto smearingsbrl=GetSmearedCF2(fitfun3, hhl1);

// auto smearinglefttophi=GetSmearedCF(smearinglefttot, hhp1);
// auto smearingrighttophi=GetSmearedCF(smearingrighttot, hhp1);


 auto smearinglefttophi=GetSmearedCF(smearinglefttot, hhl1);
 auto smearingrighttophi=GetSmearedCF(smearingrighttot, hhr1);


//  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  // Plotting



//  auto i = new TCanvas("i", "i", 650, 550);
//  ratPeak->Draw("col");
//  ratPeak->GetXaxis()->SetRangeUser(0, 500);
//  ratPeak->GetYaxis()->SetRangeUser(0, 500);
//  ratPeak->GetXaxis()->SetNdivisions(505);
//  ratPeak->GetYaxis()->SetNdivisions(505);
//  ratPeak->SetTitle(
//      ";  Ratio_{p#minus K^{+}/p#minus K^{-}}; Ratio_{p#minus#phi} (MeV/#it{c})");
//  i->Print("ratioPeak.pdf");
//  i->Print("ratioPeak.png");


//  auto j = new TCanvas("j", "j", 650, 550);
//  ratlSB->Draw("col");
//  ratlSB->GetXaxis()->SetRangeUser(0, 500);
//  ratlSB->GetYaxis()->SetRangeUser(0, 500);
//  ratlSB->GetXaxis()->SetNdivisions(505);
//  ratlSB->GetYaxis()->SetNdivisions(505);
//  ratlSB->SetTitle(
//      "; Ratio_{p#minus K^{+}/p#minus K^{-}}; Ratio_{p#minus#phi} (MeV/#it{c})");
//  j->Print("ratioSBleft.pdf");
//  j->Print("ratioSBleft.png");


//  auto k = new TCanvas("k", "k", 650, 550);
//  ratrSB->Draw("col");
//  ratrSB->GetXaxis()->SetRangeUser(0, 500);
//  ratrSB->GetYaxis()->SetRangeUser(0, 500);
//  ratrSB->GetXaxis()->SetNdivisions(505);
//  ratrSB->GetYaxis()->SetNdivisions(505);
//  ratrSB->SetTitle(
//      "; Ratio_{p#minus K^{+}/p#minus K^{-}}; Ratio_{p#minus#phi} (MeV/#it{c})");
//  k->Print("ratioSBright.pdf");
//  k->Print("ratioSBright.png");


//  auto l = new TCanvas("l", "l", 650, 550);
//  totPeak->Draw("col");
//  totPeak->GetXaxis()->SetRangeUser(0, 500);
//  totPeak->GetYaxis()->SetRangeUser(0, 500);
//  totPeak->GetXaxis()->SetNdivisions(505);
//  totPeak->GetYaxis()->SetNdivisions(505);
//  totPeak->SetTitle(
//      "; #it{k}*_{p#minus K^{+}+p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
//  l->Print("totPeak.pdf");
//  l->Print("totPeak.png");

//  auto m = new TCanvas("m", "m", 650, 550);
//  totlSB->Draw("col");
//  totlSB->GetXaxis()->SetRangeUser(0, 500);
//  totlSB->GetYaxis()->SetRangeUser(0, 500);
//  totlSB->GetXaxis()->SetNdivisions(505);
//  totlSB->GetYaxis()->SetNdivisions(505);
//  totlSB->SetTitle(
//      "; #it{k}*_{p#minus K^{+}+p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
//  m->Print("totSBleft.pdf");
//  m->Print("totSBleft.png");

//  auto nn = new TCanvas("n", "n", 650, 550);
//  totrSB->Draw("col");
//  totrSB->GetXaxis()->SetRangeUser(0, 500);
//  totrSB->GetYaxis()->SetRangeUser(0, 500);
//  totrSB->GetXaxis()->SetNdivisions(505);
//  totrSB->GetYaxis()->SetNdivisions(505);
//  totrSB->SetTitle(
//      "; #it{k}*_{p#minus K^{+}+p#minus K^{-}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
//  nn->Print("totSBright.pdf");
//  nn->Print("totSBright.png");



  int lineWidth = 3;
  DreamPlot::SetStyleGraph(fitfun, 20, kGray, 0.7);
  fitfun->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingp, 20, kBlue, 0.7);
  smearingp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingm, 20, kRed, 0.7);
  smearingm->SetLineWidth(4);
// smearing->SetLineStyle(2);
// smearing->SetMarkerColor(kRed);

  auto gc3 = new TCanvas("CFpPhi_MCTRUTH", "CFpPhi_MCTRUTH", 0, 0, 650, 550);
  gc3->SetRightMargin(right);
  gc3->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.98, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  phitruth->Draw("same");  
  smearingp->Draw("same");
  smearingm->Draw("same");
  //fitfun->Draw("L3same");

  auto gleg5 = new TLegend(0.65, 0.6, 0.65 + 0.25, 0.9);
  gleg5->SetBorderSize(0);
  gleg5->SetTextFont(42);
  gleg5->SetTextSize(gStyle->GetTextSize() * 0.7);
  gleg5->SetFillStyle(0);
  gleg5->AddEntry(phitruth,"CF MC TRUTH", "pez");
  gleg5->AddEntry(smearingp,"smeared CF K^{+}","l");
  gleg5->AddEntry(smearingm,"smeared CF K^{-}","l");
  //gleg5->AddEntry(fitfun,"Fit function","l");
  gleg5->Draw("same");
  gc3->Print("CFMCtruth.pdf");
  gc3->Print("CFMCtruth.png");

  DreamPlot::SetStyleGraph(smearingleftp, 20, kBlue, 0.5);
  smearingleftp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingleftm, 20, kRed, 0.5);
  smearingleftm->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(fitfun2, 20, kGray, 0.7);
  fitfun2->SetLineWidth(4);

  auto g2c3 = new TCanvas("CFpPhi_SBleft", "CFpPhi_SBleft", 0, 0, 650, 550);
  g2c3->SetRightMargin(right);
  g2c3->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBleft->Draw("same");  
  smearingleftp->Draw("L3same");
  smearingleftm->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto gle1g5 = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  gle1g5->SetBorderSize(0);
  gle1g5->SetTextFont(42);
  gle1g5->SetTextSize(gStyle->GetTextSize() * 0.7);
  gle1g5->SetFillStyle(0);
  gle1g5->AddEntry(SBleft,"SB left", "pez");
  gle1g5->AddEntry(smearingleftp,"smeared SB left K^{+}","l");
  gle1g5->AddEntry(smearingleftm, "smeared SB left K^{-}","l");
  //gle1g5->AddEntry(fitfun2,"Fit function", "l");
  gle1g5->Draw("same");
  g2c3->Print("CFSBleft.pdf");
  g2c3->Print("CFSBleft.png");


  DreamPlot::SetStyleGraph(smearingrightp, 20, kBlue, 0.5);
  smearingrightp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingrightm, 20, kRed, 0.5);
  smearingrightm->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(fitfun3, 20, kGray, 0.7);
  fitfun3->SetLineWidth(4);



  auto daw2 = new TCanvas("CFpPhi_SBright", "CFpPhi_SBright", 0, 0, 650, 550);
  daw2->SetRightMargin(right);
  daw2->SetTopMargin(top);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.95, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBright->Draw("same");
  smearingrightp->Draw("L3same");
  smearingrightm->Draw("L3same");
  //fitfun3->Draw("L3same");
  auto legdaw = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  legdaw->SetBorderSize(0);
  legdaw->SetTextFont(42);
  legdaw->SetTextSize(gStyle->GetTextSize() * 0.7);
  legdaw->SetFillStyle(0);
//  legdaw->SetHeader("p#minus#kern[-0.95]{ }#Lambda #chiEFT (NLO)");
  legdaw->AddEntry(SBright,"SB right", "l");
  legdaw->AddEntry(smearingrightp,"smeared SB right K^{+}","l");
  legdaw->AddEntry(smearingrightm,"smeared SB right K^{-}","l");
  //legdaw->AddEntry( fitfun3,"Fit function","l");
  legdaw->Draw("same");
  daw2->Print("CFSBright.pdf");
  daw2->Print("CFSBright.png");



  DreamPlot::SetStyleGraph(smearingleftp, 20, kBlue, 0.5);
  smearingleftp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingrightp, 20, kRed, 0.5);
  smearingrightp->SetLineWidth(lineWidth);


  auto g2c36 = new TCanvas("pPhi_SBp", "pPhi_SBp", 0, 0, 650, 550);
  g2c36->SetRightMargin(right);
  g2c36->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  //SBleft->Draw("same");  
  smearingleftp->Draw("L3same");
  smearingrightp->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto gle1g56 = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  gle1g56->SetBorderSize(0);
  gle1g56->SetTextFont(42);
  gle1g56->SetTextSize(gStyle->GetTextSize() * 0.7);
  gle1g56->SetFillStyle(0);
  //gle1g56->AddEntry(SBleft,"SB left", "pez");
  gle1g56->AddEntry(smearingleftp,"smeared SB left K^{+}","l");
  gle1g56->AddEntry(smearingrightp, "smeared SB right K^{+}","l");
  //gle1g56->AddEntry(fitfun2,"Fit function", "l");
  gle1g56->Draw("same");
  g2c36->Print("pPhi_SBp.pdf");
  g2c36->Print("pPhi_SBp.png");


  DreamPlot::SetStyleGraph(smearingleftm, 20, kBlue, 0.5);
  smearingleftm->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingrightm, 20, kRed, 0.5);
  smearingrightm->SetLineWidth(lineWidth);


  auto g2c33 = new TCanvas("pPhi_SBm", "pPhi_SBm", 0, 0, 650, 550);
  g2c33->SetRightMargin(right);
  g2c33->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  //SBleft->Draw("same");  
  smearingleftm->Draw("L3same");
  smearingrightm->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto gle1g55 = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  gle1g55->SetBorderSize(0);
  gle1g55->SetTextFont(42);
  gle1g55->SetTextSize(gStyle->GetTextSize() * 0.7);
  gle1g55->SetFillStyle(0);
  //gle1g55->AddEntry(SBleft,"SB left", "pez");
  gle1g55->AddEntry(smearingleftm,"smeared SB left K^{-}","l");
  gle1g55->AddEntry(smearingrightm, "smeared SB right K^{-}","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  gle1g55->Draw("same");
  g2c33->Print("pPhi_SBm.pdf");
  g2c33->Print("pPhi_SBm.png");


  DreamPlot::SetStyleGraph(smearingtotleft, 20, kBlue, 0.5);
  smearingtotleft->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingtotright, 20, kRed, 0.5);
  smearingtotright->SetLineWidth(lineWidth);

  auto t1 = new TCanvas("pPhi_TruthSB", "pPhi_TruthSB", 0, 0, 650, 550);
  t1->SetRightMargin(right);
  t1->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.98, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  phitruth->Draw("same");
  smearingtotleft->Draw("L3same");
  smearingtotright->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto t1leg = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  t1leg->SetBorderSize(0);
  t1leg->SetTextFont(42);
  t1leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  t1leg->SetFillStyle(0);
  t1leg->AddEntry(SBleft,"MC truth", "pez");
  t1leg->AddEntry(smearingtotleft,"CF smearing Matrix left SB","l");
  t1leg->AddEntry(smearingtotright,"CF right smearing Matrix right SB","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  t1leg->Draw("same");
  t1->Print("MCtruthSBsmearing.pdf");
  t1->Print("MCtruthSBsmearing.png");


  DreamPlot::SetStyleGraph(smearingtotleft, 20, kBlue, 0.5);
  smearingtotleft->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingtotright, 20, kRed, 0.5);
  smearingtotright->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingtot, 20, kViolet, 0.5);
  smearingtot->SetLineWidth(lineWidth);


  auto tt22 = new TCanvas("pPhi_TruthALL", "pPhi_TruthALL", 0, 0, 650, 550);
  tt22->SetRightMargin(right);
  tt22->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(1.01, 1.07);
  dummyHist->GetXaxis()->SetNdivisions(504);
  phitruth->Draw("same");
  smearingtotleft->Draw("L3same");
  smearingtotright->Draw("L3same");
  smearingtot->Draw("L3same");

  //fitfun2->Draw("L3same");
  auto tt22leg = new TLegend(0.2, 0.2, 0.2 + 0.25, 0.45);
  tt22leg->SetBorderSize(0);
  tt22leg->SetTextFont(42);
  tt22leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  tt22leg->SetFillStyle(0);
  tt22leg->AddEntry(SBleft,"MC truth", "pez");
  tt22leg->AddEntry(smearingtotleft,"CF smearing Matrix left SB","l");
  tt22leg->AddEntry(smearingtotright,"CF right smearing Matrix right SB","l");
  tt22leg->AddEntry(smearingtot,"CF smearing Matrix peak","l");

  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  tt22leg->Draw("same");
  tt22->Print("MCtruthSBsmearingALL.pdf");
  tt22->Print("MCtruthSBsmearingALL.png");



  DreamPlot::SetStyleGraph(smearingtot, 20, kViolet, 0.5);
  smearingtot->SetLineWidth(lineWidth);


  auto t2 = new TCanvas("pPhi_PeakTOT", "pPhi_PeakTOT", 0, 0, 650, 550);
  t2->SetRightMargin(right);
  t2->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(1.01, 1.07);
  dummyHist->GetXaxis()->SetNdivisions(504);
  phitruth->Draw("same");
  smearingtot->Draw("L3same");
  //smearingtotright->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto t2leg = new TLegend(0.2, 0.2, 0.2 + 0.25, 0.45);
  t2leg->SetBorderSize(0);
  t2leg->SetTextFont(42);
  t2leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  t2leg->SetFillStyle(0);
  t2leg->AddEntry(phitruth,"MC TRUTH", "pez");
  t2leg->AddEntry(smearingtot,"CF smearing Matrix peak","l");
  //t2leg->AddEntry(smearingtotright,"CF right smearing Matrix","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  t2leg->Draw("same");
  t2->Print("smearedCFpeak.pdf");
  t2->Print("smearedCFpeak.png");

  DreamPlot::SetStyleGraph(smearinglefttot, 20, kBlue, 0.5);
  smearingtot->SetLineWidth(lineWidth);

  auto t3 = new TCanvas("pPhi_SBleftTOT", "pPhi_SBleftTOT", 0, 0, 650, 550);
  t3->SetRightMargin(right);
  t3->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBleft->Draw("same");
  smearinglefttot->Draw("L3same");
  //smearingtotright->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto t3leg = new TLegend(0.45, 0.2, 0.45 + 0.25, 0.45);
  t3leg->SetBorderSize(0);
  t3leg->SetTextFont(42);
  t3leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  t3leg->SetFillStyle(0);
  t3leg->AddEntry(SBleft,"SB left", "pez");
  t3leg->AddEntry(smearinglefttot,"CF smearing Matrix left SB","l");
  //t3leg->AddEntry(smearingtotright,"CF right smearing Matrix","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  t3leg->Draw("same");
  t3->Print("smearedCFleftSB.pdf");
  t3->Print("smearedCFleftSB.png");

  DreamPlot::SetStyleGraph(smearingrighttot, 20, kRed, 0.5);
  smearingrighttot->SetLineWidth(lineWidth);

  auto t4 = new TCanvas("pPhi_SBrightTOT", "pPhi_SBrightTOT", 0, 0, 650, 550);
  t4->SetRightMargin(right);
  t4->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.95, 1.06);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBright->Draw("same");
  smearingrighttot->Draw("L3same");
  //smearingtotright->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto t4leg = new TLegend(0.2, 0.2, 0.2 + 0.25, 0.45);
  t4leg->SetBorderSize(0);
  t4leg->SetTextFont(42);
  t4leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  t4leg->SetFillStyle(0);
  t4leg->AddEntry(SBright,"SB right", "pez");
  t4leg->AddEntry(smearingrighttot,"CF smearing Matrix right SB","l");
  //t4leg->AddEntry(smearingtotright,"CF right smearing Matrix","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  t4leg->Draw("same");
  t4->Print("smearedCFrightSB.pdf");
  t4->Print("smearedCFrightSB.png");

  auto t666 = new TCanvas("pPhi_SB", "pPhi_SB", 0, 0, 650, 550);
  t666->SetRightMargin(right);
  t666->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  //SBleft->Draw("same");
  smearingtotleft->Draw("L3same");
  smearingtotright->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto t666leg = new TLegend(0.2, 0.7, 0.2 + 0.25, 0.9);
  t666leg->SetBorderSize(0);
  t666leg->SetTextFont(42);
  t666leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  t666leg->SetFillStyle(0);
  //t666leg->AddEntry(SBleft,"SB left", "pez");
  t666leg->AddEntry(smearingtotleft,"smeared SB left","l");
  t666leg->AddEntry(smearingtotright, "smeared SB right","l");
  //t666leg->AddEntry(fitfun2,"Fit function", "l");
  t666leg->Draw("same");
  t666->Print("pPhi_TOT.pdf");
  t666->Print("pPhi_TOT.png");


  DreamPlot::SetStyleGraph(smearinglefttot, 20, kBlue, 0.5);
  smearinglefttot->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingsblp, 20, kViolet, 0.5);
  smearingsblp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingsblr, 20, kRed, 0.5);
  smearingsblr->SetLineWidth(lineWidth);


  auto tt223 = new TCanvas("pPhi_SBleftALL", "pPhi_SBleftALL", 0, 0, 650, 550);
  tt223->SetRightMargin(right);
  tt223->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.1);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBleft->Draw("same");
  smearinglefttot->Draw("L3same");
  smearingsblp->Draw("L3same");
  smearingsblr->Draw("L3same");

  //fitfun2->Draw("L3same");
  auto tt223leg = new TLegend(0.45, 0.2, 0.45 + 0.25, 0.45);
  tt223leg->SetBorderSize(0);
  tt223leg->SetTextFont(42);
  tt223leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  tt223leg->SetFillStyle(0);
  tt223leg->AddEntry(SBleft,"SB left", "pez");
  tt223leg->AddEntry(smearinglefttot,"CF smearing Matrix left SB","l");
  tt223leg->AddEntry(smearingsblr,"CF smearing Matrix right SB","l");
  tt223leg->AddEntry(smearingsblp,"CF smearing Matrix peak","l");

  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  tt223leg->Draw("same");
  tt223->Print("MCSBleftsmearingALL.pdf");
  tt223->Print("MCSBleftsmearingALL.png");


  DreamPlot::SetStyleGraph(smearingrighttot, 20, kRed, 0.5);
  smearingrighttot->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingsbrp, 20, kViolet, 0.5);
  smearingsbrp->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingsbrl, 20, kBlue , 0.5);
  smearingsbrl->SetLineWidth(lineWidth);


  auto tt77 = new TCanvas("pPhi_SBrightALL", "pPhi_SBrightALL", 0, 0, 650, 550);
  tt77->SetRightMargin(right);
  tt77->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.95, 1.06);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBright->Draw("same");
  smearingrighttot->Draw("L3same");
  smearingsbrp->Draw("L3same");
  smearingsbrl->Draw("L3same");

  //fitfun2->Draw("L3same");
  auto tt77leg = new TLegend(0.2, 0.2, 0.2 + 0.25, 0.45);
  tt77leg->SetBorderSize(0);
  tt77leg->SetTextFont(42);
  tt77leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  tt77leg->SetFillStyle(0);
  tt77leg->AddEntry(SBright,"SB right", "pez");
  tt77leg->AddEntry(smearingrighttot,"CF smearing Matrix right SB","l");
  tt77leg->AddEntry(smearingsbrl,"CF smearing Matrix left SB","l");
  tt77leg->AddEntry(smearingsbrp,"CF smearing Matrix peak","l");

  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  tt77leg->Draw("same");
  tt77->Print("MCSBrightsmearingALL.pdf");
  tt77->Print("MCSBrightsmearingALL.png");



  DreamPlot::SetStyleGraph(smearingrighttophi, 20, kRed, 0.5);
  smearingrighttophi->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearinglefttophi, 20, kBlue , 0.5);
  smearinglefttophi->SetLineWidth(lineWidth);

  auto tot123 = new TCanvas("TOTsmearing", "TOTsmearing", 0, 0, 650, 550);
  tot123->SetRightMargin(right);
  tot123->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.9, 1.16);
  dummyHist->GetXaxis()->SetNdivisions(504);
  Peak->Draw("same");
  smearingrighttophi->Draw("L3same");
  smearinglefttophi->Draw("L3same");
  //smearingsbrl->Draw("L3same");

  //fitfun2->Draw("L3same");
  auto tot123leg = new TLegend(0.55, 0.65, 0.55 + 0.25, 0.85);
  tot123leg->SetBorderSize(0);
  tot123leg->SetTextFont(42);
  tot123leg->SetTextSize(gStyle->GetTextSize() * 0.7);
  tot123leg->SetFillStyle(0);
  tot123leg->AddEntry(Peak,"CF p-#phi", "pez");
  tot123leg->AddEntry(smearingrighttophi,"p-(K^{+}K^{-})#rightarrow p-#phi RIGHT SB","l");
  tot123leg->AddEntry(smearinglefttophi,"p-(K^{+}K^{-}#rightarrow p-#phi LEFT SB","l");

  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  tot123leg->Draw("same");
  tot123->Print("TOTsmearing.pdf");
  tot123->Print("TOTsmearing.png");


  DreamPlot::SetStyleGraph(smearinglefttot, 20, kBlue, 0.5);
  smearingtot->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearinglefttophi, 20, kBlue+2, 0.5);
  smearinglefttophi->SetLineWidth(lineWidth);


  auto totleft = new TCanvas("TOTsmearingleftSB", "TOTsmearingleftSB", 0, 0, 650, 550);
  totleft->SetRightMargin(right);
  totleft->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBleft->Draw("same");
  smearinglefttot->Draw("L3same");
  smearinglefttophi->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto totleftleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
  totleftleg->SetBorderSize(0);
  totleftleg->SetTextFont(42);
  totleftleg->SetTextSize(gStyle->GetTextSize() * 0.7);
  totleftleg->SetFillStyle(0);
  totleftleg->AddEntry(SBleft,"SB left", "pez");
  totleftleg->AddEntry(smearinglefttot,"p-(K^{+}K^{-})#rightarrow p-K ","l");
  totleftleg->AddEntry(smearinglefttophi,"p-K#rightarrow p-(K^{+}K^{-})","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  totleftleg->Draw("same");
  totleft->Print("TOTsmearingleftSB.pdf");
  totleft->Print("TOTsmearingleftSB.png");

  DreamPlot::SetStyleGraph(smearingrighttot, 20, kRed, 0.5);
  smearingrighttot->SetLineWidth(lineWidth);
  DreamPlot::SetStyleGraph(smearingrighttophi, 20, kRed+2, 0.5);
  smearingrighttophi->SetLineWidth(lineWidth);

  auto totlright = new TCanvas("TOTsmearingrightSB", "TOTsmearingrightSB", 0, 0, 650, 550);
  totlright->SetRightMargin(right);
  totlright->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.92, 1.02);
  dummyHist->GetXaxis()->SetNdivisions(504);
  SBright->Draw("same");
  smearingrighttot->Draw("L3same");
  smearingrighttophi->Draw("L3same");
  //fitfun2->Draw("L3same");
  auto totlrightleg = new TLegend(0.65, 0.25, 0.65 + 0.25, 0.45);
  totlrightleg->SetBorderSize(0);
  totlrightleg->SetTextFont(42);
  totlrightleg->SetTextSize(gStyle->GetTextSize() * 0.7);
  totlrightleg->SetFillStyle(0);
  totlrightleg->AddEntry(SBright,"SB right", "pez");
  totlrightleg->AddEntry(smearingrighttot,"p-(K^{+}K^{-})#rightarrow p-K ","l");
  totlrightleg->AddEntry(smearingrighttophi,"p-K#rightarrow p-(K^{+}K^{-})","l");
  //gle1g55->AddEntry(fitfun2,"Fit function", "l");
  totlrightleg->Draw("same");
  totlright->Print("TOTsmearingrightSB.pdf");
  totlright->Print("TOTsmearingrightSB.png");




//  auto histSmear = (TH2D*) infile->Get("histSmear");
//  auto histSmearFullSBrightKp = (TH2D*) infile->Get("histSmearFullSBrightKp");
//  auto histSmearFullSBrightKm = (TH2D*) infile->Get("histSmearFullSBrightKm");
//  auto histSmearFullSBleftKp = (TH2D*) infile->Get("histSmearFullSBleftKp");
//  auto histSmearFullSBleftKm = (TH2D*) infile->Get("histSmearFullSBleftKm");
//  auto histSmearFullPhiPeakKp = (TH2D*) infile->Get("histSmearFullPhiPeakKp");
//  auto histSmearFullPhiPeakKm = (TH2D*) infile->Get("histSmearFullPhiPeakKm");

  outfile->cd();

  histSmearFullPhiPeakKp->Write("histSmearFullPhiPeakKp");
  histSmearFullPhiPeakKm->Write("histSmearFullPhiPeakKm");
  smearingp->Write("SmearedPhitruthKp");
  smearingm->Write("SmearedPhitruthKm");
  phitruth->Write("PhiTruth");
  fitfun->Write("Fitfunphitruth");

//  ratPeak->Write("ratioPeak");
//  ratlSB->Write("ratioSBleft");
//  ratrSB->Write("ratioSBright");

//  totPeak->Write("totalPeak");
//  totlSB->Write("totalSBleft");
//  totrSB->Write("totalSBright");

//  smearingtot->Write("SmearingPeakTOT");
//  smearingtotleft->Write("SmearingSBleftTOT");
//  smearingtotright->Write("SmearingSBrightTOT");

//    hp1->Write("ratioPeak");
//    hl1->Write("ratioSBleft");
//    hr1->Write("ratioSBright");

//    hhp1->Write("totalPeak");
//    hhl1->Write("totalSBleft");
//    hhr1->Write("totalSBright");

  histSmearFullSBleftKp->Write("histSmearFullSBleftKp");
  histSmearFullSBleftKm->Write("histSmearFullSBleftKm");
  smearingleftp->Write("SmearedLeftKp");
  smearingleftm->Write("SmearedLeftKm");
  SBleft->Write("SBleft");
  fitfun2->Write("FitfunSBleft");

  histSmearFullSBrightKp->Write("histSmearFullSBrightKp");
  histSmearFullSBrightKm->Write("histSmearFullSBrightKm");
  smearingrightp->Write("SmearedRightKp");
  smearingrightm->Write("SmearedRightKm");
  SBright->Write("SBright");
  fitfun3->Write("FitfunSBright");
  outfile->Close();

}
