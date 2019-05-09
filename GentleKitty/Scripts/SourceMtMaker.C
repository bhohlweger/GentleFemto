#include "TAxis.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "gsl_sf_dawson.h"
#include "DLM_Histo.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "DreamPlot.h"
#include "CATSInputSigma0.h"
#include "TRandom3.h"
#include "CATSLambdaParam.h"
#include "TDatabasePDG.h"
#include "TidyCats.h"

void boundaries(TNtuple* tuple, double radVal, TGraph *grUpper,
                TGraph *grLower) {
  tuple->Draw(Form("S >> h%.2f", radVal),
              Form("std::abs(r - %.2f) < 0.01", radVal));
  auto histRad = (TH1D*) gROOT->FindObject(Form("h%.2f", radVal));
  histRad->Sumw2();

  double median;
  double q = 0.5; // 0.5 for "median"
  histRad->ComputeIntegral(); // just a precaution
  histRad->GetQuantiles(1, &median, &q);
  auto medianBin = histRad->FindBin(median);

  auto histRadCumulative = histRad->GetCumulative();
  histRadCumulative->Scale(1 / (double) histRad->GetEntries());
//  auto canRad = new TCanvas("cRad", "cRad");
//  canRad->Divide(2, 1);
//  canRad->cd(1);
//  histRad->SetTitle(Form("r = %.2f fm; S(r); Entries", radVal));
//  histRad->Draw();
//  canRad->cd(2);
//  histRadCumulative->Draw();
//  canRad->Print(Form("Rad_%.1f.pdf", radVal));

  int binMin = 0;
  int binMax = 0;
  for (int iBin = 0; iBin < histRadCumulative->GetNbinsX(); iBin++) {
    if (binMin == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) - 0.3415)) {
      binMin = iBin;
    }
    if (binMax == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) + 0.3415)) {
      binMax = iBin;
      break;
    }
  }
  auto radMin = histRadCumulative->GetXaxis()->GetBinCenter(binMin);
  auto radMax = histRadCumulative->GetXaxis()->GetBinCenter(binMax);

//  std::cout << radVal << " " << radMin << " " << radMax << "\n";
  grUpper->SetPoint(grUpper->GetN(), radVal, radMax);
  grLower->SetPoint(grLower->GetN(), radVal, radMin);
//  delete canRad;
}

void SourceMtMaker(TH1D* histmT, int pdg1, double frac1, double massEff1,
                   double teff1, int pdg2, double frac2, double massEff2,
                   double teff2) {

  const double meanmT = histmT->GetMean();

  const int nTries = 15000;
  const double mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass()
      * 1000;
  const double mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;

  std::vector<TGraph> sourceVec;
  DLM_CleverMcLevyReso* CleverMcLevy = new DLM_CleverMcLevyReso();
  CleverMcLevy->InitNumMcIter(1000000);
  CleverMcLevy->InitStability(1, 2 - 1e-6, 2 + 1e-6);
  CleverMcLevy->InitScale(35, 0.25, 2.0);
  CleverMcLevy->InitRad(512, 0, 64);
  CleverMcLevy->InitType(2);
  CleverMcLevy->InitReso(0, 1);
  CleverMcLevy->InitReso(1, 1);
  CleverMcLevy->SetUpReso(0, 0, 1. - frac1, massEff1, teff1, mass1, massPion);
  CleverMcLevy->SetUpReso(1, 0, 1. - frac2, massEff2, teff2, mass2, massPion);
  CleverMcLevy->InitNumMcIter(1000000);

  CATS cats;
  cats.SetAnaSource(CatsSourceForwarder, CleverMcLevy, 2);
  cats.SetAnaSource(1, 2.0);

  // what we currenly use
  TGraph grSourceUpper;
  TGraph grSourceLower;
  const double rCoreLower = 0.712;
  const double rCoreUpper = 0.775;
  cats.SetAnaSource(0, rCoreLower);
  for (double i = 0; i < 150; ++i) {
    grSourceLower.SetPoint(i, i * 0.1, cats.EvaluateTheSource(0, i * 0.1, 0));
  }
  cats.SetAnaSource(0, rCoreUpper);
  for (double i = 0; i < 150; ++i) {
    grSourceUpper.SetPoint(i, i * 0.1, cats.EvaluateTheSource(0, i * 0.1, 0));
  }

  // Parametrization of r_core as a function of mT (Oton)
  auto rUpper = new TF1(
      "rUpper",
      "[&](double *x, double *p){ return p[0]*pow(x[0], p[1])+p[2]; }", 0, 5,
      3);
  auto rLower = new TF1(
      "rLower",
      "[&](double *x, double *p){ return p[0]*pow(x[0], p[1])+p[2]; }", 0, 5,
      3);
  // upper
  // 3 sigma
  rUpper->FixParameter(0, 0.750327);
  rUpper->FixParameter(1, -1.63876);
  rUpper->FixParameter(2, 0.547119);

  // 1 sigma
//  rUpper->FixParameter(0, 0.759781);
//  rUpper->FixParameter(1, -1.53513);
//  rUpper->FixParameter(2, 0.504981);

  // lower
  // 3 sigma
  rLower->FixParameter(0, 0.789205);
  rLower->FixParameter(1, -1.32559);
  rLower->FixParameter(2, 0.410767);

  // 1 sigma
//  rLower->FixParameter(0, 0.775207);
//  rLower->FixParameter(1, -1.42141);
//  rLower->FixParameter(2, 0.456818);

  TGraph grmTCore;

  double rCoremT;
  double mTRandom;
  auto file = new TFile("results.root", "RECREATE");
  auto fit = new TNtuple("fit", "fit", "r:S:mT:rCore");
  auto histR = new TH1F("histR", ";#it{r}_{Gauss,eff} (fm); Entries", 1000, 0,
                        2);

  auto gaussFit =
      new TF1(
          "gaussSource",
          [&](double *x, double *p) {return
            4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
            std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
          },
          0, 10, 1);

  for (int i = 0; i < nTries; ++i) {
    gaussFit->SetParameter(0, 1.3);
    mTRandom = histmT->GetRandom();
    if (gRandom->Uniform() > 0.5) {
      rCoremT = rUpper->Eval(mTRandom);
    } else {
      rCoremT = rLower->Eval(mTRandom);
    }
    cats.SetAnaSource(0, rCoremT);
    grmTCore.SetPoint(i, mTRandom, rCoremT);

    TGraph grSource;
    grSource.SetName(Form("source_%i", i));
    for (double i = 0; i < 100; ++i) {
      double source = cats.EvaluateTheSource(0, i * 0.1, 0);
      if(source < 0) continue;
      grSource.SetPoint(i, i * 0.1, source);
      fit->Fill(i * 0.1, source, mTRandom, rCoremT);
    }
    grSource.Fit(gaussFit, "RQN0", "", 0, 2);
    histR->Fill(gaussFit->GetParameter(0));
    sourceVec.push_back(grSource);
  }

  auto grSourceLowerConf = new TGraph();
  auto grSourceUpperConf = new TGraph();
  for (double i = 0; i < 100; ++i) {
    boundaries(fit, 0.1 * i, grSourceUpperConf, grSourceLowerConf);
  }

  auto r = new TCanvas();
  grmTCore.Draw("AP");
  DreamPlot::SetStyleGraph(&grmTCore);
  grmTCore.SetMarkerColor(kRed + 2);
  grmTCore.SetTitle("; #it{m}_{T} (GeV/#it{c}^{2}; #it{r}_{core} (fm)");
  rUpper->Draw("same");
  rUpper->SetLineColor(kGreen + 2);
  rLower->Draw("same");
  rLower->SetLineColor(kGreen + 2);
  r->Print("radiusmT.pdf");

  auto rad = new TCanvas();
  DreamPlot::SetStyleHisto(histR);
  histR->Draw();
  histR->Fit("gaus");
  rad->Print("rGauss.pdf");

  auto c = new TCanvas();
  sourceVec[0].Draw("AL");
  sourceVec[0].SetTitle(";#it{r} (fm); S(#it{r}) (fm^{-1})");
  DreamPlot::SetStyleGraph(&sourceVec[0]);
  sourceVec[0].SetLineColorAlpha(kRed + 2, 0.01);
  for (size_t i = 0; i < sourceVec.size(); ++i) {
    sourceVec[i].Draw("lsame");
    sourceVec[i].SetLineColorAlpha(kRed + 2, 0.01);
  }
  grSourceLower.SetLineColor(kBlue + 2);
  grSourceLower.SetLineStyle(2);
  grSourceUpper.SetLineColor(kBlue + 2);
  grSourceUpper.SetLineStyle(2);
  grSourceLower.Draw("same");
  grSourceUpper.Draw("same");
  grSourceLowerConf->SetLineColor(kGreen+2);
  grSourceLowerConf->SetLineStyle(3);
  grSourceLowerConf->Draw("same");
  grSourceUpperConf->SetLineColor(kGreen+2);
  grSourceUpperConf->SetLineStyle(3);
  grSourceUpperConf->Draw("same");
  grSourceUpperConf->Fit(gaussFit, "RN", "", 0, 2);
  grSourceLowerConf->Fit(gaussFit, "RN", "", 0, 2);
  c->Print(Form("source_%i_%i.pdf", pdg1, pdg2));

  file->cd();
  c->Write();
  fit->Write();
  file->Close();
}

/// =====================================================================================
int main(int argc, char* argv[]) {
  DreamPlot::SetStyle();

  const char* filename = argv[1];
  const char* trigger = argv[2];
  const char* suffix = argv[3];

  auto file1 = TFile::Open(filename);
  TDirectoryFile *dir1 = (TDirectoryFile*) (file1->FindObjectAny(
      Form("%sResults%s", trigger, suffix)));
  TList *dirResults1;
  dir1->GetObject(Form("%sResults%s", trigger, suffix), dirResults1);
  auto tmpFolder = (TList*) dirResults1->FindObject("Particle0_Particle2");
  auto histmT2D = (TH2D*) tmpFolder->FindObject("SEmTDist_Particle0_Particle2");
  auto *histmT = histmT2D->ProjectionY("proj",
                                       histmT2D->GetXaxis()->FindBin(0.),
                                       histmT2D->GetXaxis()->FindBin(0.2));

  auto d = new TCanvas();
  DreamPlot::SetStyleHisto(histmT);
  histmT->SetTitle("; #it{m}_{T} (GeV/#it{c}^{2}; Entries");
  histmT->Draw();
  d->Print("mTDist.pdf");

  SourceMtMaker(histmT, 2212, 0.3578, 1361.52, 1.65, 3212, 0.3735, 1581.73,
                4.28);
}
