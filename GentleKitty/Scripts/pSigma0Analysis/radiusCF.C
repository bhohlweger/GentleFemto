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
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "DreamPlot.h"
#include "CATSInputSigma0.h"
#include "DLM_CkDecomposition.h"
#include "CATSLambdaParam.h"
#include "TDatabasePDG.h"
#include "TidyCats.h"
#include "TStyle.h"
#include "TLatex.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TH2F *gr, double radius) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->Fill(mom, radius, ck->Eval(mom));
  }
}

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TH2F *gr, double radius) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->Fill(mom, radius, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
int main(int argc, char* argv[]) {
  DreamPlot::SetStyle();
  double* radius = new double[1];
  radius[0] = 0;
  int momBins = 300;
  int kmin = -1;
  int kmax = 205;

  int radBins = 325;
  const double radmin = 0.75;
  const double radOffset = 0.01;
  std::vector<double> radVec;
  for(int i = 0; i <= radBins; ++i) {
    radVec.push_back(radmin + i * radOffset);
  }
  std::cout << "Radii: " << radVec[0] << " " << radVec[radVec.size() - 1] << "\n";


  auto outfile = new TFile("radiusCF.root", "RECREATE");
  auto grLambdachiEFT = new TH2F("grLambdachiEFT",
                                 "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                 momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grLambdachiEFT);
  auto grLambdafss2 = new TH2F("grLambdafss2",
                               "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                               momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grLambdafss2);
  auto grLambdaESC16 = new TH2F("grLambdaESC16",
                                "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grLambdaESC16);
  auto grLambdaNSC97f = new TH2F("grLambdaNSC97f",
                                 "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                 momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grLambdaNSC97f);

  auto grSigma0chiEFT = new TH2F("grSigma0chiEFT",
                                 "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                 momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grSigma0chiEFT);
  auto grSigma0ESC16 = new TH2F("grSigma0ESC16",
                                "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grSigma0ESC16);
  auto grSigma0NSC97f = new TH2F("grSigma0NSC97f",
                                 "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                                 momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grSigma0NSC97f);
  auto grSigma0fss2 = new TH2F("grSigma0fss2",
                               "; #it{k}* (MeV/#it{c}); #it{r}_{0} (fm); #it{C}(#it{k}*)",
                               momBins, kmin, kmax, radBins, radVec[0], radVec[radVec.size() - 1]);
  DreamPlot::SetStyleHisto(grSigma0fss2);

  /// p-Lambda
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /// chiEFT
  auto LambdaChiEFTHaidenbauer = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                            Lednicky_SingletTriplet);
  LambdaChiEFTHaidenbauer->SetPotPar(0, 2.91);
  LambdaChiEFTHaidenbauer->SetPotPar(1, 2.78);
  LambdaChiEFTHaidenbauer->SetPotPar(2, 1.54);
  LambdaChiEFTHaidenbauer->SetPotPar(3, 2.72);

  /// fss2
  auto Lambdafss2 = new DLM_Ck(1, 4, momBins, kmin, kmax,
                               Lednicky_SingletTriplet);
  Lambdafss2->SetPotPar(0, 2.59);
  Lambdafss2->SetPotPar(1, 2.83);
  Lambdafss2->SetPotPar(2, 1.60);
  Lambdafss2->SetPotPar(3, 3.01);

  /// ESC16
  auto LambdaESC16 = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                Lednicky_SingletTriplet);
  LambdaESC16->SetPotPar(0, 1.88);
  LambdaESC16->SetPotPar(1, 3.58);
  LambdaESC16->SetPotPar(2, 1.86);
  LambdaESC16->SetPotPar(3, 3.37);

  /// NSC97f
  auto LambdaNSC97f = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                 Lednicky_SingletTriplet);
  LambdaNSC97f->SetPotPar(0, 2.51);
  LambdaNSC97f->SetPotPar(1, 3.03);
  LambdaNSC97f->SetPotPar(2, 1.75);
  LambdaNSC97f->SetPotPar(3, 3.32);

  /// p-Sigma0
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /// chiEFT
  TidyCats *tidy = new TidyCats();
  CATS chiEFTKitty;
  tidy->GetCatsProtonSigma0(&chiEFTKitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0Haidenbauer);
  chiEFTKitty.KillTheCat();
  auto Sigma0chiEFT = new DLM_Ck(1, 0, chiEFTKitty);

  /// ESC16
  CATS ESC16Kitty;
  tidy->GetCatsProtonSigma0(&ESC16Kitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0ESC16);
  ESC16Kitty.KillTheCat();
  auto Sigma0ESC16 = new DLM_Ck(1, 0, ESC16Kitty);

  /// NSC97f
  CATS NSC97fKitty;
  tidy->GetCatsProtonSigma0(&NSC97fKitty, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
  NSC97fKitty.KillTheCat();
  auto Sigma0NSC97f = new DLM_Ck(1, 0, NSC97fKitty);

  /// fss2
  DLM_Ck* Sigma0fss2 = new DLM_Ck(1, 0, momBins, kmin, kmax,
                                  Lednicky_gauss_Sigma0);

  for (int i = 0; i < radBins; ++i) {
    radius[0] = radVec[i] + 0.5 * radOffset;

    LambdaChiEFTHaidenbauer->SetSourcePar(0, radius[0]);
    LambdaChiEFTHaidenbauer->Update();
    FillCkGraph(LambdaChiEFTHaidenbauer, grLambdachiEFT, radius[0]);

    Lambdafss2->SetSourcePar(0, radius[0]);
    Lambdafss2->Update();
    FillCkGraph(Lambdafss2, grLambdafss2, radius[0]);

    LambdaESC16->SetSourcePar(0, radius[0]);
    LambdaESC16->Update();
    FillCkGraph(LambdaESC16, grLambdaESC16, radius[0]);

    LambdaNSC97f->SetSourcePar(0, radius[0]);
    LambdaNSC97f->Update();
    FillCkGraph(LambdaNSC97f, grLambdaNSC97f, radius[0]);

    Sigma0chiEFT->SetSourcePar(0, radius[0]);
    Sigma0chiEFT->Update();
    FillWaveGraph(chiEFTKitty, grSigma0chiEFT, radius[0]);

    Sigma0ESC16->SetSourcePar(0, radius[0]);
    Sigma0ESC16->Update();
    FillWaveGraph(ESC16Kitty, grSigma0ESC16, radius[0]);

    Sigma0NSC97f->SetSourcePar(0, radius[0]);
    Sigma0NSC97f->Update();
    FillWaveGraph(NSC97fKitty, grSigma0NSC97f, radius[0]);

    Sigma0fss2->SetSourcePar(0, radius[0]);
    Sigma0fss2->Update();
    FillCkGraph(Sigma0fss2, grSigma0fss2, radius[0]);
  }

  const float minCFSigma = 0.;
  const float minCFLambda = 0.;
  const float maxCFSigma = 3.;
  const float maxCFLambda = 3.;

  grLambdachiEFT->SetMinimum(minCFLambda);
  grLambdafss2->SetMinimum(minCFLambda);
  grLambdaESC16->SetMinimum(minCFLambda);
  grLambdaNSC97f->SetMinimum(minCFLambda);

  grSigma0chiEFT->SetMinimum(minCFSigma);
  grSigma0ESC16->SetMinimum(minCFSigma);
  grSigma0NSC97f->SetMinimum(minCFSigma);
  grSigma0fss2->SetMinimum(minCFSigma);

  grLambdachiEFT->SetMaximum(maxCFLambda);
  grLambdafss2->SetMaximum(maxCFLambda);
  grLambdaESC16->SetMaximum(maxCFLambda);
  grLambdaNSC97f->SetMaximum(maxCFLambda);

  grSigma0chiEFT->SetMaximum(maxCFSigma);
  grSigma0ESC16->SetMaximum(maxCFSigma);
  grSigma0NSC97f->SetMaximum(maxCFSigma);
  grSigma0fss2->SetMaximum(maxCFSigma);

  grLambdachiEFT->GetXaxis()->SetRangeUser(0, 200);
  grLambdafss2->GetXaxis()->SetRangeUser(0, 200);
  grLambdaESC16->GetXaxis()->SetRangeUser(0, 200);
  grLambdaNSC97f->GetXaxis()->SetRangeUser(0, 200);

  grSigma0chiEFT->GetXaxis()->SetRangeUser(0, 200);
  grSigma0ESC16->GetXaxis()->SetRangeUser(0, 200);
  grSigma0NSC97f->GetXaxis()->SetRangeUser(0, 200);
  grSigma0fss2->GetXaxis()->SetRangeUser(0, 200);

  grLambdachiEFT->GetZaxis()->SetTitleOffset(0.85);
  grLambdafss2->GetZaxis()->SetTitleOffset(0.85);
  grLambdaESC16->GetZaxis()->SetTitleOffset(0.85);
  grLambdaNSC97f->GetZaxis()->SetTitleOffset(0.85);

  grSigma0chiEFT->GetZaxis()->SetTitleOffset(0.85);
  grSigma0ESC16->GetZaxis()->SetTitleOffset(0.85);
  grSigma0NSC97f->GetZaxis()->SetTitleOffset(0.85);
  grSigma0fss2->GetZaxis()->SetTitleOffset(0.85);

  grLambdachiEFT->Write();
  grLambdafss2->Write();
  grLambdaESC16->Write();
  grLambdaNSC97f->Write();

  grSigma0chiEFT->Write();
  grSigma0ESC16->Write();
  grSigma0NSC97f->Write();
  grSigma0fss2->Write();

  TLatex text;
  text.SetTextSize(gStyle->GetTextSize());
  text.SetNDC(true);


  auto c = new TCanvas();
  c->SetRightMargin(0.16);
  grLambdachiEFT->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Lambda");
  text.DrawLatex(0.625, 0.82, "#chiEFT (NLO)");
  c->Print("CFrad_Lambda_chiEFT.pdf");
  c->Clear();

  grLambdafss2->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Lambda");
  text.DrawLatex(0.625, 0.82, "fss2");
  c->Print("CFrad_Lambda_fss2.pdf");
  c->Clear();

  grLambdaESC16->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Lambda");
  text.DrawLatex(0.625, 0.82, "ESC16");
  c->Print("CFrad_Lambda_ESC16.pdf");
  c->Clear();

  grLambdaNSC97f->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Lambda");
  text.DrawLatex(0.625, 0.82, "NSC97f");
  c->Print("CFrad_Lambda_NSC97f.pdf");
  c->Clear();


  grSigma0chiEFT->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Sigma^{0}");
  text.DrawLatex(0.625, 0.82, "#chiEFT (NLO)");
  c->Print("CFrad_Sigma_chiEFT.pdf");
  c->Clear();

  grSigma0ESC16->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Sigma^{0}");
  text.DrawLatex(0.625, 0.82, "ESC16");
  c->Print("CFrad_Sigma_ESC16.pdf");
  c->Clear();

  grSigma0NSC97f->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Sigma^{0}");
  text.DrawLatex(0.625, 0.82, "NSC97f");
  c->Print("CFrad_Sigma_NSC97f.pdf");
  c->Clear();

  grSigma0fss2->Draw("colz");
  text.DrawLatex(0.625, 0.88, "p#minus#kern[-0.75]{ }#Sigma^{0}");
  text.DrawLatex(0.625, 0.82, "fss2");
  c->Print("CFrad_Sigma_fss2.pdf");
  c->Clear();


  outfile->Close();
}
