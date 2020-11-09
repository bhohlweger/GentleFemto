#include "DLM_CkModels.h"
#include "DLM_CkDecomposition.h"
#include "TApplication.h"
#include "TMath.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}

/// =====================================================================================
int main(int argc, char *argv[]) {
  TApplication *app = new TApplication("app", 0, 0);

  DreamPlot::SetStyle();

  const int momBins = 100;
  const double kmin = 0;
  const double kmax = 300;
  const double radius = 1.25;
  const double scatLen = atof(argv[1]);
  const double effRange = atof(argv[2]);

  std::cout << "1/f0: " << scatLen << " fm\n";
  std::cout << "d_0 : " << effRange << " fm\n";

  const int lineWidth = 3;

  auto grCkCATS = new TGraph();
  DreamPlot::SetStyleGraph(grCkCATS, 20, kBlue + 3, 0.8);
  grCkCATS->SetLineWidth(lineWidth);

  auto CkCATS = new DLM_Ck(1, 4, momBins, kmin, kmax, Lednicky_Singlet_InvScatLen);
  CkCATS->SetPotPar(0, scatLen);
  CkCATS->SetPotPar(1, effRange);
  CkCATS->SetSourcePar(0, radius);
  CkCATS->Update();
  FillCkGraph(CkCATS, grCkCATS);


  auto dummyHist = new TH1F("dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, kmax);
  DreamPlot::SetStyleHisto(dummyHist, 20, kWhite);
  dummyHist->GetYaxis()->SetRangeUser(0., 3);
  dummyHist->GetXaxis()->SetNdivisions(504);
  dummyHist->GetYaxis()->SetNdivisions(504);

  auto c = new TCanvas("CFComp1", "CFComp1"); //, 0, 0, 650, 550);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.025);
  dummyHist->Draw();
  grCkCATS->Draw("L3");

  TLatex text2;
  text2.SetNDC(true);
  text2.SetTextSize(gStyle->GetTextSize() * 0.9);
  text2.DrawLatex(0.75, 0.85, TString::Format("f_{0}^{-1} = %.2f fm^{-1}", scatLen));
  text2.DrawLatex(0.75, 0.78, TString::Format("d_{0}  = %.2f fm", effRange));

  c->Print("Lednicky.pdf");

  app->Run();

}
