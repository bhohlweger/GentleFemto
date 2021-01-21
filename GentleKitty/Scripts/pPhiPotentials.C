#include "TApplication.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"

#include "TidyCats.h"
#include "DreamPlot.h"
#include "CATSLambdaParam.h"
#include "CATSInput.h"
#include "TStyle.h"
#include "TLatex.h"

#include "DLM_CkDecomposition.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}

/// =====================================================================================
void FillPotentialGraph(CATS& cats, TGraph* gr) {
  for (double i=1; i < 100; ++i) {
    double radius = i / 10.;
    gr->SetPoint(i-1, radius, cats.EvaluateThePotential(0, 0, 10, radius));
  }
}

/// =====================================================================================
int main(int argc, char *argv[]) {

  const double rpphi = 1.08;
  const double kmin = 0;
  const double kmax = 300;
  const int nBins = 30;

  auto outfile = new TFile("pphiPotentials.root", "RECREATE");
  
  auto grPot1 = new TGraph();
  DreamPlot::SetStyleGraph(grPot1, 20, kBlue);
  grPot1->SetTitle(";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  auto grPot2 = (TGraph*)grPot1->Clone("grPot2");
  auto grPot3 = (TGraph*)grPot1->Clone("grPot3");
  auto grPot4 = (TGraph*)grPot1->Clone("grPot4");

  auto grPotential1 = (TGraph*)grPot1->Clone("grPotential1");
  grPotential1->SetTitle(";#it{r} (fm); #it{U}(#it{r})");
  auto grPotential2 = (TGraph*)grPotential1->Clone("grPotential2");
  auto grPotential3 = (TGraph*)grPotential1->Clone("grPotential3");
  auto grPotential4 = (TGraph*)grPotential1->Clone("grPotential4");
  
  CATS cats1, cats2;
  TidyCats* tidy = new TidyCats();
  tidy->GetCatsProtonPhi(&cats1, nBins, kmin, kmax, TidyCats::pYukawa, TidyCats::sGaussian);
  cats1.SetAnaSource(0, rpphi);
  cats1.KillTheCat();
  auto DLM_pphi1 = new DLM_Ck(1, 4, cats1);
  DLM_pphi1->Update();

  FillCkGraph(DLM_pphi1, grPot1);
  FillPotentialGraph(cats1, grPotential1);
  DLM_pphi1->SetPotPar(0, -0.5);
  DLM_pphi1->SetPotPar(1, 100);
  DLM_pphi1->Update();
  FillCkGraph(DLM_pphi1, grPot2);
  FillPotentialGraph(cats1, grPotential2);
  
  tidy->GetCatsProtonPhi(&cats2, nBins, kmin, kmax, TidyCats::pGaussian, TidyCats::sGaussian);
  cats2.SetAnaSource(0, rpphi);
  cats2.KillTheCat();
  auto DLM_pphi2 = new DLM_Ck(1, 4, cats2);
  DLM_pphi2->Update();

  FillCkGraph(DLM_pphi2, grPot3);
  FillPotentialGraph(cats2, grPotential3);

  DLM_pphi2->SetPotPar(0, 300);
  DLM_pphi2->SetPotPar(1, 10);
  DLM_pphi2->Update();
  FillCkGraph(DLM_pphi2, grPot4);
  FillPotentialGraph(cats2, grPotential4);
  
  outfile->cd();
  grPot1->Write("CkYukawa1");
  grPotential1->Write("Yukawa1");
  grPot2->Write("CkYukawa2");
  grPotential2->Write("Yukawa2");
  grPot3->Write("CkGaussian1");
  grPotential3->Write("Gaussian1");
  grPot4->Write("CkGaussian2");
  grPotential4->Write("Gaussian2");
  outfile->Close();
}
