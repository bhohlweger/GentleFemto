#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TidyCats.h"
#include "TDatabasePDG.h"

void testCats() {
  CATS ppCats;
  TidyCats* tidy = new TidyCats();
  tidy->GetCatsProtonProton(&ppCats, 300, 0, 300, TidyCats::sResonance);
  ppCats.SetAnaSource(0, 1.1);
  ppCats.SetAnaSource(1, 2.0);
  ppCats.KillTheCat();
  TGraph* graph = new TGraph(1000);
  TGraph* cf = new TGraph(300);
  graph->SetName("theSource");
  cf->SetName("cf");
  for (int iRad = 0; iRad < 1e5; iRad++) {
    double rad = iRad*1e-4;
    graph->SetPoint(iRad, rad, ppCats.EvaluateTheSource(11, rad, 0));
  }
  for (int iKstar = 0; iKstar < 300; iKstar++) {
    cf->SetPoint(iKstar, ppCats.GetMomentum(iKstar), ppCats.GetCorrFun(iKstar));
  }
  TFile* out = TFile::Open("out.root","recreate");
  out->cd();
  graph->Write();
  cf->Write();
  out->Close();
  return;
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}
