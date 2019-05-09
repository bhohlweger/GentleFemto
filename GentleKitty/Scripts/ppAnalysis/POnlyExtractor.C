#include "VariationAnalysis.h"
#include "DreamPlot.h"
#include "TCanvas.h"

int main(int argc, char *argv[]) {
  const char* filename = argv[1];
  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar", 26,
                                                      81);
  analysis->SetRadiusRanges(100, 0.9, 1.1);
  analysis->ReadFitFile(filename);
  analysis->EvalRadius();
  TH1F* CorrelationFunction = analysis->GetCorrelationFunctio(0);
  TGraphErrors* Fit = analysis->GetModel();

  auto *c2 = new TCanvas();
  c2->cd();
  CorrelationFunction->GetXaxis()->SetRangeUser(0,250.);
  CorrelationFunction->Draw();
  Fit->Draw("PESame");
  c2->SaveAs("CF.pdf");
  return 0;
}
