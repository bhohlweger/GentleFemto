#include "VariationAnalysis.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "DreamData.h"
#include "TFile.h"

int main(int argc, char *argv[]) {
  const char* filename = argv[1];
  const char* SystFile = argv[2];
  TFile* systFile = TFile::Open(SystFile,"read");
  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar", 26,
                                                      81);
  analysis->SetRadiusRanges(100, 0.9, 1.1);
  analysis->ReadFitFile(filename);
  analysis->EvalRadius();
  TGraphErrors* Fit = analysis->GetModel();
  DreamData *ProtonProton = new DreamData("ProtonProton");
  ProtonProton->SetUnitConversionData(1);
  ProtonProton->SetCorrelationFunction(analysis->GetCorrelationFunctio(0));
//  ProtonProton->SetSystematics()
  return 0;
}
