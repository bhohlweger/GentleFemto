#include "VariationAnalysis.h"
#include "DreamPlot.h"

int main(int argc, char *argv[]) {
  const char* filename = argv[1];
  VariationAnalysis* analysis = new VariationAnalysis("hCk_ReweightedppVar", 26,
                                                      81);
  analysis->ReadFitFile(filename);
  TGraphErrors* fitBand = analysis->ModelFitBands();
  TH1F* CorrelationFunction = analysis->GetCorrelationFunctio(0);
  return 0;
}
