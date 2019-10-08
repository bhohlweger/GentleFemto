#include "TROOT.h"
#include "ReadDreamFile.h"
#include "DreamPlot.h"
#include "DreamKayTee.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include "TPad.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TApplication.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

gStyle->SetOptStat(0);
gStyle->SetOptTitle(kFALSE);

DreamPlot::SetStyle();

ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
TString foldername = filename;
foldername.ReplaceAll("AnalysisResults.root", "");

DreamFile->SetQuite();
DreamFile->ReadAndProjectmTHistosBBar(filename, prefix, addon, 200.);

return 0;
}
