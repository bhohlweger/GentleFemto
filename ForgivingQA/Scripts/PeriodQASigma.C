#include "TSystemDirectory.h"
#include "TList.h"
#include "MakeHistosGreat.h"
#include "ForgivingReader.h"
#include "PeriodQA.h"
#include "DecayQA.h"
#include <iostream>

// arguments
// 1: Path to all analysis files for the periods (named: AnalysisResults_LHC16d.root, ...)
// 2: prefix (MB/HM/...)
// 3: addon (0/1/...)
int main(int argc, char* argv[]) {
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  MakeHistosGreat::SetStyle(false);

  DrawStyle styler;
  styler.drawMarker = 2;
  styler.drawColor = 8;
  styler.drawSize = 1.1;
  styler.drawLineColor = kTeal + 3;
  styler.drawSignalFitColor = kBlue + 4;
  styler.drawBackgroundFitColor = kGreen + 2;

  PeriodQA *qa = new PeriodQA();
  qa->SetDirectory(argv[1]);
  qa->SetStyler(styler);
  qa->ProcessSigmaQA(prefix, addon);
}

