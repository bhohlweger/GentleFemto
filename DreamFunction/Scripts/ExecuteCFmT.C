#include "TROOT.h"
#include "TSystem.h"

void ExecuteCFmT() {
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamDist.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamFile.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamPlot.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamData.cxx+g");
  const char* filename =
      "/Users/bernhardhohlweger/cernbox/pPb/kTmT_CF/AnalysisResults.root";
  const char* prefix = "MB";
  TString MacroNameOne = Form(
      "~/GentleFemto/DreamFunction/Scripts/GetCFvskT.C (\"%s\",\"%s\")",
      filename, prefix);

  gInterpreter->ExecuteMacro(MacroNameOne.Data());

  return;
}
