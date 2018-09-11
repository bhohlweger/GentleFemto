#include "TROOT.h"
#include "TSystem.h"

void ExecuteCFmT(const char* filename, const char* prefix) {
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamDist.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamFile.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamKayTee.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPlot.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamData.cxx+g");

  TString MacroNameOne = Form(
      "~/GentleFemto/DreamFunction/Scripts/GetCFvskT.C (\"%s\",\"%s\")",
      filename, prefix);

  gInterpreter->ExecuteMacro(MacroNameOne.Data());

  return;
}
