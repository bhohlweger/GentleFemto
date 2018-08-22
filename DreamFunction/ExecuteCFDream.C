#include "TROOT.h"
#include "TSystem.h"

void ExecuteCFDream(const char* filename, const char* prefix) {
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamDist.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamFile.cxx+g");
  TString MacroName =
      Form(
          "/home/hohlweger/GentleFemto/DreamFunction/GetCorrelations.C (\"%s\",\"%s\")",
          filename, prefix);
  gInterpreter->ExecuteMacro(MacroName.Data());
}
