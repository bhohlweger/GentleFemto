#include "TROOT.h"
#include "TSystem.h"

void ExecuteCFDream(const char* filename, const char* prefix) {
  gROOT->LoadMacro("~/Plotting/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/Plotting/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/Plotting/DreamFunction/ReadDreamFile.cxx+g");
  TString MacroName=
      Form("/home/hohlweger/Plotting/DreamFunction/GetCorrelations.C (\"%s\",\"%s\")",filename,prefix);
  gInterpreter->ExecuteMacro(MacroName.Data());
}
