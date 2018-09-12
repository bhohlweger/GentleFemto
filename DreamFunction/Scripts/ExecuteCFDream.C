#include "TROOT.h"
#include "TSystem.h"

void ExecuteCFDream(const char* filename, const char* prefix, const char* addon="") {
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamDist.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamFile.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPlot.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamData.cxx+g");
  TString MacroNameOne = Form("~/GentleFemto/DreamFunction/Scripts/GetCorrelations.C (\"%s\",\"%s\", \"%s\")", filename, prefix, addon);
  gInterpreter->ExecuteMacro(MacroNameOne.Data());

  TString foldername = filename;
  TString FileName = "AnalysisResults.root";

  foldername.Replace(foldername.First(FileName.Data()),FileName.Length(),"");
  TString MacroNameTwo = Form("~/GentleFemto/DreamFunction/Scripts/METoSEReweighting.C (\"%s\")", foldername.Data());
  gInterpreter->ExecuteMacro(MacroNameTwo.Data());
}
