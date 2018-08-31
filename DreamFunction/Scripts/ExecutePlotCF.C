#include "TROOT.h"
#include "TSystem.h"

void ExecutePlotCF() {
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamDist.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPair.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamCF.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/ReadDreamFile.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamData.cxx+g");
  gROOT->LoadMacro("~/GentleFemto/DreamFunction/DreamPlot.cxx+g");
  const char* Data_pPb = "/Users/bernhardhohlweger/cernbox/pPb/DataFullest/";
  const char* Syst_pPb = "/Users/bernhardhohlweger/cernbox/pPb/";
  const char* Sim_pPb = "";
  const char* Cats_pPb = "/Users/bernhardhohlweger/cernbox/pPb/CATS_FITS/ReweightedMEFit/20MeV/";

  const char* Data_ppHM = "/Users/bernhardhohlweger/cernbox/HM13TeV/AnalysisData/merge2046_54_55/";
  const char* Syst_ppHM = "/Users/bernhardhohlweger/cernbox/pPb/";
  const char* Sim_ppHM = "";
  const char* Cats_ppHM = "/Users/bernhardhohlweger/cernbox/HM13TeV/CATS_OUTPUT/20MeV/";

//  TString MacroNamepPb =
//      Form(
//          "~/GentleFemto/DreamFunction/Scripts/PlotCF_pPb.C (\"%s\",\"%s\",\"%s\",\"%s\")",
//          Data_pPb,Syst_pPb,Sim_pPb,Cats_pPb);
//  gInterpreter->ExecuteMacro(MacroNamepPb.Data());

  TString MacroNamepp =
      Form(
          "~/GentleFemto/DreamFunction/Scripts/PlotCF_ppHM.C (\"%s\",\"%s\",\"%s\",\"%s\")",
          Data_ppHM,Syst_ppHM,Sim_ppHM,Cats_ppHM);
  gInterpreter->ExecuteMacro(MacroNamepp.Data());
}
