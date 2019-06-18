/*
 * VariationmTAnalysis.cxx
 *
 *  Created on: Jun 18, 2019
 *      Author: schmollweger
 */
#include "VariationmTAnalysis.h"
#include "TSystemDirectory.h"
#include <iostream>

VariationmTAnalysis::VariationmTAnalysis()
    : fAnalysis(),
      fSystematic(),
      fHistname() {
  // TODO Auto-generated constructor stub

}

VariationmTAnalysis::~VariationmTAnalysis() {
  // TODO Auto-generated destructor stub
}

void VariationmTAnalysis::SetSystematic(const char* DataDir) {
  static int outputCounter = 0;
  TSystemDirectory *workdir = new TSystemDirectory("workdir", DataDir);
  TList *RootList = workdir->GetListOfFiles();
  RootList->Sort();
  TIter next(RootList);
  TObject* obj = nullptr;
  DreamSystematics Systematics(DreamSystematics::pp);
  Systematics.SetUpperFitRange(150);
  Systematics.SetBarlowUpperRange(150);
  while (obj = next()) {
    TString FileName = obj->GetName();
    if (FileName.Contains(".root")) {
      TH1F* histo = nullptr;
      TFile* File = TFile::Open(
          TString::Format("%s/%s", DataDir, FileName.Data()).Data(), "read");
      if (!File) {
        Warning(
            "VariationmTAnalysis::SetSystematic",
            TString::Format("File %s does not exist, exiting \n",
                            FileName.Data()));
        return;
      }
      TList* FileKeys = File->GetListOfKeys();
      TIter FileIter(FileKeys);
      TObject* FileObj;
      while (FileObj = FileIter()) {
        TString FileObjName = FileObj->GetName();
        if (FileObjName.Contains(fHistname) && FileObjName.Contains("MeV")) {
          histo = (TH1F*) (File->FindObjectAny(FileObjName.Data()))->Clone(
              TString::Format("%sClone", FileObjName.Data()));
          histo->SetDirectory(0);
          break;
        }
      }
      if (!histo) {
        Warning(
            "VariationmTAnalysis::SetSystematic",
            TString::Format("No Histogram found for %s in file %s. Exiting \n",
                            fHistname, FileName.Data()).Data());
        return;
      }
      if (FileName.Contains("Var0")) {
        //that's the default
        Systematics.SetDefaultHist(histo);
      } else {
        Systematics.SetVarHist(histo);
      }
      File->Close();
    }
  }
  Systematics.EvalSystematics();
  Systematics.WriteOutput(Form("%u",outputCounter));
  outputCounter++;
//  fSystematic.push_back(Systematics);
  return;
}

void VariationmTAnalysis::SetVariation(const char* VarDir) {
  return;
}

void VariationmTAnalysis::MakePlots() {
  auto c1 = new TCanvas("c24455454", "c256345654");
  c1->Divide(3, 2);
  int counter = 1;
  for (auto it : fSystematic) {
    DreamData *ProtonProton = new DreamData(Form("ProtonProton%i",counter));
    ProtonProton->SetUnitConversionData(1);
    ProtonProton->SetUnitConversionCATS(1);
    ProtonProton->SetCorrelationFunction(it.GetDefault());
    ProtonProton->SetSystematics(it.GetSystematicError(), 2);
    ProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p}", "fpe");
    ProtonProton->SetLegendName("Coulomb + Argonne #nu_{18} (fit)", "l");
    ProtonProton->SetRangePlotting(0, 200, 0.725, 3.5);
    ProtonProton->SetNDivisions(505);
    ProtonProton->SetLegendCoordinates(
        0.30, 0.65 - 0.09 * ProtonProton->GetNumberOfModels(), 0.7, 0.725);
    TPad* tmp = (TPad*)c1->cd(counter);
    ProtonProton->DrawCorrelationPlot(tmp);
    counter++;
  }
  TFile* out = TFile::Open("tmp.root", "recreate");
  out->cd();
  c1->Write();
  c1->SaveAs("mTPlots.pdf");
  out->Write();
  out->Close();

}
