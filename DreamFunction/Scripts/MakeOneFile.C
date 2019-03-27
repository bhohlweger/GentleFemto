#include "TFile.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TClass.h"
#include <iostream>
void MakeOneFile(const char* outputFile, const char* outputName,
                 const char* inputFile) {
  TFile* outfile = TFile::Open(outputFile, "UPDATE");
  TFile* infile = TFile::Open(inputFile, "READ");
  TIter next(infile->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*) next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TGraph"))
      continue;
    TGraph *Graph = (TGraph*) key->ReadObj();
    std::cout << Graph->GetName() << std::endl;
    outfile->cd();
    Graph->Write(TString::Format("%s%s", outputName, Graph->GetName()));
  }
  infile->Close();
  outfile->Close();
  return;
}

int main(int argc, char* argv[]) {
  MakeOneFile(argv[1], argv[2], argv[3]);
  return 0;
}
