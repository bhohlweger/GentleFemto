#include "TFile.h"
#include <iostream>
#include "TString.h"
void AddSystematics(const char* output, const char* inputOne,
                    const char* inputTwo) {
  //function to add the TDirectoryFiles of two files to one.
  //you might want to just use hadd.
  TFile* fileOutput = TFile::Open(output, "RECREATE");
  TFile* fileOne = TFile::Open(inputOne, "READ");
  TFile* fileTwo = TFile::Open(inputTwo, "READ");

  TList *tmpFileOne = fileOne->GetListOfKeys();
  for (int it = 0; it < tmpFileOne->GetEntries(); ++it) {
    TString name = tmpFileOne->At(it)->GetName();
    TDirectoryFile *dir =
        (TDirectoryFile*) (fileOne->FindObjectAny(name.Data()));
    if (dir) {
      TList* List = nullptr ;
      dir->GetObject(dir->GetListOfKeys()->At(0)->GetName(),List);
      fileOutput->cd();
      TDirectoryFile* DirClone = new TDirectoryFile(dir->GetName(),dir->GetName());
      DirClone->Add(List);
      DirClone->Write(dir->GetName(), TObject::kSingleKey);
    }
  }

  TList *tmpFileTwo = fileTwo->GetListOfKeys();
  for (int it = 0; it < tmpFileTwo->GetEntries(); ++it) {
    TString name = tmpFileTwo->At(it)->GetName();
    TDirectoryFile *dir =
        (TDirectoryFile*) (fileTwo->FindObjectAny(name.Data()));
    if (name.Contains("MultSelection")) {
      continue;
    }
    if (dir) {
      TList* List = nullptr ;
      dir->GetObject(dir->GetListOfKeys()->At(0)->GetName(),List);
      fileOutput->cd();
      TDirectoryFile* DirClone = new TDirectoryFile(dir->GetName(),dir->GetName());
      DirClone->Add(List);
      DirClone->Write(dir->GetName(), TObject::kSingleKey);
    }
  }
  fileOutput->Close();
}

int main(int argc, char* argv[]) {
  AddSystematics(argv[1], argv[2], argv[3]);
  return 1;
}
