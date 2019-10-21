/*
 * ForgivingReader.cxx
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */
#include "ForgivingReader.h"
#include "TFile.h"
#include <stdio.h>
#include <iostream>


ForgivingReader::ForgivingReader()
    : fInput(nullptr),
      fPrefix(""),
      fSuffix("") {
}

ForgivingReader::ForgivingReader(const char *filename, const char* prefix,
                                 const char* suffix)
    : fInput(nullptr),
      fPrefix(prefix),
      fSuffix(suffix) {
  fInput = TFile::Open(Form("%s", filename), "read");
  if (!fInput) {
    std::cout << "No input file found for " << filename << std::endl;
  }
}

ForgivingReader::~ForgivingReader() {
  // TODO Auto-generated destructor stub
}

TList* ForgivingReader::GetListInDir(const char* pathToList) {
  TList *outList = nullptr;
  if (!fInput) {
    std::cout << "No input file set! \n";
    return nullptr;
  }
  fInput->GetObject(pathToList, outList);
  if (!outList) {
    std::cerr << "Error GetListInDir: Does not exist \n";
    std::cout << pathToList << std::endl;
  }
  return outList;
}

TList* ForgivingReader::GetListInList(TList* inList,
                                      std::vector<const char*> pathToList) {
  TList* outList = inList;
  for (auto it : pathToList) {
    outList = (TList*)outList->FindObject(it);
  }
  if (outList == inList) {
    std::cerr << "Error GetListInDir: Out List equal In List \n";
  }
  return outList;
}

TH1* ForgivingReader::Get1DHistInList(TList* inList, const char* histname) {
  return (TH1*)inList->FindObject(histname);
}

TH2* ForgivingReader::Get2DHistInList(TList* inList, const char* histname) {
  return (TH2*)inList->FindObject(histname);
}
