/*
 * ReadDreamFile.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "ReadDreamFile.h"
#include <iostream>
ReadDreamFile::ReadDreamFile(int nPart1,int nPart2)
:fNPart1(nPart1)
,fNPart2(nPart2)
,fSE(nullptr)
,fSEMult(nullptr)
,fME(nullptr)
,fMEMult(nullptr)
{
}

ReadDreamFile::~ReadDreamFile()
{
  // TODO Auto-generated destructor stub
}


void ReadDreamFile::SetAnalysisFile(const char* PathAnalysisFile,
                                    const char* Prefix)
{
  TFile* _file0=TFile::Open(PathAnalysisFile,"READ");
  const float normleft = 200;
  const float normright = 400;

  TDirectoryFile *dirResults=
      (TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults",Prefix)));
  TList *Results;
  dirResults->GetObject(Form("%sResults",Prefix),Results);
  TList *PartList;

  fSE=new TH1F**[fNPart1];
  fSEMult=new TH2F**[fNPart1];
  fME=new TH1F**[fNPart1];
  fMEMult=new TH2F**[fNPart1];

  for (int iPart1=0;iPart1<fNPart1;++iPart1)
  {
    fSE[iPart1]=new TH1F*[fNPart2];
    fSEMult[iPart1]=new TH2F*[fNPart2];
    fME[iPart1]=new TH1F*[fNPart2];
    fMEMult[iPart1]=new TH2F*[fNPart2];

    for (int iPart2=iPart1;iPart2<fNPart2;++iPart2)
    {
      TString FolderName=Form("Particle%i_Particle%i",iPart1,iPart2);
      PartList=(TList*)Results->FindObject(FolderName.Data());
      fSE[iPart1][iPart2]=nullptr;
      fSE[iPart1][iPart2]=
          (TH1F*)PartList->FindObject(Form("SEDist_%s",FolderName.Data()));
      if (!fSE[iPart1][iPart2]) {
        std::cout << "SE Histogramm missing from " << FolderName.Data()
                        <<std::endl;
      }
      fSE[iPart1][iPart2]->Sumw2();

      fSEMult[iPart1][iPart2]=nullptr;
      fSEMult[iPart1][iPart2]=
          (TH2F*)PartList->FindObject(Form("SEMultDist_%s",FolderName.Data()));
      if (!fSEMult[iPart1][iPart2]) {
        std::cout << "SEMult Histogramm missing from " << FolderName.Data()
                        <<std::endl;
      }
      fSEMult[iPart1][iPart2]->Sumw2();

      fME[iPart1][iPart2]=nullptr;
      fME[iPart1][iPart2]=
          (TH1F*)PartList->FindObject(Form("MEDist_%s",FolderName.Data()));
      if (!fME[iPart1][iPart2]) {
        std::cout << "ME Histogramm missing from " << FolderName.Data()
                        <<std::endl;
      }
      fME[iPart1][iPart2]->Sumw2();

      fMEMult[iPart1][iPart2]=nullptr;
      fMEMult[iPart1][iPart2]=
          (TH2F*)PartList->FindObject(Form("MEMultDist_%s",FolderName.Data()));
      if (!fMEMult[iPart1][iPart2]) {
        std::cout << "ME Mult Histogramm missing from " << FolderName.Data()
                        <<std::endl;
      }
      fMEMult[iPart1][iPart2]->Sumw2();
    }
  }
}

DreamPair* ReadDreamFile::GetPairDistributions(
    int iPart1,
    int iPart2,
    const char* name)
{
  //user needs to ensure deletion
  if (iPart2 < iPart1)
  {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamPair* pair = new DreamPair(name);
  pair->SetSEDist(fSE[iPart1][iPart2]);
  pair->SetSEMultDist(fSEMult[iPart1][iPart2]);
  pair->SetMEDist(fME[iPart1][iPart2]);
  pair->SetMEMultDist(fMEMult[iPart1][iPart2]);
  return pair;
}
