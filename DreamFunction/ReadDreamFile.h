/*
 * ReadDreamFile.h
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#ifndef DREAMFUNCTION_READDREAMFILE_H_
#define DREAMFUNCTION_READDREAMFILE_H_
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include "DreamKayTee.h"
#include "DreamDist.h"
class ReadDreamFile {
 public:
  ReadDreamFile(int nPart1, int nPart2);
  virtual ~ReadDreamFile();
  void SetAnalysisFile(const char* AnalysisFile, const char* prefix, const char* Addon="");
  void ReadkTHistos(const char* AnalysisFile, const char* prefix, const char* addon = "");
  void ReadmTHistos(const char* AnalysisFile, const char* prefix, const char* addon = "");
  DreamDist* GetPairDistributions(int iPart1, int iPart2, const char* name);
  DreamKayTee* GetkTPairDistributions(int iPart1, int iPart2,int iAPart1, int iAPart2);
  DreamKayTee* GetmTPairDistributions(int iPart1, int iPart2,int iAPart1, int iAPart2);
  const int fNPart1;
  const int fNPart2;
 private:
  TH1F*** fSE;
  TH2F*** fSEMult;
  TH2F*** fSEkT;
  TH2F*** fSEmT;
  TH1F*** fME;
  TH2F*** fMEMult;
  TH2F*** fMEkT;
  TH2F*** fMEmT;
};

#endif /* DREAMFUNCTION_READDREAMFILE_H_ */
