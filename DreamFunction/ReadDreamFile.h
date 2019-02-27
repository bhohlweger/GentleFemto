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
#include "DreamdEtadPhi.h"

class ReadDreamFile {
 public:
  ReadDreamFile(int nPart1, int nPart2);
  virtual ~ReadDreamFile();
  void SetAnalysisFile(const char* AnalysisFile, const char* prefix,
                       const char* Addon = "");
  void SetSigmaAnalysisFile(const char* PathAnalysisFile, const char* suffix);
  void ExtractResults(const TList *Results);
  void ReadkTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "");
  void ReadmTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "");
  void ReaddEtadPhiAtRadHists(const unsigned int nMaxMix,
                              const char* AnalysisFile, const char* prefix,
                              const char* Addon = "");
  void ReaddEtadPhiHists(const unsigned int NBinsmT, const char* AnalysisFile,
                         const char* prefix, const char* Addon = "");
  DreamDist* GetPairDistributions(int iPart1, int iPart2, const char* name);
  DreamKayTee* GetkTPairDistributions(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamKayTee* GetmTPairDistributions(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamdEtadPhi* GetdEtadPhiDistribution(int iPart1, int iPart2, int iAPart1,
                                         int iAPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiAtRadDistribution(int iPart1, int iPart2, int iMix1,
                                              int iAPart1, int iAPart2,
                                              int iMix2, int iRad,
                                              bool smallkStar);
  const int fNPart1;
  const int fNPart2;
 private:
  TH1F*** fSE;
  TH2F*** fSEMult;
  TH2F*** fSEkT;
  TH2F*** fSEmT;
  TH2F**** fSEdEtadPhimT;
  TH2F***** fSEdEtadPhiAtRad;
  TH2F***** fSEdEtadPhiAtRadSmallkStar;
  TH2F*** fSEdEtadPhi;
  TH1F*** fME;
  TH2F*** fMEMult;
  TH2F*** fMEkT;
  TH2F*** fMEmT;
  TH2F**** fMEdEtadPhimT;
  TH2F*** fMEdEtadPhi;
  TH2F***** fMEdEtadPhiAtRad;
  TH2F***** fMEdEtadPhiAtRadSmallkStar;
};

#endif /* DREAMFUNCTION_READDREAMFILE_H_ */
