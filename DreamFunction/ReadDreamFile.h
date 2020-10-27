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
  void SetAnalysisFile(const char* PathAnalysisFile, const char* Path,
                       const char* Prefix, const char* Addon);
  void SetAnalysisFileSample(const char* AnalysisFile, const char* prefix,
                       const char* Addon = "");
  void SetAnalysisFileAncestors(const char* AnalysisFile, const char* prefix,
                       const char* Addon = "");
  void ExtractResults(const TList *Results);
  void ExtractResultsAncestors(const TList *Results);
  void ReadkTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "");
  void ReadmTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "");
  void ReadmTHistosAncestors(const char* AnalysisFile, const char* prefix,
                    const char* addon = "");
  void ReadAndProjectmTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "", double kcut = 200.);
  void ReadAndProjectmTHistosAncestors(const char* AnalysisFile, const char* prefix,
                    const char* addon = "", double kcut = 200.);
  void ReadAndProjectkTHistos(const char* AnalysisFile, const char* prefix,
                    const char* addon = "", double kcut = 200.);
  void ReaddEtadPhiAtRadHists(const unsigned int nMaxMix,
                              const char* AnalysisFile, const char* prefix,
                              const char* Addon = "");
  void ReaddEtadPhiHists(const unsigned int NBinsmT, const char* AnalysisFile,
                         const char* prefix, const char* Addon = "");
  void ReaddEtadPhiHistsAncestors(const unsigned int NBinsmT, const char* AnalysisFile,
                         const char* prefix, const char* Addon = "");
  void ReadmTMultHistos(const char* AnalysisFile, const char* prefix,
			const char* addon, const int nmTBins);
  DreamDist* GetPairDistributions(int iPart1, int iPart2, const char* name);
  DreamDist* GetPairDistributionsCommon(int iPart1, int iPart2, const char* name);
  DreamDist* GetPairDistributionsNonCommon(int iPart1, int iPart2, const char* name);

  DreamKayTee* GetkTPairDistributions(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamKayTee* GetmTPairDistributions(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamKayTee* GetmTMultPairDistributions(int iPart1, int iPart2, int iAPart1,
					  int iAPart2, const int nmTBins);
  DreamKayTee* GetmTPairDistributionsCommon(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamKayTee* GetmTPairDistributionsCommon(int iPart1, int iPart2);
  DreamKayTee* GetmTPairDistributionsNonCommon(int iPart1, int iPart2, int iAPart1,
                                      int iAPart2);
  DreamKayTee* GetmTPairDistributionsNonCommon(int iPart1, int iPart2);
  DreamKayTee* GetmTPairDistributionsBBar(int iPart1, int iPart2);
  DreamdEtadPhi* GetdEtadPhiDistribution(int iPart1, int iPart2, int iAPart1,
                                         int iAPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiDistributionCommon(int iPart1, int iPart2, int iAPart1,
                                         int iAPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiDistributionNonCommon(int iPart1, int iPart2, int iAPart1,
                                         int iAPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiDistributionSingle(int iPart1, int iPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiDistributionSingleCommon(int iPart1, int iPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiDistributionSingleNonCommon(int iPart1, int iPart2, int imT = 0);
  DreamdEtadPhi* GetdEtadPhiAtRadDistribution(int iPart1, int iPart2, int iMix1,
                                              int iAPart1, int iAPart2,
                                              int iMix2, int iRad,
                                              bool smallkStar);
  void SetQuite() { fQuiet = true;}
  const int fNPart1;
  const int fNPart2;
 private:
  bool fQuiet;
  TH1F*** fSE;
  TH2F*** fSEMult;
  TH2F*** fSEkT;
  TH2F*** fSEmT;
  TH2F**** fSEmTMult;
  TH1F*** fSEmTProj;
  TH1F*** fProjmT;
  TH1F*** fMeanmT;
  TH2F**** fSEdEtadPhimT;
  TH2F***** fSEdEtadPhiAtRad;
  TH2F***** fSEdEtadPhiAtRadSmallkStar;
  TH2F*** fSEdEtadPhi;
  TH1F*** fME;
  TH2F*** fMEMult;
  TH2F*** fMEkT;
  TH2F*** fMEmT;
  TH2F**** fMEmTMult;
  TH2F**** fMEdEtadPhimT;
  TH2F*** fMEdEtadPhi;
  TH2F***** fMEdEtadPhiAtRad;
  TH2F***** fMEdEtadPhiAtRadSmallkStar;
  TH1F*** fSECommon;
  TH1F*** fSENonCommon;
  TH2F*** fSEMultCommon;
  TH2F*** fSEMultNonCommon;
  TH2F*** fSEmTCommon;
  TH2F*** fSEmTNonCommon;
  TH2F*** fSEdEtadPhiCommon;
  TH2F*** fSEdEtadPhiNonCommon;
};

#endif /* DREAMFUNCTION_READDREAMFILE_H_ */
