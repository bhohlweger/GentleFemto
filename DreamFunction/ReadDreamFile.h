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

#include "DreamDist.h"
class ReadDreamFile {
 public:
  ReadDreamFile(int nPart1,int nPart2);
  virtual ~ReadDreamFile();
  void SetAnalysisFile(const char* AnalysisFile, const char* prefix);
  DreamDist* GetPairDistributions(int iPart1,int iPart2, const char* name);
  const int     fNPart1;
  const int     fNPart2;
  TH1F***  fSE;
  TH2F***  fSEMult;
  TH1F***  fME;
  TH2F***  fMEMult;

};

#endif /* DREAMFUNCTION_READDREAMFILE_H_ */
