/*
 * PlayWithCats.h
 *
 *  Created on: Jan 14, 2019
 *      Author: schmollweger
 */

#ifndef GENTLEKITTY_PLAYWITHCATS_H_
#define GENTLEKITTY_PLAYWITHCATS_H_
#include "TH1F.h"
#include "TFile.h"

class PlayWithCats {
 public:
  PlayWithCats();
  virtual ~PlayWithCats();
  void PlotPotentials();
  void PlotPotentialSum();
  void ExtractUncertaintyData(const char* inPath);
  void ExtractUncertaintyFit(const char* inFile);
  void GenerateDefault();
  void GenerateCoulombOnly();
  void CloseFile() {
    fOutFile->Close();
  }
  ;
//  void ShiftBinning() {
//    fStartBin++;
//    fMomBins++;
//    fkMin-=fDelta;
//  }
 private:
  TH1F* fCFHistData;
  TH1F* fCFHistErrorFit;
  TH1F* fCFHistDefaultData;
  TH1F* fCFHistDefaultFit;
  TH1F* fCFCoulombOnly;
  double fQCDTime;
  const double fMass_p;
  const double fMass_Xim;
  int fMomBins;
  float fkMin;
  float fkMax;
  const float fRadius;
  TFile *fOutFile;
};

#endif /* GENTLEKITTY_PLAYWITHCATS_H_ */
