/*
 * SideBandFit.h
 *
 *  Created on: Nov 6, 2018
 *      Author: hohlweger
 */

#ifndef GENTLEKITTY_SIDEBANDFIT_H_
#define GENTLEKITTY_SIDEBANDFIT_H_
#include "ReadDreamFile.h"

class SideBandFit {
 public:
  SideBandFit();
  virtual ~SideBandFit();
  void SetSideBandFile(const char* path, const char* suffixUp,
                       const char* suffixDown);
  void SideBandCFs(bool doQA);
  void SetNormalizationRange(float left, float right) {
    fnormleft = left;
    fnormright = right;
  }
  ;
  void SetRebin(int iRebin) {
    fRebin = iRebin;
  }
  ;
  void SetSidebandCF(TH1F* hCK) {
    fSideBandCFs.push_back(hCK);
  }
  TH1F* AddCF(TH1F* CF1, TH1F* CF2, TH1F* CF3, TH1F* CF4, const char* name);
  TH1F* AddCF(TH1F* CF1, TH1F* CF2, const char* name);
  TH1F* GetSideBands(int iHist) {
    return (iHist < fSideBandCFs.size() ? fSideBandCFs.at(iHist) : nullptr);
  }
  ;
  void WriteOutput(const char* outputPath);
  static double Parameterization(const double& Momentum,
                                 const double* SourcePar, const double* PotPar);
  static double ParameterizationROOT(double* Momentum, double* PotPar);
  void FitSideBands(TH1F* cfSide, double* potPar);
 private:
  ReadDreamFile* fAnalysisFileUp;
  ReadDreamFile* fAnalysisFileDown;
  std::vector<DreamPair*> fSideBands;
  std::vector<TH1F*> fSideBandCFs;
  float fnormleft;  // in MeV
  float fnormright;  // in MeV
  int fRebin;
};
#endif /* GENTLEKITTY_SIDEBANDFIT_H_ */
