/*
 * CATSInput.h
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */

#ifndef GENTLEKITTY_CATSINPUT_H_
#define GENTLEKITTY_CATSINPUT_H_
#include "TString.h"
#include "TH2F.h"
#include "ReadDreamFile.h"
#include <vector>
class CATSInput {
 public:

  CATSInput();
  virtual ~CATSInput();
  void SetCalibBaseDir(const char* path) {
    fNameBasedir.Clear();
    fNameBasedir.Append(path);
  }
  ;
  void SetMomResFileName(const char* filename, int Fraction_Res = 2,
                         double UnitConv_Res = 1) {
    fFraction_Res = Fraction_Res;
    fUnitConv_Res = UnitConv_Res;
    fNameMomResFile.Clear();
    fNameMomResFile.Append(filename);
  }
  ;
  void ReadResFile();
  TH2F* GetResFile(int iPair) const {
    return fRes.size() > iPair ? fRes[iPair] : nullptr;
  }
  ;
  void SetSigmaFileName(const char* filename, int Fraction_Sig = 1,
                        double UnitConv_Sig = 1) {
    fFraction_Sig = Fraction_Sig;
    fUnitConv_Sig = UnitConv_Sig;
    fNameSigmaFile.Clear();
    fNameSigmaFile.Append(filename);
  }
  ;
  void ReadSigmaFile();
  TH2F* GetSigmaFile(int iPair) const {
    return fSigma.size() > iPair ? fSigma[iPair] : nullptr;
  }
  ;
  void ReadCorrelationFile(const char* path, const char* prefix = "MB",
                           const char* suffix = "");
  void ObtainCFs(int rebin, float normleft, float normright);
  DreamCF* ObtainCFSyst(int rebin, const char* name, DreamDist* ppDist,
                        DreamDist* ApApDist, DreamDist* ppFake = nullptr,
                        DreamDist* ApApFake = nullptr);
  TH1F* GetCF(TString pair, TString hist);
  void AddSystematics(TString SysFile, TH1F* hist);
  void SetNormalization(float normleft, float normright) {
    fnormalizationLeft = normleft;
    fnormalizationRight = normright;
  }
  void SetFixedkStarMinBin(bool doIt, float kMin = 0.) {
    fFixBinningExternal = doIt;
    fFixkMin = kMin;
  }
 protected:
  ReadDreamFile* fDreamFile;
  float fnormalizationLeft;
  float fnormalizationRight;
  DreamCF* fCF_pp;
 private:
  bool fFixBinningExternal;
  float fFixkMin;
  TString fNameBasedir;
  TString fNameMomResFile;
  TString fNameSigmaFile;
  int fFraction_Res;
  int fFraction_Sig;
  double fUnitConv_Res;
  double fUnitConv_Sig;
  std::vector<TH2F*> fRes;
  std::vector<TH2F*> fSigma;
  DreamCF* fCF_pL;
  DreamCF* fCF_LL;
  DreamCF* fCF_pXi;
};

#endif /* GENTLEKITTY_CATSINPUT_H_ */
