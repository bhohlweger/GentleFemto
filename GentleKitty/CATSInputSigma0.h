#ifndef GENTLEKITTY_CATSINPUTSIGMA0_H_
#define GENTLEKITTY_CATSINPUTSIGMA0_H_
#include "CATSInput.h"
#include "TGraphAsymmErrors.h"

class CATSInputSigma0 : public CATSInput {
 public:

  CATSInputSigma0();
  virtual ~CATSInputSigma0();

  void ReadSigma0CorrelationFile(const char* path, const char* trigger,
                                 const char* suffixChar);
  void ObtainCFs(int rebin, float normleft, float normright, int rebinSyst = 1,
                 bool isAllCF = true);
  void CountPairs(const char* path, const char* trigger,
                  const char* suffixChar);
  TH1F* GetCF(TString pair, TString hist);
  TGraphAsymmErrors* GetCFGr(TString pair, TString hist);
  unsigned int GetFemtoPairs(float kMin, float kMax, TString pair);
  unsigned int GetNProtons() const {
    return fnProtons;
  }
  unsigned int GetNAntiProtons() const {
    return fnAntiProtons;
  }
  unsigned int GetNProtonTotal() const {
    return fnProtons + fnAntiProtons;
  }
  unsigned int GetNSigma0() const {
    return fnSigma0;
  }
  float GetSigma0Purity() const {
    return fPuritySigma0;
  }
  float GetSigma0PurityPt() const {
    return fPuritySigma0pt;
  }


 private:
  DreamCF* fCF_pSigma;
  DreamCF* fCF_SidebandUp;
  DreamCF* fCF_SidebandLow;
  unsigned int fnProtons;
  unsigned int fnAntiProtons;
  unsigned int fnSigma0;
  float fPuritySigma0;
  float fPuritySigma0pt;
};

#endif /* GENTLEKITTY_CATSINPUT_H_ */
