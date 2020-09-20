#ifndef GENTLEKITTY_CATSINPUTPHI_H_
#define GENTLEKITTY_CATSINPUTPHI_H_
#include "CATSInput.h"
#include "TGraphAsymmErrors.h"

class CATSInputPhi : public CATSInput {
 public:

  CATSInputPhi();
  virtual ~CATSInputPhi();

  void ReadPhiCorrelationFile(const char* path, const char* trigger,
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
  unsigned int GetNPhi() const {
    return fnPhi;
  }
  float GetPhiPurity() const {
    return fPurityPhi;
  }
  float GetPhiPurityPt() const {
    return fPurityPhipt;
  }


 private:
  DreamCF* fCF_pPhi;
  DreamCF* fCF_SidebandUp;
  DreamCF* fCF_SidebandLow;
  unsigned int fnProtons;
  unsigned int fnAntiProtons;
  unsigned int fnPhi;
  float fPurityPhi;
  float fPurityPhipt;
};

#endif /* GENTLEKITTY_CATSINPUT_H_ */
