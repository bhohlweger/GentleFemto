#ifndef GENTLEKITTY_CATSINPUTSIGMA0_H_
#define GENTLEKITTY_CATSINPUTSIGMA0_H_
#include "CATSInput.h"

class CATSInputSigma0 : public CATSInput {
 public:

  CATSInputSigma0();
  virtual ~CATSInputSigma0();

  void ReadSigma0CorrelationFile(const char* path, const char* trigger,
                                 const char* suffixChar);
  void ObtainCFs(int rebin, float normleft, float normright, int rebinSyst = 1,
                 bool isAllCF = true);
  TH1F* GetCF(TString pair, TString hist);
 private:
  DreamCF* fCF_pSigma;
  DreamCF* fCF_SidebandUp;
  DreamCF* fCF_SidebandLow;
};

#endif /* GENTLEKITTY_CATSINPUT_H_ */
