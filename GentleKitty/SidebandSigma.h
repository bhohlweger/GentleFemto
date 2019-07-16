#ifndef GENTLEKITTY_SidebandSigmaSIGMA_H_
#define GENTLEKITTY_SidebandSigmaSIGMA_H_
#include "ReadDreamFile.h"
#include "SideBandFit.h"

class SidebandSigma : public SideBandFit {
 public:
  SidebandSigma();
  virtual ~SidebandSigma();
  void SetSideBandFile(const char* path, const char* trigger,
                       const char* suffix);
  void SideBandCFs();
  TGraphAsymmErrors* GetSideBandGraph(int iHist) {
    return (iHist < fSidebandCFGr.size() ? fSidebandCFGr.at(iHist) : nullptr);
  }

 private:
  ReadDreamFile* fAnalysisFile;
  std::vector<TGraphAsymmErrors*> fSidebandCFGr;
};
#endif
