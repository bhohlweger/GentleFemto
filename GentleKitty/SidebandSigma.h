#ifndef GENTLEKITTY_SidebandSigmaSIGMA_H_
#define GENTLEKITTY_SidebandSigmaSIGMA_H_
#include "ReadDreamFile.h"
#include "SideBandFit.h"

class SidebandSigma : public SideBandFit {
 public:
  SidebandSigma();
  virtual ~SidebandSigma();
  void SetSideBandFile(const char* path, const char* suffix);
  void SideBandCFs();

 private:
  ReadDreamFile* fAnalysisFile;
};
#endif
