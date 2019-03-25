#ifndef FORGIVINGQA_PERIODQA_H_
#define FORGIVINGQA_PERIODQA_H_
#include "TList.h"
#include "TH2F.h"
#include "ForgivingReader.h"
#include "ForgivingFitter.h"
#include "MakeHistosGreat.h"

class PeriodQA {
 public:
  PeriodQA();

  void SetDirectory(const char* dir);
  void ProcessQA(const char* prefix, const char* addon);
  void ProcessSigmaQA(const char* prefix, const char* addon);

  TH1F* PeriodQAHist(const char* name, const char* title);

 private:
  TString fDirectory;
  std::vector<TString> fPeriods;
};

#endif //FORGIVINGQA_PERIODQA_H_
