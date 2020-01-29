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
  void SetStyler(DrawStyle styler) { fStyler = styler; }
  void ProcessQA(const char* prefix, const char* addon);
  void ProcessSigmaQA(const char* prefix, const char* addon);

  TH1F* PeriodQAHist(const char* name, const char* title);

 private:
  double GetErrorNPart(double nEvt, double nPart, double nPartErr) const;
  TString fDirectory;
  std::vector<TString> fPeriods;
  DrawStyle fStyler;
};

inline
double PeriodQA::GetErrorNPart(double nEvt, double nPart,
                               double nPartErr) const {
  double nEvtSq = nEvt * nEvt;
  return std::sqrt(
      nPartErr * nPartErr / nEvtSq + nPart * nPart / (nEvtSq * nEvt));
}

#endif //FORGIVINGQA_PERIODQA_H_
