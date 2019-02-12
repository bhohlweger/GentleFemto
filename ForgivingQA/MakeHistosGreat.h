/*
 * MakeHistosGreat.h
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_MAKEHISTOSGREAT_H_
#define FORGIVINGQA_MAKEHISTOSGREAT_H_
#include "TH1.h"
#include "TH2.h"
#include <vector>
class MakeHistosGreat {
 public:
  MakeHistosGreat();
  virtual ~MakeHistosGreat();
  void FormatHistogram(TH1* hist, unsigned int marker, unsigned int color);
  void FormatHistogram(TH2 *histo);
  void DrawAndStore(std::vector<TH1*> hist, const char* outname,
                    const char* drawOption = "");
  void DrawLogYAndStore(std::vector<TH1*> hist, const char* outname,
                        const char* drawOption = "");
  void DrawAndStore(std::vector<TH2*> hist, const char* outname,
                    const char* drawOption = "");
  static void SetStyle(bool title = true);
  void SetTightMargin(bool set = false) {
    fTightMargin = set;
  }
  ;
 private:
  bool fTightMargin;
};

#endif /* FORGIVINGQA_MAKEHISTOSGREAT_H_ */
