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
class MakeHistosGreat {
 public:
  MakeHistosGreat();
  virtual ~MakeHistosGreat();
  void FormatHistogram(TH1* hist, unsigned int marker, unsigned int color);
  void DrawAndStore(TH1* hist, const char* outname);
  void DrawAndStore(TH2* hist, const char* outname);
  static void SetStyle(bool graypalette = false, bool title = true);
};

#endif /* FORGIVINGQA_MAKEHISTOSGREAT_H_ */
