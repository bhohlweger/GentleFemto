/*
 * MakeHistosGreat.h
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_MAKEHISTOSGREAT_H_
#define FORGIVINGQA_MAKEHISTOSGREAT_H_
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "ForgivingFitter.h"
#include <vector>

struct DrawStyle {
  int drawMarker = 0;
  unsigned int drawColor = 1;
  int drawLineColor =  kOrange - 1;
  int drawSignalFitColor =  kBlue;
  int drawBackgroundFitColor =  kBlue;
  float drawSize = 1.;
};

class MakeHistosGreat {
 public:
  MakeHistosGreat();
  MakeHistosGreat(const char* outfile);
  virtual ~MakeHistosGreat();
  void FormatHistogram(TH1* hist, DrawStyle styler);
  void FormatHistogram(TH1* hist, unsigned int marker, unsigned int color,
                       float size = 1);
  void FormatSmallHistogram(TH1* hist, DrawStyle styler);
  void FormatSmallHistogram(TH1* hist, unsigned int marker, unsigned int color,
                            float size = 1);
  void FormatHistogram(TH2 *histo);
  void DrawAndStore(std::vector<TH1*> hist, const char* outname,
                    const char* drawOption = "");
  void DrawOnPad(std::vector<TH1*> hist, TPad* TPain, const char* drawOption =
                     "");
  void DrawLogYAndStore(std::vector<TH1*> hist, const char* outname,
                        const char* drawOption = "");
  void DrawAndStore(std::vector<TH2*> hist, const char* outname,
                    const char* drawOption = "");
  void DrawLogZAndStore(std::vector<TH2*> hist, const char* outname,
                        const char* drawOption = "");
  static void SetStyle(bool title = true);
  void SetTightMargin(bool set = false) {
    fTightMargin = set;
  }
  ;
  void DrawLatexLabel(float pTMin, float pTMax, ForgivingFitter* fit, TPad* pad,
                      const char* part, float xPos, float yPos, float offset = 0.07);
  void DrawPerformance(ForgivingFitter* fit, TPad* pad, const char* part,
                       float xPos, float yPos, float pTmin = -1, float pTmax =
                           -1);
  void DrawPublication(ForgivingFitter* fit, TPad* pad, const char* part,
                       float xPos, float yPos, float pTmin = -1, float pTmax =
                           -1);
  void DrawLine(TPad* pad, float xMin, float xMax, float yMin, float yMax, int color =  kOrange - 1);

  void DumpToFile(TCanvas *c, TH1* hist, const char* name);
 private:
  bool fTightMargin;
  TFile *fOutfile;
};

inline void MakeHistosGreat::FormatHistogram(TH1* hist, DrawStyle styler) {
  FormatHistogram(hist, styler.drawMarker, styler.drawColor, styler. drawSize);
}

inline void MakeHistosGreat::FormatSmallHistogram(TH1* hist, DrawStyle styler) {
  FormatSmallHistogram(hist, styler.drawMarker, styler.drawColor, styler. drawSize);
}

inline void MakeHistosGreat::DumpToFile(TCanvas *c, TH1* hist, const char* name) {
  if (fOutfile) {
    if (c) {
      c->Write(name);
    }
    if (hist) {
      hist->Write(name);
    }
  }
}

#endif /* FORGIVINGQA_MAKEHISTOSGREAT_H_ */
