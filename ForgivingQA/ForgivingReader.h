/*
 * ForgivingReader.h
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_FORGIVINGREADER_H_
#define FORGIVINGQA_FORGIVINGREADER_H_
#include "TFile.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TH2.h"
class ForgivingReader {
 public:
  ForgivingReader();
  ForgivingReader(const char *filename, const char* prefix, const char* suffix =
                      "");
  virtual ~ForgivingReader();
  TList* GetListInDir(const char* pathToList);
  TList* GetListInList(TList* inList, std::vector<const char*> pathToList);
  TH1* Get1DHistInList(TList* inList, const char* histname);
  TH2* Get2DHistInList(TList* inList, const char* histname);
  TList* GetQA() {
    TString QAStd = Form("%sQA%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",QAStd.Data(),QAStd.Data()));
  }
  ;
  TList* GetEventCuts() {
    TString QAStd = Form("%sEvtCuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",QAStd.Data(),QAStd.Data()));
  }
  ;
  TList* GetTrackCuts() {
    TString TCStd = Form("%sTrackCuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",TCStd.Data(),TCStd.Data()));
  }
  ;
  TList* GetAntiTrackCuts() {
    TString ATCStd = Form("%sAntiTrackCuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",ATCStd.Data(),ATCStd.Data()));
  }
  ;
  TList* Getv0Cuts() {
    TString v0Std = Form("%sv0Cuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",v0Std.Data(),v0Std.Data()));
  }
  ;
  TList* GetAntiv0Cuts() {
    TString Av0Std = Form("%sAntiv0Cuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",Av0Std.Data(),Av0Std.Data()));
  }
  ;
  TList* GetCascadeCuts() {
    TString cascStd = Form("%sCascadeCuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",cascStd.Data(),cascStd.Data()));
  }
  ;
  TList* GetAntiCascadeCuts() {
    TString AcascStd = Form("%sAntiCascadeCuts%s", fPrefix, fSuffix);
    return GetListInDir(Form("%s/%s",AcascStd.Data(),AcascStd.Data()));
  }
  ;
 private:
  TFile* fInput;
  const char* fPrefix;
  const char* fSuffix;
};

#endif /* FORGIVINGQA_FORGIVINGREADER_H_ */
