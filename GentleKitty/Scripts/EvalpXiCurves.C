#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TNtuple.h"
#include "CATSInput.h"
#include "DreamPlot.h"
#include <iostream>
#include "TMath.h"
static const int tupleLength = 18;
//ONLY WORKS FOR 20 MEV BINNING
void EvalpXiCurves(const char* pair, const char* cfpath,
                   const char* varFolder) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->ReadCorrelationFile(cfpath);
  CATSinput->ObtainCFs(5, 240, 340);
  TString HistName = Form("hCk_Reweighted%sMeV_0", pair);
  TH1F* CF_Histo = CATSinput->GetCF(pair, HistName);
  CF_Histo->GetXaxis()->SetRangeUser(0, 500);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_Histo->SetStats(false);
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(CF_Histo);
  const char* beg = "GraphFile";
  const char* ext = ".root";

  char* dir = gSystem->ExpandPathName(varFolder);
  void* dirp = gSystem->OpenDirectory(dir);
  TString MygraphName = "";
  if (TString(pair) == "pp") {
    MygraphName += "Fit";
  } else if (TString(pair) == "pXi") {
    MygraphName += "pXimGraph";
  }
  const char* entry;
  std::vector<TGraph*> Graph;
  Int_t n = 0;
  TString str;
  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->cd(0);
  c1->SetCanvasSize(1920, 1280);
  c1->SetMargin(0.15, 0.05, 0.2, 0.05);
  c1->cd();
  CF_Histo->Draw();
  TNtuple *outCoulomb = new TNtuple("coulomb", "coulomb",
                                    "delta1:delta2:delta3:delta4:delta5:delta6:"
                                    "delta7:delta8:delta9:delta10:delta11:"
                                    "delta12:delta13:delta14:delta15:delta16:"
                                    "delta17:delta18");
  TNtuple *outCoulombStrong = new TNtuple(
      "coulombStrong", "coulombStrong",
      "delta1:delta2:delta3:delta4:delta5:delta6:"
      "delta7:delta8:delta9:delta10:delta11:"
      "delta12:delta13:delta14:delta15:delta16:"
      "delta17:delta18");
  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    str = entry;
//    std::cout << str << std::endl;
    if (str.EndsWith(ext) && str.BeginsWith(beg)) {
      TFile* file = TFile::Open(gSystem->ConcatFileName(dir, entry));
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*) next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TGraph"))
          continue;
        TGraph *Graph = (TGraph*) key->ReadObj();
        TString graphName = Form("%s", Graph->GetName());
        float tupler[tupleLength];
        if (graphName.Contains(MygraphName.Data())) {
          TString nameCopy = str.Data();
          nameCopy.Replace(0, nameCopy.Index("Var") + 3, "");
          nameCopy.Replace(nameCopy.Index(".root"), 5, "");
          for (int iPnt = 0; iPnt < tupleLength; ++iPnt) {
            if (graphName.Contains("COULOMB")) {
              tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
                  - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
              if (iPnt == tupleLength - 1) {
                Graph->SetLineColor(3);
                outCoulomb->Fill(tupler);
              }
            } else {
              tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
                  - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
              if (iPnt == tupleLength - 1) {
                Graph->SetLineColor(2);
                outCoulombStrong->Fill(tupler);
              }
            }
          }
          Graph->Draw("L3same");
        }
      }
      file->Close();
    }
  }
  c1->SaveAs(Form("%s/CF_%s_model.pdf", varFolder, pair));
  CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
  c1->SaveAs(Form("%s/CF_%s_model_zoom.pdf", varFolder, pair));
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "RECREATE");
  output->cd();
  c1->Write();
  outCoulomb->Write();
  outCoulombStrong->Write();
  output->Close();
  delete CATSinput;
}

void EvalError(const char* cfpath, const char* varFolder) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->ReadCorrelationFile(cfpath);
  CATSinput->ObtainCFs(5, 240, 340);
  TString HistName = Form("hCk_Reweighted%sMeV_0", "pXi");
  TH1F* CF_Histo = CATSinput->GetCF("pXi", HistName);
  CF_Histo->GetXaxis()->SetRangeUser(0, 500);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_Histo->SetStats(false);
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(CF_Histo);
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "update");
  TNtuple* outCoulomb = (TNtuple*) output->Get("coulomb");
  TNtuple* outCoulombStrong = (TNtuple*) output->Get("coulombStrong");
  TGraph grCoulombUp;
  TGraph grCoulombDefault;
  TGraph grCoulombLow;
  TGraph grCoulombStrongUp;
  TGraph grCoulombStrongDefault;
  TGraph grCoulombStrongLow;

  for (int iBranches = 1; iBranches < tupleLength + 1; ++iBranches) {
    float kVal = CF_Histo->GetBinCenter(iBranches);
    float CkExp = CF_Histo->GetBinContent(iBranches);

    outCoulomb->Draw(Form("delta%u >> h", iBranches));
    TH1F* histCoulomb = (TH1F*) gROOT->FindObject("h");

    double CoulombbinLimLow = histCoulomb->GetXaxis()->GetBinLowEdge(
        histCoulomb->FindFirstBinAbove(1, 1));
    double CoulombbinLimUp = histCoulomb->GetXaxis()->GetBinUpEdge(
        histCoulomb->FindLastBinAbove(1, 1));

    double DeltaCoulomb = (CoulombbinLimUp - CoulombbinLimLow)
        / (2 * TMath::Sqrt(12));

    double CoulombdefVal = (CoulombbinLimUp + CoulombbinLimLow) / 2.;

    grCoulombDefault.SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal);
    grCoulombUp.SetPoint(iBranches - 1, kVal,
                         CkExp - CoulombdefVal + DeltaCoulomb);
    grCoulombLow.SetPoint(iBranches - 1, kVal,
                          CkExp - CoulombdefVal - DeltaCoulomb);

    delete histCoulomb;

//
    outCoulombStrong->Draw(Form("delta%u >> h", iBranches));
    TH1F* histCoulombStrong = (TH1F*) gROOT->FindObject("h");

    double CoulombStrongbinLimLow =
        histCoulombStrong->GetXaxis()->GetBinLowEdge(
            histCoulombStrong->FindFirstBinAbove(1, 1));
    double CoulombStrongbinLimUp = histCoulombStrong->GetXaxis()->GetBinUpEdge(
        histCoulombStrong->FindLastBinAbove(1, 1));

    double DeltaCoulombStrong = (CoulombStrongbinLimUp - CoulombStrongbinLimLow)
        / (2 * TMath::Sqrt(12));
    double CoulombStrongdefVal =
        (CoulombStrongbinLimUp + CoulombStrongbinLimLow) / 2.;

    grCoulombStrongDefault.SetPoint(iBranches - 1, kVal,
                                    CkExp - CoulombStrongdefVal);
    grCoulombStrongUp.SetPoint(
        iBranches - 1, kVal, CkExp - CoulombStrongdefVal + DeltaCoulombStrong);
    grCoulombStrongLow.SetPoint(
        iBranches - 1, kVal, CkExp - CoulombStrongdefVal - DeltaCoulombStrong);
    delete histCoulombStrong;

  }
  float chisqCoulombDefault = 0;
  float chisqCoulombUp = 0;
  float chisqCoulombDown = 0;
  float chisqCoulombStrongDefault = 0;
  float chisqCoulombStrongUp = 0;
  float chisqCoulombStrongDown = 0;
  int ndf = -1;
  for (int iBin = 1; iBin < CF_Histo->GetXaxis()->FindBin(200); ++iBin) {
    ndf++;
    double kVal = CF_Histo->GetBinCenter(iBin);
    double CkVal = CF_Histo->GetBinContent(iBin);
    double CkErrStat = CF_Histo->GetBinError(iBin);

    double chiCoulombDefault = (CkVal - grCoulombDefault.Eval(kVal)) / CkErrStat;
    chisqCoulombDefault += chiCoulombDefault * chiCoulombDefault;

    double chiCoulombUp = (CkVal - grCoulombUp.Eval(kVal)) / CkErrStat;
    chisqCoulombUp += chiCoulombUp * chiCoulombUp;

    double chiCoulombDown = (CkVal - grCoulombLow.Eval(kVal)) / CkErrStat;
    chisqCoulombDown += chiCoulombDown * chiCoulombDown;

    double chiCoulombStrongDefault = (CkVal - grCoulombStrongDefault.Eval(kVal))
        / CkErrStat;
    chisqCoulombStrongDefault += chiCoulombStrongDefault * chiCoulombStrongDefault;

    double chiCoulombStrongUp = (CkVal - grCoulombStrongUp.Eval(kVal))
        / CkErrStat;
    chisqCoulombStrongUp += chiCoulombStrongUp * chiCoulombStrongUp;

    double chiCoulombStrongDown = (CkVal - grCoulombStrongLow.Eval(kVal))
        / CkErrStat;
    chisqCoulombStrongDown += chiCoulombStrongDown * chiCoulombStrongDown;
  }
  std::cout << "ndf " << ndf << std::endl;
  std::cout << "Chisq Coulomb Default " << chisqCoulombDefault << std::endl;
  std::cout << "Chisq Coulomb Up " << chisqCoulombUp << std::endl;
  std::cout << "Chisq Coulomb Down " << chisqCoulombDown << std::endl;
  std::cout << std::endl;
  std::cout << "Chisq Coulomb + Strong Default " << chisqCoulombStrongDefault << std::endl;
  std::cout << "Chisq Coulomb + Strong Up " << chisqCoulombStrongUp << std::endl;
  std::cout << "Chisq Coulomb + Strong Down " << chisqCoulombStrongDown << std::endl;

  double pvalXiCoulombDefault = TMath::Prob(chisqCoulombDefault, round(ndf));
  double nSigmaXiCoulombDefault = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiCoulombDefault);

  double pvalXiCoulombUp = TMath::Prob(chisqCoulombUp, round(ndf));
  double nSigmaXiCoulombUp = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiCoulombUp);

  double pvalXiCoulombDown = TMath::Prob(chisqCoulombDown, round(ndf));
  double nSigmaXiCoulombDown = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiCoulombDown);

  double pvalXiCoulombStrongDefault = TMath::Prob(chisqCoulombStrongDefault, round(ndf));
  double nSigmaXiCoulombStrongDefault = TMath::Sqrt(2)
      * TMath::ErfcInverse(pvalXiCoulombStrongDefault);
  double pvalXiCoulombStrongUp = TMath::Prob(chisqCoulombStrongUp, round(ndf));
  double nSigmaXiCoulombStrongUp = TMath::Sqrt(2)
      * TMath::ErfcInverse(pvalXiCoulombStrongUp);
  double pvalXiCoulombStrongDown = TMath::Prob(chisqCoulombStrongDown, round(ndf));
  double nSigmaXiCoulombStrongDown = TMath::Sqrt(2)
      * TMath::ErfcInverse(pvalXiCoulombStrongDown);
  std::cout << "=============================================\n";
  std::cout << "==============SIGMA VALUES BABY==============\n";
  std::cout << "=============================================\n";
  std::cout << "////// nSigma Coulomb Default = " << nSigmaXiCoulombDefault << "//////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "///////// nSigma Coulomb Up= " << nSigmaXiCoulombUp << "/////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "//////// nSigma Coulomb Down= " << nSigmaXiCoulombDown << "////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "//nSigma Coulomb + Strong Default = " << nSigmaXiCoulombStrongDefault << "//"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "/////nSigma Coulomb + Strong Up= " << nSigmaXiCoulombStrongUp << "/////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "////nSigma Coulomb + Strong Down = " << nSigmaXiCoulombStrongDown << "///"
            << std::endl;
  std::cout << "=============================================\n";
  output->cd();
  grCoulombUp.SetName("pXimGraphUpperLim_COULOMB");
  grCoulombUp.SetLineColor(3);
  grCoulombUp.Write();
  grCoulombDefault.SetName("pXimGraphDefault_COULOMB");
  grCoulombDefault.SetLineColor(3);
  grCoulombDefault.Write();
  grCoulombLow.SetName("pXimGraphLowerLim_COULOMB");
  grCoulombLow.SetLineColor(3);
  grCoulombLow.Write();
  grCoulombStrongUp.SetName("pXimGraphUpperLim");
  grCoulombStrongUp.SetLineColor(2);
  grCoulombStrongUp.Write();
  grCoulombStrongDefault.SetName("pXimGraphDefault");
  grCoulombStrongDefault.SetLineColor(2);
  grCoulombStrongDefault.Write();
  grCoulombStrongLow.SetName("pXimGraphLowerLim");
  grCoulombStrongLow.SetLineColor(2);
  grCoulombStrongLow.Write();
  CF_Histo->Write();



  output->Close();

  // N sigma stuff

  delete CATSinput;
}

void SidebandCurves(const char* varFolder) {
  TGraph* grSidebandUp;
  TGraph* grSidebandDefault;
  TGraph* grSidebandDown;
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "update");
  TFile* FileUp = TFile::Open(Form("%s/GraphFile_pXi_iter1_Var155.root", varFolder), "update");
  TFile* FileDefault = TFile::Open(Form("%s/GraphFile_pXi_iter1_Var13.root", varFolder), "update");
  TFile* FileDown = TFile::Open(Form("%s/GraphFile_pXi_iter1_Var17.root", varFolder), "update");
  grSidebandUp = (TGraph*)FileUp->Get("SideBandStrongWithLambda");
  grSidebandDefault = (TGraph*)FileDefault->Get("SideBandStrongWithLambda");
  grSidebandDown = (TGraph*)FileDown->Get("SideBandStrongWithLambda");
  output->cd();
  grSidebandUp->SetName("pXimGraphSidebandUp");
  grSidebandUp->Write();
  grSidebandDefault->SetName("pXimGraphSidebandDefault");
  grSidebandDefault->Write();
  grSidebandDown->SetName("pXimGraphSidebandDown");
  grSidebandDown->Write();
  output->Close();
  FileUp->Close();
  FileDefault->Close();
  FileDown->Close();
}

void CombineIntoOneFile(const char* PathToppFolder,
                        const char* PathTopXiFolder) {
  TFile* pXi = TFile::Open(Form("%s/outfile.root", PathTopXiFolder), "READ");
  TFile* pp = TFile::Open(
      Form("%s/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root",
           PathToppFolder),
      "update");
  TIter nextpXi(pXi->GetListOfKeys());
  TKey *keypXi;
  while ((keypXi = (TKey*) nextpXi())) {
    TClass *clpXi = gROOT->GetClass(keypXi->GetClassName());
    if (!clpXi->InheritsFrom("TGraph"))
      continue;
    TGraph *Graph = (TGraph*) keypXi->ReadObj();
    std::cout << Graph->GetName() << std::endl;
    pp->cd();
    Graph->Write(Graph->GetName());
  }
  pXi->Close();
  pp->Close();

}

int main(int argc, char *argv[]) {
  EvalpXiCurves(argv[1], argv[2], argv[3]);
  EvalError(argv[2], argv[3]);
  SidebandCurves(argv[3]);
  CombineIntoOneFile(argv[4], argv[3]);
  return 0;
}
