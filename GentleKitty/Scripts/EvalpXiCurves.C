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
void EvalpXiCurves(const char* cfpath, const char* varFolder) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
//  CATSinput->SetFixedkStarMinBin(true, 0.008);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  TString InputFile = TString::Format("%s/AnalysisResults.root", cfpath);
  DreamFile->SetAnalysisFile(InputFile.Data(), "MB");
  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXi = CATSinput->ObtainCFSyst(5, "pXi", pXi, ApAXi);
  TString HistpXiName = "hCk_ReweightedpXiMeV_1";
//  CATSinput->ReadCorrelationFile(cfpath);
//  CATSinput->ObtainCFs(5, 240, 340);
//  TString HistName = Form("hCk_Reweighted%sMeV_0", pair);
//  TH1F* CF_Histo = CATSinput->GetCF(pair, HistName);
  TH1F* CF_Histo = CFpXi->FindCorrelationFunction(HistpXiName.Data());
  CF_Histo->GetYaxis()->SetRangeUser(0.9, 2.5);
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
  MygraphName += "pXimGraph";
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
  TNtuple *outFit = new TNtuple("pXiFit", "pXiFit",
                                "delta1:delta2:delta3:delta4:delta5:delta6:"
                                "delta7:delta8:delta9:delta10:delta11:"
                                "delta12:delta13:delta14:delta15:delta16:"
                                "delta17:delta18");
  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    str = entry;
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
          for (int iPnt = 0; iPnt < tupleLength; ++iPnt) {
            tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
                - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
            if (iPnt == tupleLength - 1) {
              Graph->SetLineColor(3);
              outFit->Fill(tupler);
            }
          }
          Graph->Draw("L3same");
        }
      }
      file->Close();
    }
  }
  c1->SaveAs(Form("%s/CF_%s_model.pdf", varFolder, "pXi"));
  CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
  c1->SaveAs(Form("%s/CF_%s_model_zoom.pdf", varFolder, "pXi"));
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "RECREATE");
  output->cd();
  c1->Write();
  outFit->Write();
  output->Close();
  delete CATSinput;
}

void EvalError(const char* cfpath, const char* varFolder) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
//  CATSinput->SetFixedkStarMinBin(true, 0.008);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  TString InputFile = TString::Format("%s/AnalysisResults.root", cfpath);
  DreamFile->SetAnalysisFile(InputFile.Data(), "MB");
  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXi = CATSinput->ObtainCFSyst(5, "pXi", pXi, ApAXi);
  TString HistpXiName = "hCk_ReweightedpXiMeV_1";
  TH1F* CF_Histo = CFpXi->FindCorrelationFunction(HistpXiName.Data());
  CF_Histo->GetXaxis()->SetRangeUser(0, 500);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_Histo->SetStats(false);
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(CF_Histo);
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "update");
  TNtuple* outFit = (TNtuple*) output->Get("pXiFit");
  TGraph grUp;
  TGraph grDefault;
  TGraph grLow;

  for (int iBranches = 1; iBranches < tupleLength + 1; ++iBranches) {
    float kVal = CF_Histo->GetBinCenter(iBranches);
    float CkExp = CF_Histo->GetBinContent(iBranches);

    outFit->Draw(Form("delta%u >> h", iBranches));
    TH1F* hist = (TH1F*) gROOT->FindObject("h");

    double binLimLow = hist->GetXaxis()->GetBinLowEdge(
        hist->FindFirstBinAbove(1, 1));
    double binLimUp = hist->GetXaxis()->GetBinUpEdge(
        hist->FindLastBinAbove(1, 1));

    double Delta = (binLimUp - binLimLow) / (2 * TMath::Sqrt(12));

    double CoulombdefVal = (binLimUp + binLimLow) / 2.;

    grDefault.SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal);
    grUp.SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal + Delta);
    grLow.SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal - Delta);

    delete hist;

  }
  float chisqDefault = 0;
  float chisqUp = 0;
  float chisqDown = 0;
  int ndf = -1;
  for (int iBin = 1; iBin < CF_Histo->GetXaxis()->FindBin(200); ++iBin) {
    ndf++;
    double kVal = CF_Histo->GetBinCenter(iBin);
    double CkVal = CF_Histo->GetBinContent(iBin);
    double CkErrStat = CF_Histo->GetBinError(iBin);

    double chiDefault = (CkVal - grDefault.Eval(kVal)) / CkErrStat;
    chisqDefault += chiDefault * chiDefault;

    double chiUp = (CkVal - grUp.Eval(kVal)) / CkErrStat;
    chisqUp += chiUp * chiUp;

    double chiDown = (CkVal - grLow.Eval(kVal)) / CkErrStat;
    chisqDown += chiDown * chiDown;
  }
  std::cout << "ndf " << ndf << std::endl;
  std::cout << "Chisq Default " << chisqDefault << std::endl;
  std::cout << "Chisq Up " << chisqUp << std::endl;
  std::cout << "Chisq Down " << chisqDown << std::endl;
  std::cout << std::endl;


  double pvalXiDefault = TMath::Prob(chisqDefault, round(ndf));
  double nSigmaXiDefault = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiDefault);

  double pvalXiUp = TMath::Prob(chisqUp, round(ndf));
  double nSigmaXiUp = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiUp);

  double pvalXiDown = TMath::Prob(chisqDown, round(ndf));
  double nSigmaXiDown = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiDown);

  std::cout << "=============================================\n";
  std::cout << "==============SIGMA VALUES BABY==============\n";
  std::cout << "=============================================\n";
  std::cout << "///////// nSigma Default = " << nSigmaXiDefault << "/////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "//////////// nSigma Up= " << nSigmaXiUp << "////////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "/////////// nSigma Down= " << nSigmaXiDown << "///////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "=============================================\n";
  output->cd();
  grUp.SetName("pXimGraphUpperLim");
  grUp.SetLineColor(3);
  grDefault.SetName("pXimGraphDefault");
  grDefault.SetLineColor(3);
  grLow.SetName("pXimGraphLowerLim");
  grLow.SetLineColor(3);

  auto cErr = new TCanvas("cErr","cErr");
  CF_Histo->Draw();
  grUp.Draw("SAME");
  grDefault.Draw("SAME");
  grLow.Draw("SAME");



  grUp.Write();
  grDefault.Write();
  grLow.Write();

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
  TFile* FileUp = TFile::Open(
      Form("%s/GraphFile_pXi_iter1_Var155.root", varFolder), "update");
  TFile* FileDefault = TFile::Open(
      Form("%s/GraphFile_pXi_iter1_Var13.root", varFolder), "update");
  TFile* FileDown = TFile::Open(
      Form("%s/GraphFile_pXi_iter1_Var17.root", varFolder), "update");
  grSidebandUp = (TGraph*) FileUp->Get("SideBandStrongWithLambda");
  grSidebandDefault = (TGraph*) FileDefault->Get("SideBandStrongWithLambda");
  grSidebandDown = (TGraph*) FileDown->Get("SideBandStrongWithLambda");
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
  std::cout << "EvalpXiCurves \n";
  EvalpXiCurves(argv[1], argv[2]);
  std::cout << "EvalError \n";
  EvalError(argv[1], argv[2]);
  std::cout << "SidebandCurves \n";
  SidebandCurves(argv[2]);
  std::cout << "CombineIntoOneFile \n";
  CombineIntoOneFile(argv[3], argv[2]);
  return 0;
}
