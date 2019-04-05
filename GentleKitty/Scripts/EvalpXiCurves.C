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
void EvalpXiCurves(const char* cfpath, const char* prefix,
                   const char* varFolder,float kStarMin) {
  CATSInput *CATSinput = new CATSInput("PXi");
  CATSinput->SetNormalization(0.240, 0.340);
  CATSinput->SetFixedkStarMinBin(true, kStarMin);
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  TString InputFile = TString::Format("%s/AnalysisResults.root", cfpath);
  DreamFile->SetAnalysisFile(InputFile.Data(), prefix);
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
  CF_Histo->GetXaxis()->SetRangeUser(kStarMin*1000, 500);
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
  CF_Histo->Write();
  output->Close();
  delete CATSinput;
}

void EvalError(const char* cfpath, const char* prefix, const char* varFolder, const char* sysPath, float kStarMin) {
  CATSInput *CATSinput = new CATSInput("PXi");
  CATSinput->SetNormalization(0.240, 0.340);
  CATSinput->SetFixedkStarMinBin(true, kStarMin);
  CATSinput->SetCalibBaseDir(sysPath);
  std::cout << "sysPath: " << sysPath << std::endl;
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  TString InputFile = TString::Format("%s/AnalysisResults.root", cfpath);
  DreamFile->SetAnalysisFile(InputFile.Data(), prefix);
  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXi = CATSinput->ObtainCFSyst(5, "pXi", pXi, ApAXi);
  TString HistpXiName = "hCk_ReweightedpXiMeV_1";
  TH1F* CF_Histo = CFpXi->FindCorrelationFunction(HistpXiName.Data());
  CF_Histo->GetXaxis()->SetRangeUser(0, 200);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  CF_Histo->GetYaxis()->CenterTitle(true);
  CF_Histo->SetStats(false);
  TH1F* CF_Draw = (TH1F*) CF_Histo->Clone(
      TString::Format("%sDraw", CF_Histo->GetName()));
  TH1F* Dummy = (TH1F*) CF_Histo->Clone("Dummy");
  Dummy->Reset();
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(CF_Histo);
  std::cout << "Error Before: " << CF_Histo->GetBinError(1) << std::endl;
  CATSinput->AddSystematics("C2totalsysPXi.root",CF_Histo);
  std::cout << "Error After: " << CF_Histo->GetBinError(1) << std::endl;
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "update");
  TNtuple* outFit = (TNtuple*) output->Get("pXiFit");
  TGraph grUp;
  TGraph grDefault;
  TGraph grLow;
  TGraph grChiPerPointUp;
  TGraph grChiPerPointDefault;
  TGraph grChiPerPointLow;

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
    double CkErrSyst = 0;
    double CkErr = TMath::Sqrt(CkErrStat*CkErrStat+CkErrSyst*CkErrSyst);
    double chiDefault = (CkVal - grDefault.Eval(kVal)) / CkErr;
    chisqDefault += chiDefault * chiDefault;

    double chiUp = (CkVal - grUp.Eval(kVal)) / CkErr;
    chisqUp += chiUp * chiUp;

    double chiDown = (CkVal - grLow.Eval(kVal)) / CkErr;
    chisqDown += chiDown * chiDown;

    grChiPerPointDefault.SetPoint(iBin - 1, kVal, chiDefault);
    grChiPerPointUp.SetPoint(iBin - 1, kVal, chiUp);
    grChiPerPointLow.SetPoint(iBin - 1, kVal, chiDown);
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
  TGraph SigmaGraph = TGraph();
  SigmaGraph.SetPoint(0,0,nSigmaXiDown);
  SigmaGraph.SetPoint(1,1,nSigmaXiDefault);
  SigmaGraph.SetPoint(2,2,nSigmaXiUp);
  SigmaGraph.SetName("pXiGraphSigma");
  output->cd();
  grUp.SetName("pXimGraphUpperLim");
  grUp.SetLineColor(3);
  grDefault.SetName("pXimGraphDefault");
  grDefault.SetLineColor(3);
  grLow.SetName("pXimGraphLowerLim");
  grLow.SetLineColor(3);
  auto cErr = new TCanvas("cErr", "cErr");
  //xlow ylow xup yup
  cErr->cd();
  TPad *p1 = new TPad("p1", "p1", 0.05, 0.3, 0.95, 0.99);
  p1->SetTopMargin(0.01);
  p1->SetBottomMargin(0);
  p1->SetFillStyle(4000);
  p1->SetFrameFillStyle(4000);
  p1->Draw();
  TPad *p2 = new TPad("p2", "p2", 0.05, 0.01, 0.95, 0.3);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.45);
  p2->Draw();

  p1->cd();
  CF_Draw->GetXaxis()->SetLabelSize(0);
  CF_Draw->Draw();
  grUp.Draw("SAME");
  grDefault.Draw("SAME");
  grLow.Draw("SAME");

  p2->cd();
  Dummy->GetXaxis()->SetTitleSize(0.1);
  Dummy->GetXaxis()->SetTitleOffset(0.9);
  Dummy->GetXaxis()->SetLabelSize(0.1);

  Dummy->GetYaxis()->SetTitle("#chi_{loc} (k*)");
  Dummy->GetYaxis()->SetRangeUser(-5., 5.);
  Dummy->GetYaxis()->SetNdivisions(2);
  Dummy->GetYaxis()->SetTitleSize(0.1);
  Dummy->GetYaxis()->SetLabelSize(0.1);
  Dummy->GetYaxis()->CenterTitle(true);
  Dummy->GetYaxis()->SetTitleOffset(0.4);
  Dummy->GetYaxis()->SetLabelOffset(0.02);

  Dummy->Draw("Y+");
  grChiPerPointDefault.Draw("L3SAME");
  grChiPerPointUp.Draw("L3SAME");
  grChiPerPointLow.Draw("L3SAME");

  cErr->SaveAs(TString::Format("%s/ChiSqQA.pdf", varFolder));
  grUp.Write();
  grDefault.Write();
  grLow.Write();
  SigmaGraph.Write();

  grChiPerPointDefault.SetName("PerPointChiSqDefault");
  grChiPerPointUp.SetName("PerPointChiSqUp");
  grChiPerPointLow.SetName("PerPointChiSqLow");

  grChiPerPointDefault.Write();
  grChiPerPointUp.Write();
  grChiPerPointLow.Write();

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
  grSidebandUp->SetName("pXiSidebandUp");
  grSidebandUp->Write();
  grSidebandDefault->SetName("pXiSidebandDefault");
  grSidebandDefault->Write();
  grSidebandDown->SetName("pXiSidebandDown");
  grSidebandDown->Write();
  output->Close();
  FileUp->Close();
  FileDefault->Close();
  FileDown->Close();
}

void CombineIntoOneFile(const char* PathTopXiFolder, const char* GraphName,
                        const char* PathToCombinedFolder, bool sideBand) {
  TFile* pXi = TFile::Open(Form("%s/outfile.root", PathTopXiFolder), "READ");
  TFile* pp = TFile::Open(
      Form("%s/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root",
           PathToCombinedFolder),
      "update");
  TIter nextpXi(pXi->GetListOfKeys());
  TKey *keypXi;
  while ((keypXi = (TKey*) nextpXi())) {
    TClass *clpXi = gROOT->GetClass(keypXi->GetClassName());
    if (!clpXi->InheritsFrom("TGraph"))
      continue;
    TGraph *Graph = (TGraph*) keypXi->ReadObj();
    TString graphName = Graph->GetName();
    if (graphName.Contains("Graph")) {
      pp->cd();
      Graph->Write(TString::Format("%s%s", GraphName, Graph->GetName()));
    }
  }
  pp->cd();
  if (sideBand) {
    std::cout << "Sideband Time!!! \n";
    ((TH1F*) pXi->Get("hCk_ReweightedpXiMeV_1"))->Write();
    ((TH1F*) pXi->Get("pXiSidebandDefault"))->Write();
    ((TH1F*) pXi->Get("pXiSidebandDown"))->Write();
    ((TH1F*) pXi->Get("pXiSidebandUp"))->Write();
  }
  pXi->Close();
  pp->Close();

}

int main(int argc, char *argv[]) {
  //argv[1] =  cfpath (without AnalysisResults.root)
  //argv[2] =  prefix
  //argv[3] =  System (0 = pPb, 1 = pp MB, 2 = pp HM)
  //argv[4] =  varfolder
  //argv[5] =  GraphoutName
  //argv[6] =  Path to the combined file
  //argv[7] =  Store Sidebands & CK? > 0
  bool sidebands = false;
  const char* SideBandArg= argv[7];
  TString Sidebandu = Form("%s", SideBandArg);
  std::cout << "you said " << Sidebandu.Data();
  if (Sidebandu != "") {
    std::cout << " so I am doing the Sideband \n";
    sidebands = true;
  } else {
    sidebands = false;
    std::cout << "so I am not doing the Sideband \n";
  }
  int system = atoi(argv[3]);
  TString CalibBaseDir;
  double kStarMin;
  if (system == 0) {
    CalibBaseDir = "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
    kStarMin = 0.008;
  } else if (system == 1) {
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
  } else if (system == 2) {
    CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
    kStarMin = 0.;
  }
  std::cout << "EvalpXiCurves \n";
  EvalpXiCurves(argv[1], argv[2], argv[4],kStarMin);
  std::cout << "EvalError \n";
  EvalError(argv[1], argv[2], argv[4],CalibBaseDir.Data(),kStarMin);
  std::cout << "SidebandCurves \n";
  SidebandCurves(argv[4]);
  std::cout << "CombineIntoOneFile \n";
  CombineIntoOneFile(argv[4], argv[5], argv[6], sidebands);
  return 0;
}
