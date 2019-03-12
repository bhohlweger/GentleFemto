//#include <XiOnly.C>
#include <iostream>
#include <fstream>

#include "CATStools.h"
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TFile.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

void RUN2_SYSTEMATICS_MEDIAN(const char* InputFolder, int Numiter,
                             const char* OutDirName, const int system) {
  const char* FileName = Form("%s/OutFile_CutVarAdd_Iter%u.root", InputFolder,
                              Numiter);
  gStyle->SetOptStat(0);

  TFile* inFile = TFile::Open(FileName, "READ");
  TNtuple* sysVarTree = (TNtuple*) inFile->Get("ntResult");
  if (!sysVarTree) {
    std::cout << "no Tree loaded\n";
  }
  const float histRangeLow = (system == 0) ? 1.29 : 1.1;
  const float histRangeUp = (system == 0) ? 1.35 : 1.3;
  std::cout << "histRangeLow: " << histRangeLow << '\t' << "histRangeUp: "
            << histRangeUp << std::endl;
//  const int nBins = (system == 0) ? 1000 : 500;
  const int nBins = 100;
  auto histRad = new TH1D("hRad", "hRad", nBins, histRangeLow, histRangeUp);
  auto histRadIter = new TH2D("hRadIter", "hRadIter", nBins, histRangeLow,
                              histRangeUp, 360, 0, 360);
  sysVarTree->Draw("Radius_pp>>hRad");
  auto mean = histRad->GetMean();

  int n = histRad->GetXaxis()->GetNbins();
  std::vector<double> x(n);
  histRad->GetXaxis()->GetCenter(&x[0]);
  const double * y = histRad->GetArray();
  // exclude underflow/overflows from bin content array y
  auto median = TMath::Median(n, &x[0], &y[1]);
  auto medianBin = histRad->FindBin(median);
  std::cout << "Median: " << median << "\t Mean: " << mean << std::endl;

  auto histRadCumulative = histRad->GetCumulative();
  auto canRad = new TCanvas("cRad", "cRad", 1200, 800);
  canRad->Divide(3, 3);
  canRad->cd(1);
  histRad->SetTitle("; r_{0} (fm); Entries");
  histRad->Draw();
  canRad->cd(2);
  histRadCumulative->Scale(1 / (double) histRad->GetEntries());
  histRadCumulative->Draw();

  int binMin = 0;
  int binMax = 0;
  for (int iBin = 0; iBin < histRadCumulative->GetNbinsX(); iBin++) {
    if (binMin == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) - 0.32)) {
      binMin = iBin;
    }
    if (binMax == 0
        && (histRadCumulative->GetBinContent(iBin)
            > histRadCumulative->GetBinContent(medianBin) + 0.32)) {
      binMax = iBin;
      break;
    }
  }
  auto radMin = histRadCumulative->GetXaxis()->GetBinCenter(binMin);
  auto radMax = histRadCumulative->GetXaxis()->GetBinCenter(binMax);
  std::cout << "RadMin: " << radMin << "\tRadMax: " << radMax << '\n';
  canRad->cd(3);
  sysVarTree->Draw("IterID:Radius_pp>>hRadIter");
  canRad->cd(4);
  auto hIterMean = histRadIter->ProjectionY(
      "hIterMean", histRadIter->GetXaxis()->FindBin(mean),
      histRadIter->GetXaxis()->FindBin(mean));
  hIterMean->Draw();
  int uIterMean;
  hIterMean->GetBinWithContent(1., uIterMean);
  std::cout << "uIterMean: " << uIterMean << std::endl;

  canRad->cd(5);
  auto hIterMedian = histRadIter->ProjectionY(
      "hIterMedian", histRadIter->GetXaxis()->FindBin(median),
      histRadIter->GetXaxis()->FindBin(median));
  hIterMedian->Draw();
  int uIterMedian;
  hIterMedian->GetBinWithContent(1., uIterMedian);
  std::cout << "uIterMedian: " << uIterMedian << std::endl;

  canRad->cd(6);
  auto hIterUp = histRadIter->ProjectionY(
      "hIterUp", histRadIter->GetXaxis()->FindBin(radMax),
      histRadIter->GetXaxis()->FindBin(radMax));
  hIterUp->Draw();
  int uIterUp;
  hIterUp->GetBinWithContent(1., uIterUp);
  std::cout << "uIterUp: " << uIterUp << std::endl;

  canRad->cd(7);
  auto hIterLow = histRadIter->ProjectionY(
      "hIterLow", histRadIter->GetXaxis()->FindBin(radMin),
      histRadIter->GetXaxis()->FindBin(radMin));
  hIterLow->Draw();
  int uIterLow;
  hIterLow->GetBinWithContent(1., uIterLow);
  std::cout << "uIterLow: " << uIterLow << std::endl;

  canRad->SaveAs(Form("%s/canRad.pdf", OutDirName));

  TFile* outFile = TFile::Open(
      Form("%s/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root", OutDirName),
      "RECREATE");

  TFile* fileMean = TFile::Open(
      Form("%s/GraphFile_CutVarAdd_Iter%u_uIter%i.root", InputFolder, Numiter,
           uIterMean));
  if (!fileMean) {
    std::cout << "Missing file Mean \n";
    return;
  }
  TGraph* FileGrDefault_pp = (TGraph*) fileMean->Get(
      Form("FitResult_pp_%u", uIterMean));
//  auto CkDef = (TH1F*)fileMean->Get("");
  TGraph GrDefault_pp(*FileGrDefault_pp);
  GrDefault_pp.SetName("ppGraphDefault");
  outFile->cd();
  GrDefault_pp.Write();
  fileMean->Close();

  TFile* fileSysLow = TFile::Open(
      Form("%s/GraphFile_CutVarAdd_Iter%u_uIter%i.root", InputFolder, Numiter,
           uIterLow));
  if (!fileSysLow) {
    std::cout << "Missing file Sys Low \n";
    return;
  }
  TGraph* FileGrLow_pp = (TGraph*) fileSysLow->Get(
      Form("FitResult_pp_%u", uIterLow));
  TGraph GrLow_pp(*FileGrLow_pp);
  GrLow_pp.SetName("ppGraphLowerLim");
  outFile->cd();
  GrLow_pp.Write();
  fileSysLow->Close();

  TFile* fileSysUp = TFile::Open(
      Form("%s/GraphFile_CutVarAdd_Iter%u_uIter%i.root", InputFolder, Numiter,
           uIterUp));
  if (!fileSysUp) {
    std::cout << "Missing file SysUp \n";
    return;
  }
  TGraph* FileGrUp_pp = (TGraph*) fileSysUp->Get(
      Form("FitResult_pp_%u", uIterUp));
  TGraph GrUp_pp(*FileGrUp_pp);
  GrUp_pp.SetName("ppGraphUpperLim");
  outFile->cd();
  GrUp_pp.Write();
  fileSysUp->Close();

  float uIterIDDefault;
  float vFemReg_pp;
  float vModpL;
  float vFrac_pp_pL;
  float vFrac_pL_pSigma0;
  float vFrac_pL_pXim;
  float HaveWeABaseLine;

  float rDefault_pp;
  float rErr_pp;
  float pa_pp, pb_pp;
  float iNorm;

  sysVarTree->SetBranchAddress("IterID", &uIterIDDefault);
  sysVarTree->SetBranchAddress("vFemReg_pp", &vFemReg_pp);
  sysVarTree->SetBranchAddress("vModpL", &vModpL);
  sysVarTree->SetBranchAddress("vFrac_pp_pL", &vFrac_pp_pL);
  sysVarTree->SetBranchAddress("vFrac_pL_pSigma0", &vFrac_pL_pSigma0);
  sysVarTree->SetBranchAddress("vFrac_pL_pXim", &vFrac_pL_pXim);
  sysVarTree->SetBranchAddress("BLSlope", &HaveWeABaseLine);
  sysVarTree->SetBranchAddress("IterID", &uIterIDDefault);
  sysVarTree->SetBranchAddress("Radius_pp", &rDefault_pp);
  sysVarTree->SetBranchAddress("RadiusErr_pp", &rErr_pp);
  sysVarTree->SetBranchAddress("pa_pp", &pa_pp);
  sysVarTree->SetBranchAddress("pb_pp", &pb_pp);
  sysVarTree->SetBranchAddress("pb_pp", &pb_pp);
//  sysVarTree->SetBranchAddress("iNorm", &iNorm);

  for (int iEntry = 0; iEntry < sysVarTree->GetEntries(); iEntry++) {
    sysVarTree->GetEntry(iEntry);
    if (vFemReg_pp == 1 && vFrac_pp_pL == 1 && vFrac_pL_pSigma0 == 1
        && vFrac_pL_pXim == 1 && vModpL == 2
        && HaveWeABaseLine == (int) false) {
      break;
    }
  }
  auto errLow = rDefault_pp - radMin;
  auto errUp = radMax - rDefault_pp;

  const double BL_a = pa_pp;
  const double BL_b = pb_pp;

  float rLower = rDefault_pp
      - TMath::Sqrt(
          rErr_pp * rErr_pp + (0.2 * rDefault_pp) * (0.2 * rDefault_pp)
              + errLow * errLow);
  float rUp = rDefault_pp + TMath::Sqrt(rErr_pp * rErr_pp + errLow * errLow);
  std::cout << "Default radius\t" << rDefault_pp << "\nStat. Error\t" << rErr_pp
            << "\nSyst. Error low\t" << errLow << "\nSyst. Error up\t" << errUp
            << "\nLower radius\t" << rLower << "\nUpper radius\t" << rUp
            << '\n';
  TNtuple* outTuple = new TNtuple(
      "outTuple", "outTuple",
      "Rad_pp:RadStat_pp:RadSysLow_pp:RadSysUp_pp:RadLow_pXi:RadUp_pXi");
  float outArray[6];
  outArray[0] = rDefault_pp;
  outArray[1] = rErr_pp;
  outArray[2] = errLow;
  outArray[3] = errUp;
  outArray[4] = rLower;
  outArray[5] = rUp;
  outTuple->Fill(outArray);
  outFile->cd();
  outTuple->Write();

  std::ofstream radiusOut;
  radiusOut.open(Form("%s/radius.dat", OutDirName));
  radiusOut << rDefault_pp << " " << rErr_pp << " " << errLow << " " << errUp
            << "\n";
  radiusOut.close();

  std::ofstream baseline;
  baseline.open(Form("%s/baseline.dat", OutDirName));
  baseline << BL_a << " " << BL_b << "\n";
  baseline.close();

  auto canRad2 = new TCanvas();
  histRad->Rebin(2);
  histRad->Draw();

  auto histRadLimits = (TH1F*) histRad->Clone("histRadLimits");
  histRadLimits->Reset();
  for (int i = 0; i < histRad->GetNbinsX(); ++i) {
    if (histRad->GetBinCenter(i) < radMin || histRad->GetBinCenter(i) > radMax)
      continue;
    histRadLimits->SetBinContent(i, histRad->GetBinContent(i));
  }
  histRadLimits->SetFillColor(kGray + 1);
  histRadLimits->Draw("same");

  auto lineDefault = new TLine(rDefault_pp, 0, rDefault_pp,
                               histRad->GetMaximum());
  lineDefault->SetLineColor(kRed + 2);
  lineDefault->SetLineWidth(2);
  lineDefault->Draw("same");

  auto lineLow = new TLine(radMin, 0, radMin, histRad->GetMaximum());
  lineLow->SetLineColor(kGreen + 2);
  lineLow->SetLineWidth(2);
  lineLow->Draw("same");

  auto lineUp = new TLine(radMax, 0, radMax, histRad->GetMaximum());
  lineUp->SetLineColor(kGreen + 2);
  lineUp->SetLineWidth(2);
  lineUp->Draw("same");

  canRad2->SaveAs(Form("%s/radius.pdf", OutDirName));

//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rLower, 11,
//			outFile, "UpperLim", true);
//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rDefault_pp, 12,
//			outFile, "Default", true);
//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rUp, 13, outFile,
//			"LowerLim", true);

  outFile->Close();
}

int main(int argc, char *argv[]) {
  RUN2_SYSTEMATICS_MEDIAN(argv[1], atoi(argv[2]), argv[3], atoi(argv[4]));
  return 0;
}
