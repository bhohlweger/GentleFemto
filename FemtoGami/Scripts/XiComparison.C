#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TApplication.h"

TColor myColor1;
std::vector<int> fColors = {
  kBlack,         //  0
  kRed + 1,//1
  kBlue + 2,//2
  kGreen + 3,//3
  kMagenta - 8,//4
  kOrange - 7,//5
  kCyan + 2,//6
  kYellow + 2,//7
  kWhite,//8
  kGreen - 5,//9
  myColor1.GetColor(255,127,0),//10
  myColor1.GetColor(31,120,180),//11
  myColor1.GetColor(178,223,138),//12
  kBlue + 3//13
};
std::vector<int> fMarkers = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond,
			     kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

void SetStyleGraph(TGraph *histo, int marker, int color) {
  // histo->GetXaxis()->SetLabelSize(0.045);
  // histo->GetXaxis()->SetTitleSize(0.05);
  // histo->GetXaxis()->SetLabelOffset(0.01);
  // histo->GetXaxis()->SetTitleOffset(1.2);
  // histo->GetXaxis()->SetLabelFont(42);
  // histo->GetYaxis()->SetLabelSize(0.045);
  // histo->GetYaxis()->SetTitleSize(0.05);
  // histo->GetYaxis()->SetLabelOffset(0.01);
  // histo->GetYaxis()->SetTitleOffset(1.0);
  histo->SetMarkerSize(1.4);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(fMarkers[marker]);
  histo->SetMarkerColor(fColors[color]);
  histo->SetLineColor(fColors[color]);
  histo->SetFillColor(fColors[color]);
}

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  gStyle->SetEndErrorSize(5);
  const char* WorkDir = argv[1];
  TApplication app("TheApp", &argc, argv);
  TString cfName = TString::Format("%s/debug_Var0.root", WorkDir).Data();
  TString cfgraphName = TString::Format("%s/CFinput_Var0.root", WorkDir).Data();
  
  TFile* cfFile = TFile::Open(cfName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
						   "InputCF_ResGami_woBL_ResGami_GenuineGami");
  TFile* cfgraphFile = TFile::Open(cfgraphName, "read");
  TGraphAsymmErrors* cf_graph = (TGraphAsymmErrors*) cfgraphFile->FindObjectAny(
										"Graph_from_hCk_RebinnedpXiVar0_0MeV");

  TGraphErrors* coulomb = (TGraphErrors*) cfFile->FindObjectAny("Coulomb");
  TGraphErrors* halRad = (TGraphErrors*) cfFile->FindObjectAny("HalAndRad");
  TGraphErrors* halOnly = (TGraphErrors*) cfFile->FindObjectAny("HalOnly");

  for (int iBin = 1; iBin <= cf_graph->GetN(); ++iBin) {
    double kStar = cf_default->GetBinCenter(iBin);
    double Ck = cf_default->GetBinContent(iBin);
    double x,y;
    cf_graph->GetPoint(iBin-1, x, y);
    cf_graph->SetPoint(iBin-1, x, Ck);
  }

  
  //Yuki Values:
  int iPoint = 0; 
  TGraphErrors* YukiQCD = new TGraphErrors(); 
  YukiQCD->SetPoint(iPoint++,18.4337,	2.59384);
  YukiQCD->SetPoint(iPoint++,18.7952,	2.57127);
  YukiQCD->SetPoint(iPoint++,19.1566,	2.52611);
  YukiQCD->SetPoint(iPoint++,19.8795,	2.48917);
  YukiQCD->SetPoint(iPoint++,20.6024,	2.41323);
  YukiQCD->SetPoint(iPoint++,22.7711,	2.27982);
  YukiQCD->SetPoint(iPoint++,26.0241,	2.08689);
  YukiQCD->SetPoint(iPoint++,28.1928,	2.00068);
  YukiQCD->SetPoint(iPoint++,31.8072,	1.85496);
  YukiQCD->SetPoint(iPoint++,35.7831,	1.74208);
  YukiQCD->SetPoint(iPoint++,43.012 ,	1.56762);
  YukiQCD->SetPoint(iPoint++,49.5181,	1.46294);
  YukiQCD->SetPoint(iPoint++,57.1084,	1.36442);
  YukiQCD->SetPoint(iPoint++,64.3373,	1.28233);
  YukiQCD->SetPoint(iPoint++,71.5663,	1.22486);
  YukiQCD->SetPoint(iPoint++,79.1566,	1.17765);
  YukiQCD->SetPoint(iPoint++,87.4699,	1.13044);
  YukiQCD->SetPoint(iPoint++,96.506 ,	1.09555);
  YukiQCD->SetPoint(iPoint++,103.373,	1.07092);
  YukiQCD->SetPoint(iPoint++,113.133,	1.04835);
  YukiQCD->SetPoint(iPoint++,122.892,	1.03398);
  YukiQCD->SetPoint(iPoint++,135.542,	1.01961);
  YukiQCD->SetPoint(iPoint++,145.301,	1.00935);
  YukiQCD->SetPoint(iPoint++,168.434,	1.00114);
  YukiQCD->SetPoint(iPoint++,229.88 ,	0.99704);
  YukiQCD->SetPoint(iPoint++,299.277,	0.99909);
  
  SetStyleGraph(cf_graph, 2, 0);
  cf_graph->GetXaxis()->SetRangeUser(0, 300);
  cf_graph->GetYaxis()->SetRangeUser(0.9, 3.5); 
  cf_graph->SetMarkerSize(24);
  cf_graph->SetMarkerColor(kBlack);
  cf_graph->SetLineColor(kBlack); 
  cf_graph->SetMarkerSize(1.0);

  SetStyleGraph(halRad, 3, 10);
  halRad->SetFillColorAlpha(fColors[10], 0.7);
  halRad->SetLineColorAlpha(fColors[10], 0.7);
  SetStyleGraph(YukiQCD, 4, 2); 
  TCanvas* c1 = new TCanvas();
  cf_graph->Draw("APEZ");
  halRad->Draw("same L3");
  YukiQCD->Draw("same L");
  cf_graph->Draw("PEZSAME");

  TLegend* leg = new TLegend(0.5,0.5,0.8,0.88);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize());

  leg->AddEntry(cf_graph, "Data", "lep");
  leg->AddEntry(halRad, "HAL CATS", "fl");
  leg->AddEntry(YukiQCD, "HAL YUKI", "l");

  leg->Draw("same");

  c1->SaveAs("pxicomp.pdf");
  
  
  app.Run(); 
  
  
  return 0; 
}
