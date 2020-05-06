#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TApplication.h"

void OmDraw(TPad* OmPad, const char* WorkDir) {
  
  TString cfData = TString::Format("%s/GenOutput.root", WorkDir).Data();
  TString cfModel = TString::Format("%s/bands.root", WorkDir).Data();
 
  TFile* cfgraphFile = TFile::Open(cfData, "read");
  cfgraphFile->ls();
  TGraphErrors* cf_graph_stat_input
    = (TGraphErrors*)cfgraphFile->Get("gGen");
  TGraphAsymmErrors* cf_graph_stat = new TGraphAsymmErrors();
  double x,dx,y;

  for (int iBin = 0; iBin < cf_graph_stat_input->GetN(); ++iBin) {
    cf_graph_stat_input->GetPoint(iBin, x, y);
    double ErrY = cf_graph_stat_input->GetErrorY(iBin); 
    cf_graph_stat->SetPoint(iBin, x, y);
    dx = 5.; 
    if(iBin==0) dx = 4.3;
    if(iBin==1) dx = 5.8;
    if(iBin==2) dx = 5.3;
    
    cf_graph_stat->SetPointError(iBin, dx, dx, ErrY, ErrY); 
  }
  
  std::cout << cf_graph_stat << std::endl; 
  TGraphAsymmErrors* cf_graph_syst
    = (TGraphAsymmErrors*) cfgraphFile->Get("gGensyst");
  std::cout << cf_graph_syst << std::endl;

  double xx,yy;
  //move the first points to the <SE> k* distr mean (see slide 12 of omegamanIV)
  cf_graph_stat->GetPoint(0,xx,yy);
  cf_graph_stat->SetPoint(0,17.2,yy);
  cf_graph_syst->SetPoint(0,17.2,yy);

  cf_graph_stat->GetPoint(1,xx,yy);
  cf_graph_stat->SetPoint(1,36.2,yy);
  cf_graph_syst->SetPoint(1,36.2,yy);
  
  
  TGraphAsymmErrors* cf_graphWidth = new TGraphAsymmErrors(*cf_graph_stat);
  std::cout << cf_graphWidth << std::endl;
  
  //here we assume bins with a width of 20 MeV/c, starting from 5 MeV/c
  double kStar; 
  double binWidth = 10; 
  for (int iBin = 0; iBin < cf_graph_stat->GetN() ; ++iBin) {
    cf_graph_stat->GetPoint(iBin, x, y);
    kStar = 15 + iBin*20;
    //std::cout << "kStar: " << kStar << std::endl;
    double xErrLeft = (x - kStar + binWidth)*0.95;
    double xErrRight = (kStar - x + binWidth)*0.95;
    //std::cout << "x: " << x << " xErrLeft: " << xErrLeft << " xErrRigth: " << xErrRight << std::endl;  
    cf_graphWidth->SetPoint(iBin, x, y);
    cf_graphWidth->SetPointError(iBin, xErrLeft, xErrRight, 0., 0.);
  }

  std::cout << "Done with the Data \n"; 
  
  TFile* ModelFile = TFile::Open(cfModel, "read");
  double xErrHigh, xErrLow, yErrHigh, yErrLow;
  double yMean, yErr; 
  TGraphAsymmErrors* coulomb_input = (TGraphAsymmErrors*) ModelFile->FindObjectAny("CoulombBand");
  TGraphErrors* coulomb = new TGraphErrors();
  coulomb->SetPoint(0, 0.3, 12.);
  coulomb->SetPointError(0, 0, 0.2); 
  for (int iBin = 0; iBin < coulomb_input->GetN(); ++iBin) {
    coulomb_input->GetPoint(iBin, x, y); 
    yErrHigh = y + coulomb_input->GetErrorYhigh(iBin);
    yErrLow = y - coulomb_input->GetErrorYlow(iBin);    
    yMean = (yErrHigh + yErrLow)/2.;
    yErr = yMean - yErrLow;
    //    std::cout << " Coulomb x: " << x << " yMean: " << yMean << " yErr: " << yErr << std::endl; 
    coulomb->SetPoint(iBin+1, x, yMean);
    coulomb->SetPointError(iBin+1, 0, yErr); 
  }
  TGraphAsymmErrors* hal5S2Rad_input = (TGraphAsymmErrors*) ModelFile->FindObjectAny("LatticeBand5S2");
  TGraphErrors* hal5S2Rad = new TGraphErrors();
  for (int iBin = 0; iBin < hal5S2Rad_input->GetN(); ++iBin) {
    hal5S2Rad_input->GetPoint(iBin, x, y); 
    yErrHigh = y + hal5S2Rad_input->GetErrorYhigh(iBin);
    yErrLow = y - hal5S2Rad_input->GetErrorYlow(iBin);    
    yMean = (yErrHigh + yErrLow)/2.;
    yErr = yMean - yErrLow;

    hal5S2Rad->SetPoint(iBin, x, yMean);
    hal5S2Rad->SetPointError(iBin, 0, yErr); 
  }
  
  TGraphAsymmErrors* hal5S2Only_input = (TGraphAsymmErrors*) ModelFile->FindObjectAny("LatticeBandStat5S2");
  TGraphErrors* hal5S2Only = new TGraphErrors();
  for (int iBin = 0; iBin < hal5S2Only_input->GetN(); ++iBin) {
    hal5S2Only_input->GetPoint(iBin, x, y); 
    yErrHigh = y + hal5S2Only_input->GetErrorYhigh(iBin);
    yErrLow = y - hal5S2Only_input->GetErrorYlow(iBin);    
    yMean = (yErrHigh + yErrLow)/2.;
    yErr = yMean - yErrLow;

    hal5S2Only->SetPoint(iBin, x, yMean);
    hal5S2Only->SetPointError(iBin, 0, yErr); 
  }
  TGraphAsymmErrors* hal3S15S2Rad_input = (TGraphAsymmErrors*) ModelFile->FindObjectAny("LatticeBand5S23S1");
  TGraphErrors* hal3S15S2Rad = new TGraphErrors();
  for (int iBin = 0; iBin < hal3S15S2Rad_input->GetN(); ++iBin) {
    hal3S15S2Rad_input->GetPoint(iBin, x, y); 
    yErrHigh = y + hal3S15S2Rad_input->GetErrorYhigh(iBin);
    yErrLow = y - hal3S15S2Rad_input->GetErrorYlow(iBin);    
    yMean = (yErrHigh + yErrLow)/2.;
    yErr = yMean - yErrLow;

    hal3S15S2Rad->SetPoint(iBin, x, yMean);
    hal3S15S2Rad->SetPointError(iBin, 0, yErr); 
  } 
  TGraphAsymmErrors* hal3S15S2Only_input = (TGraphAsymmErrors*) ModelFile->FindObjectAny("LatticeBandStat5S23S1");
  TGraphErrors* hal3S15S2Only = new TGraphErrors();
  for (int iBin = 0; iBin < hal3S15S2Only_input->GetN(); ++iBin) {
    hal3S15S2Only_input->GetPoint(iBin, x, y); 
    yErrHigh = y + hal3S15S2Only_input->GetErrorYhigh(iBin);
    yErrLow = y - hal3S15S2Only_input->GetErrorYlow(iBin);    
    yMean = (yErrHigh + yErrLow)/2.;
    yErr = yMean - yErrLow;

    hal3S15S2Only->SetPoint(iBin, x, yMean);
    hal3S15S2Only->SetPointError(iBin, 0, yErr); 
  }
  //TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("#bf{ALICE} data");
  //LegNames.push_back("Coulomb + HAL QCD elastic");
  //LegNames.push_back("Coulomb + HAL QCD elastic + inelastic");
  //LegNames.push_back("Coulomb");
  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  //LegOptions.push_back("f");
  //LegOptions.push_back("f");
  //LegOptions.push_back("f");

  //c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = OmPad->DrawFrame(0, 0.5, 305, 7.5);
  const char * texPtY = "#it{C}(#it{k}*)";
  const char * texPtX = "#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(2.);
  h->GetYaxis()->SetTitleOffset(2.);
  h->GetXaxis()->SetNdivisions(806);
  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true); 
  //TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) OmPad->cd(0);
  pad->SetRightMargin(0.01);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.);
  pad->SetBottomMargin(0.12);
  pad->Draw();
  pad->cd();
  float fXmin = 0;
  float fXmax = 305;
  float fTextXMin = 0.35;
  float ymaxL = 0.81;
  DreamData *Data = new DreamData(Form("DataOm"));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
  //  Data->SetCorrelationFunction(cf_default);
  Data->SetCorrelationGraph(cf_graph_stat);
  Data->SetSystematics(cf_graph_syst, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(false);
  Data->FemtoModelFitBands(hal5S2Only, kOrange + 1, 10, 0, -4000, false, false);
  Data->FemtoModelFitBands(hal5S2Rad, kGray + 1, 0.5, false);

  Data->FemtoModelFitBands(hal3S15S2Only, kAzure + 1, 10, 0, -4000, false, false);
  Data->FemtoModelFitBands(hal3S15S2Rad, kGray + 1, 0.5, false);

  Data->FemtoModelFitBands(coulomb, kGreen + 1, 1, 3, -4000, false, false);
  

  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_graph_stat->GetYaxis()->GetXmax() * 1.5);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.25, 0.6, 0.55, 0.9, false);
  // Data->SetLegendCoordinates(0.25, 0.56, 0.55, 0.86);
  Data->SetInletRangePlotting(100, 280, 0.73, 1.31);
  Data->SetInletCoordinates(0.12, 0.25, 0.96, 0.8);
  Data->SetAxisOffsetInlet(4.0, 3.0); 
  // Data->SetInletCoordinates(0.35, 0.23, 0.97, 0.64);
  Data->DrawCorrelationPlot(pad);

  pad->cd();
  cf_graphWidth->SetLineWidth(1);
  cf_graphWidth->SetLineColorAlpha(kBlack, 0.9);
  cf_graphWidth->Draw("same []");
  h->Draw("AXIS SAME");
  
  Data->GetInset()->cd();
  cf_graphWidth->Draw("same []");

  OmPad->cd();
  TLatex text;
  text.SetTextFont(43);
  text.SetTextSizePixels(35); 
  text.SetNDC();
  text.SetTextColor(1);
  //  text.SetText
  //text.SetTextSize(gStyle->GetTextSize() * 1.2);
  text.DrawLatex(0.01, 0.95, "#bf{b}");
  text.DrawLatex(0.87, 0.87, "p-#Omega^{-}"); 
  
  
  //h->Draw("same"); 
  //out->cd();
  //c1->Write();
  //c1->SaveAs(Form("CF_pOmega.pdf"));
  //SystError->Write();
  //out->Close();
  //app.Run();
  return;

}

void XiDraw(TPad* XiPad, const char* WorkDir) { 
  //DreamPlot::SetStyle();
  //gStyle->SetEndErrorSize(5);

  //TApplication app("TheApp", &argc, argv);
  TString cfName = TString::Format("%s/debug_Var0.root", WorkDir).Data();
  TString cfgraphName = TString::Format("%s/CFinput_Var0.root", WorkDir).Data();
  TString sysDataName =
    TString::Format("%s/Systematics_pXi.root", WorkDir).Data();
  TString sysNormName = TString::Format("%s/Systematics_pXiNorm.root", WorkDir)
    .Data();
  TString sysLamName = TString::Format("%s/Systematics_pXiLam.root", WorkDir)
    .Data();
  TString sysResName = TString::Format("%s/Systematics_pXiRes.root", WorkDir)
    .Data();

  TFile* cfFile = TFile::Open(cfName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
						   "InputCF_ResGami_woBL_ResGami_GenuineGami");
  TFile* cfgraphFile = TFile::Open(cfgraphName, "read");
  TGraphAsymmErrors* cf_graph = (TGraphAsymmErrors*) cfgraphFile->FindObjectAny(
										"Graph_from_hCk_RebinnedpXiVar0_0MeV");
  if (!cf_default) {
    std::cout << "Default  data not found \n";
    cfFile->ls();
    return;
  }
  cf_default->SetName("DefaultMeV");
  if (!cf_graph) {
    std::cout << "Default  graph not found \n";
    cfFile->ls();
    return;
  }
  TGraphAsymmErrors* cf_graphWidth = new TGraphAsymmErrors(*cf_graph);
  TGraphErrors* coulomb = (TGraphErrors*) cfFile->FindObjectAny("Coulomb");
  TGraphErrors* halRad = (TGraphErrors*) cfFile->FindObjectAny("HalAndRad");
  TGraphErrors* halOnly = (TGraphErrors*) cfFile->FindObjectAny("HalOnly");

  TGraphErrors* pOm_coulomb = (TGraphErrors*)coulomb->Clone("pOmFakeCoulomb");
  TGraphErrors* pOm_hal5S2Only = (TGraphErrors*)coulomb->Clone("pOmFakehal1");
  TGraphErrors* pOm_hal3S15S2Only = (TGraphErrors*)coulomb->Clone("pOmFakehal2");

  //  TGraphErrors* esc = (TGraphErrors*) cfFile->FindObjectAny("ESC");

  TFile* sysDataFile = TFile::Open(sysDataName, "read");
  TF1* systDataErr = (TF1*) sysDataFile->FindObjectAny("SystError");
  if (!systDataErr) {
    std::cout << "systDataErr not found \n";
    sysDataFile->ls();
    return;
  }
  TFile* sysNormFile = TFile::Open(sysNormName, "read");
  TH1F* systNormErr = (TH1F*) sysNormFile->FindObjectAny("SystErrRel");
  if (!systNormErr) {
    std::cout << "systNormErr not found \n";
    sysNormFile->ls();
    return;
  }
  TFile* sysLamFile = TFile::Open(sysLamName, "read");
  TF1* systLamErr = (TF1*) sysLamFile->FindObjectAny("SystError");
  if (!systLamErr) {
    std::cout << "systLamErr not found \n";
    sysLamFile->ls();
    return;
  }
  TFile* sysResFile = TFile::Open(sysResName, "read");
  TH1F* systResErr = (TH1F*) sysResFile->FindObjectAny("SystErrRel");
  if (!systResErr) {
    std::cout << "systResErr not found \n";
    sysResFile->ls();
    return;
  }

  TH1F* SystError = (TH1F*) cf_default->Clone("Sytematics");
  SystError->Reset();

  for (int iBin = 1; iBin <= SystError->GetNbinsX(); ++iBin) {
    double kStar = cf_default->GetBinCenter(iBin);
    double Ck = cf_default->GetBinContent(iBin);
    double x,y;
    cf_graph->GetPoint(iBin-1, x, y);
    cf_graph->SetPoint(iBin-1, x, Ck);

    double xErrLeft = (x - kStar + cf_default->GetBinWidth(iBin) / 2.)*0.95;
    double xErrRight = (kStar - x + cf_default->GetBinWidth(iBin) / 2.)*0.95;
    cf_graphWidth->SetPoint(iBin-1, x, Ck);
    cf_graphWidth->SetPointError(iBin-1, xErrLeft, xErrRight, 0., 0.);

    double errSystData = systDataErr->Eval(kStar);
    double errSystNorm = systNormErr->GetBinContent(iBin);
    double errSystLam = systLamErr->Eval(kStar);
    double errSystMom = systResErr->GetBinContent(iBin);
    double totErr = TMath::Sqrt(
				errSystData * errSystData + errSystNorm * errSystNorm
				+ errSystLam * errSystLam + errSystMom * errSystMom);
    SystError->SetBinContent(iBin, totErr);
  }
  
  //TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("#bf{ALICE} data");
  LegNames.push_back("Coulomb ");
  LegNames.push_back("Coulomb + p-#Xi^{-} HAL QCD");
  LegNames.push_back("Coulomb + p-#Omega^{-} HAL QCD elastic");
  LegNames.push_back("Coulomb + p-#Omega^{-} HAL QCD elastic + inelastic");

  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  LegOptions.push_back("f");
  LegOptions.push_back("f");
  LegOptions.push_back("f");
  LegOptions.push_back("f");

  
  //  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = XiPad->DrawFrame(0, 0.7, 305, 3.7);
  const char * texPtY = "#it{C}(#it{k}*)";
  const char * texPtX = "";//"#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(0);
  //  h->GetXaxis()->CenterTitle(true); 
  h->GetYaxis()->SetTitleOffset(2.);
  h->GetXaxis()->SetNdivisions(806);
  h->GetYaxis()->CenterTitle(true); 
  //TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) XiPad->cd(0);
  pad->SetRightMargin(0.01);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.12);
  pad->SetBottomMargin(0);
  pad->Draw();
  pad->cd();
  float fXmin = 0;
  float fXmax = 305;
  float fTextXMin = 0.35;
  float ymaxL = 0.81;
  DreamData *Data = new DreamData(Form("Data"));
  Data->SetMultiHisto(false);
  Data->SetUnitConversionData(1);
  Data->SetUnitConversionCATS(1);
  //  Data->SetCorrelationFunction(cf_default);
  Data->SetCorrelationGraph(cf_graph);
  Data->SetSystematics(SystError, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(false);
  Data->FemtoModelFitBands(coulomb, kGreen + 1, 1, 3, -4000, true, false);
  Data->FemtoModelFitBands(halOnly, kPink + 1, 10, 0, -4000, true, false);
  Data->FemtoModelFitBands(halRad, kGray + 1, 0.5, false);
  Data->FemtoModelFitBands(pOm_hal3S15S2Only, kOrange + 1, 10, 3, -4000, true, false);
  Data->FemtoModelFitBands(pOm_hal3S15S2Only, kAzure + 1, 10, 0, -4000, true, false);
  

  //  Data->FemtoModelFitBands(esc,     11, 8,  0, -4000, true);
  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_default->GetMaximum() * 1.5);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.25, 0.405, 0.55, 0.78);
  // Data->SetLegendCoordinateos(0.25, 0.635, 0.55, 0.86);
  Data->DrawCorrelationPlot(pad);

  pad->cd();
  cf_graphWidth->SetLineWidth(1);
  cf_graphWidth->SetLineColorAlpha(kBlack, 0.9);
  cf_graphWidth->Draw("same []");
  h->Draw("AXIS SAME");

  TLatex text;
  text.SetTextFont(43);
  text.SetTextSizePixels(35); 
  text.SetNDC();
  text.SetTextColor(1);
  //  text.SetText
  //text.SetTextSize(gStyle->GetTextSize() * 1.2);
  text.DrawLatex(0.01, 0.83, "#bf{a}");
  text.DrawLatex(0.87, 0.75, "p-#Xi^{-}"); 
  
  // out->cd();
  // c1->Write();
  // c1->SaveAs(Form("CF_pXi.pdf"));
  // SystError->Write();
  // out->Close();
  // app.Run();
  return;
}


int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  gStyle->SetEndErrorSize(5);
  const char* WorkDirXi = argv[1];
  const char* WorkDirOm = argv[2];

  TApplication app("TheApp", &argc, argv);
  TCanvas* c1 = new TCanvas("thePlot", "thePlot", 0, 0, 800, 1200);
  c1->cd(); 
  TPad *XiPad = new TPad("XiPad", "XiPad", 0.0, 0.5, 1.0, 1.0);
  XiPad->SetTopMargin(0.01);
  XiPad->SetRightMargin(0.05);
  XiPad->SetBottomMargin(0.);
  XiPad->SetLeftMargin(0.28);
  XiPad->SetFillStyle(4000);
  XiPad->Draw();
  XiPad->cd();
  XiDraw(XiPad, WorkDirXi);

  c1->cd();
  TPad *OmPad = new TPad("OmPad", "OmPad", 0.0, 0.0, 1.0, 0.5);
  OmPad->SetTopMargin(0.01);
  OmPad->SetRightMargin(0.05);
  OmPad->SetBottomMargin(0.28);
  OmPad->SetLeftMargin(0.28);
  OmPad->SetFillStyle(4000);
  OmPad->Draw();
  OmPad->cd();
  OmDraw(OmPad, WorkDirOm);
  c1->SaveAs("Combined.pdf");
  std::cout <<" now just running .... \n" ;
  app.Run();
  return 0; 
}
