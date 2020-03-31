#include "VariationAnalysis.h"
#include "TFile.h"
#include "DreamData.h"
#include "DreamPlot.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TApplication.h"

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  gStyle->SetEndErrorSize(5);
  const char* WorkDir = argv[1];
  TApplication app("TheApp", &argc, argv);
  
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
  TCanvas* c1;
  std::vector<const char*> LegNames;
  LegNames.push_back("p-#Omega^{-} #bf{ALICE} data");
  LegNames.push_back("Coulomb + HAL QCD");
  LegNames.push_back("Coulomb + HAL QCD + inelastic");
  LegNames.push_back("Coulomb");
  std::vector<const char*> LegOptions;
  LegOptions.push_back("fpe");
  LegOptions.push_back("f");
  LegOptions.push_back("f");
  LegOptions.push_back("f");

  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = c1->DrawFrame(0, 0.5, 305, 7.5);
  const char * texPtY = "#it{C}(#it{k}*)";
  const char * texPtX = "#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetNdivisions(806);
  TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.08);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.1);
  pad->SetBottomMargin(0.115);
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
  Data->SetCorrelationGraph(cf_graph_stat);
  Data->SetSystematics(cf_graph_syst, 2);
  Data->SetLegendName(LegNames, LegOptions);
  Data->SetDrawAxis(false);
  Data->FemtoModelFitBands(hal5S2Only, kOrange + 1, 10, 0, -4000, true, false);
  Data->FemtoModelFitBands(hal5S2Rad, kGray + 1, 0.5, false);

  Data->FemtoModelFitBands(hal3S15S2Only, kAzure + 1, 10, 0, -4000, true, false);
  Data->FemtoModelFitBands(hal3S15S2Rad, kGray + 1, 0.5, false);

  Data->FemtoModelFitBands(coulomb, kGreen + 1, 1, 3, -4000, true, false);
  

  Data->SetRangePlotting(fXmin, fXmax, 0.6, cf_graph_stat->GetYaxis()->GetXmax() * 1.5);  //ranges
  Data->SetNDivisions(505);

  float legXmin = fTextXMin - 0.02;
  Data->SetLegendCoordinates(0.35, 0.56, 0.65, 0.86);
  Data->SetInletRangePlotting(100, 280, 0.73, 1.31);
  Data->SetInletCoordinates(0.30, 0.18, 0.85, 0.55);
  Data->DrawCorrelationPlot(pad);

  pad->cd();
  cf_graphWidth->SetLineWidth(1);
  cf_graphWidth->SetLineColorAlpha(kBlack, 0.9);
  cf_graphWidth->Draw("same []");
  h->Draw("AXIS SAME");
  
  Data->GetInset()->cd();
  cf_graphWidth->Draw("same []");

  //h->Draw("same"); 
  out->cd();
  c1->Write();
  c1->SaveAs(Form("CF_pOmega.pdf"));
  //SystError->Write();
  out->Close();
  app.Run();
  return 0;

}

