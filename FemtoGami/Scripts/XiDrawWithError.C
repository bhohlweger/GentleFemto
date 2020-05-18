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

    //TApplication app("TheApp", &argc, argv);
  TString cfName = TString::Format("%s/debug_Var0.root", WorkDir).Data();
  TString modelName = TString::Format("%s/CFsXi.root", WorkDir).Data();
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
  TFile* modelFile = TFile::Open(modelName, "read");
  TH1F* cf_default = (TH1F*) cfFile->FindObjectAny(
						   "InputCF_ResGami_woBL_ResGami_GenuineGami");
  TFile* cfgraphFile = TFile::Open(cfgraphName, "read");
  TGraphAsymmErrors* cf_graph = (TGraphAsymmErrors*) cfgraphFile->FindObjectAny(
										"Graph_from_hCk_RebinnedpXiVar0_0MeV");
  if (!cf_default) {
    std::cout << "cf_default InputCF_ResGami_woBL_ResGami_GenuineGami not found \n";
    cfFile->ls();
    return 0;
  }
  cf_default->SetName("DefaultMeV");
  if (!cf_graph) {
    std::cout << "Default  graph not found \n";
    cfFile->ls();
    return 0;
  }
  TGraphAsymmErrors* cf_graphWidth = new TGraphAsymmErrors(*cf_graph);
  TGraphErrors* coulomb = (TGraphErrors*) cfFile->FindObjectAny("Coulomb");
  TGraphErrors* halRad = (TGraphErrors*) modelFile->FindObjectAny("gr_HALRadSymSum");
  TGraphErrors* halOnly = (TGraphErrors*) modelFile->FindObjectAny("HAL_SymSum");
  if (!coulomb) {
    std::cout << "Coulomb missing \n";
    cfFile->ls();
  }
  if (!halRad) {
    std::cout << "halRad missing \n";
    modelFile->ls();
  }
  if (!halOnly) {
    std::cout << "Halonly missing\n";
    modelFile->ls();
  }
  TGraphErrors* pOm_coulomb = (TGraphErrors*)coulomb->Clone("pOmFakeCoulomb");
  TGraphErrors* pOm_hal5S2Only = (TGraphErrors*)coulomb->Clone("pOmFakehal1");
  TGraphErrors* pOm_hal3S15S2Only = (TGraphErrors*)coulomb->Clone("pOmFakehal2");

//  TGraphErrors* esc = (TGraphErrors*) cfFile->FindObjectAny("ESC");

  TFile* sysDataFile = TFile::Open(sysDataName, "read");
  TF1* systDataErr = (TF1*) sysDataFile->FindObjectAny("SystError");
  if (!systDataErr) {
    std::cout << "systDataErr not found \n";
    sysDataFile->ls();
    return 0;
  }
  TFile* sysNormFile = TFile::Open(sysNormName, "read");
  TH1F* systNormErr = (TH1F*) sysNormFile->FindObjectAny("SystErrRel");
  if (!systNormErr) {
    std::cout << "systNormErr not found \n";
    sysNormFile->ls();
    return 0;
  }
  TFile* sysLamFile = TFile::Open(sysLamName, "read");
  TF1* systLamErr = (TF1*) sysLamFile->FindObjectAny("SystError");
  if (!systLamErr) {
    std::cout << "systLamErr not found \n";
    sysLamFile->ls();
    return 0;
  }
  TFile* sysResFile = TFile::Open(sysResName, "read");
  TH1F* systResErr = (TH1F*) sysResFile->FindObjectAny("SystErrRel");
  if (!systResErr) {
    std::cout << "systResErr not found \n";
    sysResFile->ls();
    return 0;
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
  //calculate nSigma Values
  double chiSqCoulomb;
  double chiSqHAL;
  int ndf = 0;
  
  double ks = 0;
  double Ck = 0; 
  double kCoulomb = 0; 
  double CkCoulomb =0;
  double kHAL = 0; 
  double CkHAL =0;
  cf_graph->GetPoint(ndf, ks, Ck);
  while (ks < 133) { 
    //total uncertainty in this bin
    cf_graph->GetPoint(ndf, ks, Ck);
    //Systematic uncertainty
    double systErr = SystError->GetBinContent(ndf+1)*Ck; 
    //Statistical uncertainty data
    double statErr = cf_graph->GetErrorY(ndf); 
    //uncertainty model
    int binModel = TMath::Nint(ks) - 1;
    coulomb->GetPoint(binModel, kCoulomb, CkCoulomb);
    if (TMath::Abs(TMath::Nint(ks)-kCoulomb) > 1e-3) {
      std::cout << "Something is off with the binning... kStarData = " << ks << " kStarCoulomb = " << kCoulomb << std::endl;
      return -999; 
    }
    double CoulombErr = coulomb->GetErrorY(binModel); 

    halRad->GetPoint(binModel, kHAL, CkHAL); 

    if (TMath::Abs(TMath::Nint(ks)-kHAL) > 1e-3) {
      std::cout << "Something is off with the binning... kStarData = " << ks << " kStarHAL = " << kHAL << std::endl;
      return -999; 
    }

    double HALErr = halRad->GetErrorY(binModel); 

    //total uncertainty in this bin

    double TotUncertCoulomb = TMath::Sqrt(systErr*systErr+statErr*statErr+CoulombErr*CoulombErr);
    double TotUncertHAL = TMath::Sqrt(systErr*systErr+statErr*statErr+HALErr*HALErr);
    //Add sfuff to chisq
    chiSqCoulomb+= (Ck-CkCoulomb)*(Ck-CkCoulomb)/(TotUncertCoulomb*TotUncertCoulomb); 
    chiSqHAL+= (Ck-CkHAL)*(Ck-CkHAL)/(TotUncertHAL*TotUncertHAL);

    std::cout << "ks = " << ks << " Ck = " << Ck << " stat Ck= " << statErr << " syst CK = " <<  systErr
	      << " Coulomb  = " << CkCoulomb << " Coulomb Err = "<< CoulombErr
	      << " ChisqCoulomb = " << (Ck-CkCoulomb) << " Coulomb tot = " << TotUncertCoulomb 
	      << " HAL = " << CkHAL << " HAL Err = " << HALErr
	      << " ChisqHAL = " << (Ck-CkHAL) << " HAL tot = " << TotUncertHAL << std::endl; 
    ndf++; 
  }

  double pvalCoulomb = TMath::Prob(chiSqCoulomb, round(ndf));
  double nSigmaCoulomb = TMath::Sqrt(2) * TMath::ErfcInverse(pvalCoulomb);

  double pvalHAL = TMath::Prob(chiSqHAL, round(ndf));
  double nSigmaHAL = TMath::Sqrt(2) * TMath::ErfcInverse(pvalHAL);

  std::cout << "=============================================\n";
  std::cout << "==============SIGMA VALUES BABY==============\n";
  std::cout << "=============================================\n";
  std::cout << "////////// nSigma Coulomb = " << nSigmaCoulomb << "//////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "//////////// nSigma HAL = " << nSigmaHAL << "////////////"
            << std::endl;
  std::cout << "=============================================\n";
  std::cout << "=============================================\n";
  
  TCanvas* c1;
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

  
  c1 = new TCanvas("c2", "c2", 0, 0, 800, 600);
  TH1 * h = c1->DrawFrame(0, 0.7, 305, 3.7);
  const char * texPtY = "#it{C}(#it{k}*)";
  const char * texPtX = "#it{k}* (MeV/#it{c})";
  h->SetXTitle(texPtX);
  h->SetYTitle(texPtY);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetXaxis()->CenterTitle(true); 
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetNdivisions(806);
  h->GetYaxis()->CenterTitle(true); 
  TFile* out = TFile::Open(Form("out.root"), "recreate");
  TPad* pad;
  float LatexX = 0.;
  pad = (TPad*) c1->cd(0);
  pad->SetRightMargin(0.01);
  pad->SetLeftMargin(0.1);
  pad->SetTopMargin(0.01);
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
  Data->SetLegendCoordinates(0.25, 0.525, 0.55, 0.9);
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
  text.DrawLatex(0.01, 0.95, "#bf{a}");
  text.DrawLatex(0.87, 0.87, "p-#Xi^{-}"); 
  
  out->cd();
  c1->Write();
  c1->SaveAs(Form("CF_pXi.pdf"));
  SystError->Write();
  out->Close();
  app.Run();
}
