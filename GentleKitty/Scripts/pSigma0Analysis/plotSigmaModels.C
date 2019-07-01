#include "TAxis.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "gsl_sf_dawson.h"
#include "DLM_Histo.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "TGraph.h"
#include "DLM_WfModel.h"
#include "DLM_CkModels.h"
#include "DreamPlot.h"
#include "CATSInputSigma0.h"
#include "DLM_CkDecomposition.h"
#include "CATSLambdaParam.h"
#include "TDatabasePDG.h"
#include "TidyCats.h"
#include "TStyle.h"

/// =====================================================================================
double Lednicky_gauss_Sigma0(const double& Momentum, const double* SourcePar,
                             double& singlet, double& triplet) {
  // This model tries to calculate the Sigma0 correlation function, one can
  // simply say trying because nothing is known about this kind of interaction

  // Parameter definitions:
  //*****************************************
  //*****************************************
  //             isospin 1/2 (==0)
  //________________________________________
  const TComplex scatt_length_1s0_0(-1.1, 0.);
  const TComplex scatt_length_3s1_0(-1.1, 4.3);
  const TComplex effrange_1s0_0(-1.5, 0.);
  const TComplex effrange_3s1_0(-2.2, -2.4);

  const double hbarc = 0.197326971844;
  // Oli's original code was in GeV and fm though, so we must convert the
  // Momentum to GeV!
  const double MomentumGeV = Momentum / 1000.;

  //
  //             isospin 3/2 (==1)
  //________________________________________
  const TComplex scatt_length_1s0_1(2.51, 0.);
  const TComplex scatt_length_3s1_1(-0.73, 0.);
  const TComplex effrange_1s0_1(4.92, 0.);
  const TComplex effrange_3s1_1(-1.22, 0.);

  const Double_t mass_proton = 0.938272;
  const Double_t mass_neutron = 0.939565;
  const Double_t mass_sigmaplus = 1.189377;
  const Double_t mass_sigma0 = 1.192642;

  Double_t mu1 = mass_proton * mass_sigma0 / (mass_proton + mass_sigma0);
  Double_t mu2 = mass_neutron * mass_sigmaplus
      / (mass_neutron + mass_sigmaplus);
  //*****************************************
  //*****************************************

  Double_t k1 = MomentumGeV / hbarc;  // momentum in the elastic channel
  Double_t k2 = TMath::Sqrt(
      2 * mu2
          * (MomentumGeV * MomentumGeV / (2 * mu1) + mass_proton + mass_sigma0
              - mass_neutron - mass_sigmaplus));
  k2 /= hbarc;  // momentum in the inelastic channel
  //(k2>k1)

  Double_t RG = SourcePar[0];  // Gaussian radius

  // Define inverse K-matrices for different spin configurations in isopin
  // basis:

  // spin singlet:
  TComplex Kmin_1s0_0 = 1. / scatt_length_1s0_0
      + 0.5 * effrange_1s0_0 * k1 * k1;
  TComplex Kmin_1s0_1 = 1. / scatt_length_1s0_1
      + 0.5 * effrange_1s0_1 * k1 * k1;
  // spin triplet:
  TComplex Kmin_3s1_0 = 1. / scatt_length_3s1_0
      + 0.5 * effrange_3s1_0 * k1 * k1;
  TComplex Kmin_3s1_1 = 1. / scatt_length_3s1_1
      + 0.5 * effrange_3s1_1 * k1 * k1;

  // transform to transition basis with help of Clebsch-Gordan coefficients:
  TComplex Kmin_11_1s0 = 2. / 3. * Kmin_1s0_1 + 1. / 3. * Kmin_1s0_0;
  TComplex Kmin_11_3s1 = 2. / 3. * Kmin_3s1_1 + 1. / 3. * Kmin_3s1_0;

  TComplex Kmin_22_1s0 = 1. / 3. * Kmin_1s0_1 + 2. / 3. * Kmin_1s0_0;
  TComplex Kmin_22_3s1 = 1. / 3. * Kmin_3s1_1 + 2. / 3. * Kmin_3s1_0;

  TComplex Kmin_12_1s0 = TMath::Sqrt(2.) / 3. * (Kmin_1s0_1 - Kmin_1s0_0);
  TComplex Kmin_12_3s1 = TMath::Sqrt(2.) / 3. * (Kmin_3s1_1 - Kmin_3s1_0);

  // Determinant of scattering amplitude matrix:
  TComplex D_1s0 = (Kmin_11_1s0 + TComplex(0., -k1))
      * (Kmin_22_1s0 + TComplex(0., -k2)) - TComplex::Power(Kmin_12_1s0, 2.);
  TComplex D_3s1 = (Kmin_11_3s1 + TComplex(0., -k1))
      * (Kmin_22_3s1 + TComplex(0., -k2)) - TComplex::Power(Kmin_12_3s1, 2.);

  // What is really needed are the scattering amplitudes:

  TComplex f11_1s0 = (Kmin_22_1s0 + TComplex(0., -k2)) / D_1s0;
  TComplex f11_3s1 = (Kmin_22_3s1 + TComplex(0., -k2)) / D_3s1;

  TComplex f12_1s0 = -(Kmin_12_1s0) / D_1s0;
  TComplex f12_3s1 = -(Kmin_12_3s1) / D_3s1;

  Double_t fF1 = gsl_sf_dawson(2 * k1 * RG) / (2 * k1 * RG);
  Double_t fF2 = (1. - TMath::Exp(-pow(2 * k1 * RG, 2))) / (2 * k1 * RG);

  // elastic part of the correlation function(1 -> 1):
  //*********************************************
  Double_t corr_1s0 = 0.5 * f11_1s0.Rho2() / (RG * RG)
      + 2. * f11_1s0.Re() / (TMath::Sqrt(TMath::Pi()) * RG) * fF1
      - f11_1s0.Im() / RG * fF2;
  corr_1s0 *= 0.25;

  Double_t corr_3s1 = 0.5 * f11_3s1.Rho2() / (RG * RG)
      + 2. * f11_3s1.Re() / (TMath::Sqrt(TMath::Pi()) * RG) * fF1
      - f11_3s1.Im() / RG * fF2;
  corr_3s1 *= 0.75;
  //*********************************************

  // inelastic part of the correlation function(2 -> 1):
  //*********************************************
  Double_t corr_1s0_inel = 0.5 * mu2 / mu1 * fabs(f12_1s0.Rho2() / (RG * RG));
  corr_1s0_inel *= 0.25;

  Double_t corr_3s1_inel = 0.5 * mu2 / mu1 * fabs(f12_3s1.Rho2() / (RG * RG));
  corr_3s1_inel *= 0.75;

  //*********************************************

  // correction term to non spheric distortions on the short range scale
  // Define matrices d_ij = 2Re d(K_ij)/dk²

  TComplex dKmin_11_1s0_dk = 2. / 3. * 0.5 * effrange_1s0_1
      + 1. / 3. * 0.5 * effrange_1s0_0;
  TComplex dKmin_11_3s1_dk = 2. / 3. * 0.5 * effrange_3s1_1
      + 1. / 3. * 0.5 * effrange_3s1_0;

  TComplex dKmin_22_1s0_dk = 1. / 3. * 0.5 * effrange_1s0_1
      + 2. / 3. * 0.5 * effrange_1s0_0;
  TComplex dKmin_22_3s1_dk = 1. / 3. * 0.5 * effrange_3s1_1
      + 2. / 3. * 0.5 * effrange_3s1_0;

  TComplex dKmin_12_1s0_dk = TMath::Sqrt(2.) / 3. * 0.5
      * (effrange_1s0_1 - effrange_1s0_0);
  TComplex dKmin_12_3s1_dk = TMath::Sqrt(2.) / 3. * 0.5
      * (effrange_3s1_1 - effrange_3s1_0);

  TComplex factor1_1s0 = f11_1s0 * TComplex::Conjugate(f12_1s0);
  TComplex factor1_3s1 = f11_3s1 * TComplex::Conjugate(f12_3s1);

  Double_t corr_1s0_distort = f11_1s0.Rho2() * 2. * dKmin_11_1s0_dk.Re()
      + f12_1s0.Rho2() * 2. * dKmin_22_1s0_dk.Re()
      + 2. * factor1_1s0.Re() * 2. * dKmin_12_1s0_dk.Re();
  corr_1s0_distort *= -0.25
      / (4 * TMath::Sqrt(TMath::Pi()) * TMath::Power(RG, 3.));

  Double_t corr_3s1_distort = f11_3s1.Rho2() * 2. * dKmin_11_3s1_dk.Re()
      + f12_3s1.Rho2() * 2. * dKmin_22_3s1_dk.Re()
      + 2. * factor1_3s1 * 2. * dKmin_12_3s1_dk.Re();
  corr_3s1_distort *= -0.75
      / (4 * TMath::Sqrt(TMath::Pi()) * TMath::Power(RG, 3.));

  Double_t corr_fin = 1. + corr_1s0 + corr_3s1 + corr_1s0_inel + corr_3s1_inel
      + corr_1s0_distort + corr_3s1_distort;

  const float scaleSinglet = 1. / 0.25;
  const float scaleTriplet = 1. / 0.75;

  singlet = 1 + scaleSinglet * corr_1s0 + scaleSinglet * corr_1s0_inel
      + scaleSinglet * corr_1s0_distort;
  triplet = 1 + scaleTriplet * corr_3s1 + scaleTriplet * corr_3s1_inel
      + scaleTriplet * corr_3s1_distort;

  return corr_fin;
}

/// =====================================================================================
void SourcePlay() {
  const double rCore = 0.748428;

  const double massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  const double massSigma0 = TDatabasePDG::Instance()->GetParticle(3212)->Mass()
      * 1000;
  double massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;

  DLM_CleverMcLevyReso* CleverMcLevyReso = new DLM_CleverMcLevyReso();
  CleverMcLevyReso->InitStability(1, 2 - 1e-6, 2 + 1e-6);
  CleverMcLevyReso->InitScale(35, 0.25, 2.0);
  CleverMcLevyReso->InitRad(512, 0, 64);
  CleverMcLevyReso->InitType(2);
  CleverMcLevyReso->InitReso(0, 1);
  CleverMcLevyReso->InitReso(1, 1);
  CleverMcLevyReso->SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, massProton,
                              massPion);
  CleverMcLevyReso->SetUpReso(1, 0, 1. - 0.3735, 1581.73, 4.28, massSigma0,
                              massPion);
  CleverMcLevyReso->InitNumMcIter(1000000);

  CATS cats;
  cats.SetAnaSource(CatsSourceForwarder, CleverMcLevyReso, 2);
  cats.SetAnaSource(0, rCore);  // this is the radius for this mT
  cats.SetAnaSource(1, 2.0);

  auto grSource = new TGraph();
  DreamPlot::SetStyleGraph(grSource, 20, kBlue + 3);
  grSource->SetLineWidth(2);
  grSource->SetTitle(";#it{r} (fm); S(#it{r}) (fm^{-1})");

  for (double i = 0; i < 150; ++i) {
    grSource->SetPoint(i, i * 0.1, cats.EvaluateTheSource(0, i * 0.1, 0));
  }

  auto gaussFit =
      new TF1(
          "gaus",
          [&](double *x, double *p) {return
            4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
            std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
          },
          0, 10, 1);
  gaussFit->SetParameter(0, 1.3);
  gaussFit->SetNpx(1000);
  gaussFit->SetLineColor(kBlue + 2);
  gaussFit->SetLineWidth(2);
  gaussFit->SetLineStyle(2);

  auto c = new TCanvas();
  grSource->Draw("AL");
  grSource->GetXaxis()->SetRangeUser(0, 12);
  grSource->Fit(gaussFit, "", "RQ", 0, 6);
  auto leg = new TLegend(0.4, 0.63, 0.85, 0.85);
  leg->SetTextFont(42);
  leg->AddEntry(grSource,
                "p#minus#Sigma^{0} #LT #it{m}_{T} #GT = 2.07 GeV/#it{c}^{2}",
                "l");
  leg->AddEntry((TObject*) nullptr, Form("#it{r}_{core} = %.3f fm", rCore), "");
  leg->AddEntry(
      gaussFit,
      Form("Gauss fit #it{r}_{G, eff} = %.3f fm", gaussFit->GetParameter(0)),
      "l");
  leg->Draw("same");

  c->Print("ResonanceSource.pdf");
}

/// =====================================================================================
void plotProtonLambda(double* radius, float momBins, float kmin, float kmax) {

  auto grTotalLednicky = new TGraph();
  grTotalLednicky->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalLednicky, 20, kBlack);
  grTotalLednicky->SetLineWidth(2);
  auto grTotalHaidenbauer = new TGraph();
  grTotalHaidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalHaidenbauer, 20, kBlack);
  grTotalHaidenbauer->SetLineWidth(2);
  auto grTotalESC16 = new TGraph();
  grTotalESC16->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalESC16, 20, kBlack);
  grTotalESC16->SetLineWidth(2);
  auto grTotalNSC97f = new TGraph();
  grTotalNSC97f->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalNSC97f, 20, kBlack);
  grTotalNSC97f->SetLineWidth(2);

  auto TotalHaidenbauer = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                  Lednicky_SingletTriplet);
  TotalHaidenbauer->SetPotPar(0, 2.91);
  TotalHaidenbauer->SetPotPar(1, 2.78);
  TotalHaidenbauer->SetPotPar(2, 1.54);
  TotalHaidenbauer->SetPotPar(3, 2.72);
  TotalHaidenbauer->SetSourcePar(0, radius[0]);
  TotalHaidenbauer->Update();
  for (unsigned int i = 0; i < TotalHaidenbauer->GetNbins(); ++i) {
    const float mom = TotalHaidenbauer->GetBinCenter(0, i);
    grTotalHaidenbauer->SetPoint(i, mom, TotalHaidenbauer->Eval(mom));
  }

  auto TotalLednicky = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                  Lednicky_SingletTriplet);
  TotalLednicky->SetPotPar(0, 2.59);
  TotalLednicky->SetPotPar(1, 2.83);
  TotalLednicky->SetPotPar(2, 1.60);
  TotalLednicky->SetPotPar(3, 3.01);
  TotalLednicky->SetSourcePar(0, radius[0]);
  TotalLednicky->Update();
  for (unsigned int i = 0; i < TotalLednicky->GetNbins(); ++i) {
    const float mom = TotalLednicky->GetBinCenter(0, i);
    grTotalLednicky->SetPoint(i, mom, TotalLednicky->Eval(mom));
  }

  auto TotalESC16 = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                  Lednicky_SingletTriplet);
  TotalESC16->SetPotPar(0, 1.88);
  TotalESC16->SetPotPar(1, 3.58);
  TotalESC16->SetPotPar(2, 1.86);
  TotalESC16->SetPotPar(3, 3.37);
  TotalESC16->SetSourcePar(0, radius[0]);
  TotalESC16->Update();
  for (unsigned int i = 0; i < TotalESC16->GetNbins(); ++i) {
    const float mom = TotalESC16->GetBinCenter(0, i);
    grTotalESC16->SetPoint(i, mom, TotalESC16->Eval(mom));
  }

  auto TotalNSC97f = new DLM_Ck(1, 4, momBins, kmin, kmax,
                                Lednicky_SingletTriplet);
  TotalNSC97f->SetPotPar(0, 2.51);
  TotalNSC97f->SetPotPar(1, 3.03);
  TotalNSC97f->SetPotPar(2, 1.75);
  TotalNSC97f->SetPotPar(3, 3.32);
  TotalNSC97f->SetSourcePar(0, radius[0]);
  TotalNSC97f->Update();
  for (unsigned int i = 0; i < TotalNSC97f->GetNbins(); ++i) {
    const float mom = TotalNSC97f->GetBinCenter(0, i);
    grTotalNSC97f->SetPoint(i, mom, TotalNSC97f->Eval(mom));
  }

  auto compareCF = new TCanvas();
  grTotalESC16->SetLineColor(kGreen + 2);
  grTotalHaidenbauer->SetLineColor(kAzure);
  grTotalLednicky->SetLineColor(kRed + 1 );
  grTotalNSC97f->SetLineColor(kOrange + 2);
  grTotalHaidenbauer->Draw("AL");
  grTotalHaidenbauer->GetXaxis()->SetRangeUser(0, kmax);
  grTotalHaidenbauer->GetYaxis()->SetRangeUser(0.6, 3.5);
  grTotalESC16->Draw("lsame");
  grTotalLednicky->Draw("lsame");
  grTotalNSC97f->Draw("lsame");
  auto leg4 = new TLegend(0.65, 0.52, 0.5, 0.84);
  leg4->SetBorderSize(0);
  leg4->SetTextFont(42);
  leg4->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg4->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg4->AddEntry(grTotalLednicky, "fss2", "l");
  leg4->AddEntry(grTotalHaidenbauer, "#chiEFT (NLO)", "l");
  leg4->AddEntry(grTotalESC16, "ESC16", "l");
  leg4->AddEntry(grTotalNSC97f, "NSC97f", "l");
  leg4->Draw("same");
  compareCF->Print("allCF_protonLambda.pdf");

  delete compareCF;
}

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TGraph *gr, std::vector<float> weights) {
  int channel = 0;
  for (auto &it : weights) {
    kitty.SetChannelWeight(channel++, it);

  }
  kitty.KillTheCat();
  FillWaveGraph(kitty, gr);
}

/// =====================================================================================
int main(int argc, char* argv[]) {
  DreamPlot::SetStyle();
  double* radius = new double[1];
  radius[0] = 1.154;
  int momBins = 40;
  int kmin = -4.99;
  int kmax = 395.01;

  SourcePlay();
  plotProtonLambda(radius, momBins, kmin, kmax);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Calibration
  TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
  CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
  CATSinput->ReadResFile();
  CATSinput->SetSigmaFileName("Sample6_MeV_compact.root");
  CATSinput->ReadSigmaFile();

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Haidenbauer
  TidyCats* tidy = new TidyCats();
  CATS Kitty;
  tidy->GetCatsProtonSigma0(&Kitty, momBins, kmin, kmax, TidyCats::sGaussian,
                            TidyCats::pSigma0Haidenbauer);
  Kitty.KillTheCat();
  auto Ck_Haidenbauer = new DLM_Ck(1, 0, Kitty);
  Ck_Haidenbauer->SetSourcePar(0, radius[0]);
  Ck_Haidenbauer->Update();

  // ESC16
  CATS KittyESC16;
  tidy->GetCatsProtonSigma0(&KittyESC16, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0ESC16);
  KittyESC16.KillTheCat();
  auto Ck_ESC16 = new DLM_Ck(1, 0, KittyESC16);
  Ck_ESC16->SetSourcePar(0, radius[0]);
  Ck_ESC16->Update();

  // NSC97f
  CATS KittyNSC97f;
  tidy->GetCatsProtonSigma0(&KittyNSC97f, momBins, kmin, kmax,
                            TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
  KittyNSC97f.KillTheCat();
  auto Ck_NSC97f = new DLM_Ck(1, 0, KittyNSC97f);
  Ck_NSC97f->SetSourcePar(0, radius[0]);
  Ck_NSC97f->Update();

  // Haidenbauer with resonances
  CATS ResonantKitty;
  tidy->GetCatsProtonSigma0(&ResonantKitty, momBins, kmin, kmax,
                            TidyCats::sResonance, TidyCats::pSigma0Haidenbauer);
  ResonantKitty.KillTheCat();
  auto Ck_ResonantHaidenbauer = new DLM_Ck(1, 0, ResonantKitty);
  Ck_ResonantHaidenbauer->SetSourcePar(0, 0.72);
  Ck_ResonantHaidenbauer->Update();

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Smearing
  auto grTotalHaidenbauerSmeared = new TGraph();
  grTotalHaidenbauerSmeared->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalHaidenbauerSmeared, 20, kGreen + 2);
  grTotalHaidenbauerSmeared->SetLineWidth(2);
  grTotalHaidenbauerSmeared->SetLineStyle(2);
  auto grTotalLednickySmeared = new TGraph();
  grTotalLednickySmeared->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalLednickySmeared, 20, kGreen + 2);
  grTotalLednickySmeared->SetLineWidth(2);
  grTotalLednickySmeared->SetLineStyle(2);

  DLM_CkDecomposition CkDec_Haidenbauer("pSigma0Haidenbauer", 1,
                                        *Ck_Haidenbauer,
                                        CATSinput->GetSigmaFile(1));

  DLM_Ck* Ck_Lednicky = new DLM_Ck(1, 0, momBins, kmin, kmax,
                                   Lednicky_gauss_Sigma0);
  Ck_Lednicky->SetSourcePar(0, radius[0]);
  Ck_Lednicky->Update();
  DLM_CkDecomposition CkDec_Lednicky("pSigma0Lednicky", 1, *Ck_Lednicky,
                                     CATSinput->GetSigmaFile(1));

  for (unsigned int i = 0; i < Ck_Haidenbauer->GetNbins(); ++i) {
    const float mom = Ck_Haidenbauer->GetBinCenter(0, i);
    grTotalHaidenbauerSmeared->SetPoint(i, mom, CkDec_Haidenbauer.EvalCk(mom));
    grTotalLednickySmeared->SetPoint(i, mom, CkDec_Lednicky.EvalCk(mom));
  }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Lambda parameters
  auto grTotalHaidenbauerLambda = new TGraph();
  grTotalHaidenbauerLambda->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalHaidenbauerLambda, 20, kAzure + 2);
  grTotalHaidenbauerLambda->SetLineWidth(2);
  grTotalHaidenbauerLambda->SetLineStyle(2);
  auto grTotalLednickyLambda = new TGraph();
  grTotalLednickyLambda->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalLednickyLambda, 20, kAzure + 2);
  grTotalLednickyLambda->SetLineWidth(2);
  grTotalLednickyLambda->SetLineStyle(2);
  const double protonPurity = 0.991213;
  const double protonPrimary = 0.874808;
  const double protonLambda = 0.0876342;

  // proton secondary contribution systematic variations
  const double protonSecondary = protonLambda / (1. - protonPrimary);

  const double sigmaPurity = 0.289;
  const double sigmaPrimary = 1.;
  const Particle sigma0(sigmaPurity, sigmaPrimary, { { 0 } });

  const Particle proton(
      protonPurity,
      protonPrimary,
      { { (1. - protonPrimary) * protonSecondary, (1. - protonPrimary)
          * (1 - protonSecondary) } });

  CATSLambdaParam lambdaParams(proton, sigma0);

  CkDec_Haidenbauer.AddContribution(
      0,
      lambdaParams.GetLambdaParam(CATSLambdaParam::Primary,
                                  CATSLambdaParam::Fake, 0, 0),
      DLM_CkDecomposition::cFake);
  CkDec_Haidenbauer.Update();

  CkDec_Lednicky.AddContribution(
      0,
      lambdaParams.GetLambdaParam(CATSLambdaParam::Primary,
                                  CATSLambdaParam::Fake, 0, 0),
      DLM_CkDecomposition::cFake);
  CkDec_Lednicky.Update();
  for (unsigned int i = 0; i < Ck_Haidenbauer->GetNbins(); ++i) {
    const float mom = Ck_Haidenbauer->GetBinCenter(0, i);
    grTotalHaidenbauerLambda->SetPoint(i, mom, CkDec_Haidenbauer.EvalCk(mom));
    grTotalLednickyLambda->SetPoint(i, mom, CkDec_Lednicky.EvalCk(mom));
  }

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Correlation functions

  auto grTotalLednicky = new TGraph();
  grTotalLednicky->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalLednicky, 20, kBlack);
  grTotalLednicky->SetLineWidth(2);
  auto grSingletLednicky = new TGraph();
  grSingletLednicky->SetLineWidth(2);
  grSingletLednicky->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSingletLednicky, 20, kRed + 2);
  auto grTripletLednicky = new TGraph();
  grTripletLednicky->SetLineWidth(2);
  grTripletLednicky->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grTripletLednicky, 20, kBlue + 2);
  auto grTotalHaidenbauer = new TGraph();
  grTotalHaidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  grTotalHaidenbauer->SetLineWidth(2);
  DreamPlot::SetStyleGraph(grTotalHaidenbauer, 20, kBlack);
  auto grTotalResonantHaidenbauer = new TGraph();
  grTotalResonantHaidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  grTotalResonantHaidenbauer->SetLineWidth(2);
  DreamPlot::SetStyleGraph(grTotalResonantHaidenbauer, 20, kGreen + 2);
  auto grSingletHaidenbauer = new TGraph();
  grSingletHaidenbauer->SetLineWidth(2);
  grSingletHaidenbauer->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSingletHaidenbauer, 20, kRed + 2);
  auto grTripletHaidenbauer = new TGraph();
  grTripletHaidenbauer->SetLineWidth(2);
  grTripletHaidenbauer->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grTripletHaidenbauer, 20, kBlue + 2);

  auto grpS0Haidenbauer = new TGraph();
  grpS0Haidenbauer->SetLineWidth(2);
  grpS0Haidenbauer->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grpS0Haidenbauer, 20, kBlack);
  auto grpS0pLHaidenbauer = new TGraph();
  grpS0pLHaidenbauer->SetLineWidth(2);
  grpS0pLHaidenbauer->SetLineStyle(3);
  DreamPlot::SetStyleGraph(grpS0pLHaidenbauer, 20, kGreen + 2);
  auto grpS0nSplusHaidenbauer = new TGraph();
  grpS0nSplusHaidenbauer->SetLineWidth(2);
  grpS0nSplusHaidenbauer->SetLineStyle(3);
  DreamPlot::SetStyleGraph(grpS0nSplusHaidenbauer, 20, kBlue + 2);
  auto grpLambdaHaidenbauer = new TGraph();
  grpLambdaHaidenbauer->SetLineWidth(2);
  grpLambdaHaidenbauer->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grpLambdaHaidenbauer, 20, kGreen + 2);
  auto grnSPlusHaidenbauer = new TGraph();
  grnSPlusHaidenbauer->SetLineWidth(2);
  grnSPlusHaidenbauer->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grnSPlusHaidenbauer, 20, kBlue + 2);


  auto grTotalESC16 = new TGraph();
  grTotalESC16->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  DreamPlot::SetStyleGraph(grTotalESC16, 20, kBlack);
  grTotalESC16->SetLineWidth(2);
  auto grSingletESC16 = new TGraph();
  grSingletESC16->SetLineWidth(2);
  grSingletESC16->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSingletESC16, 20, kRed + 2);
  auto grTripletESC16 = new TGraph();
  grTripletESC16->SetLineWidth(2);
  grTripletESC16->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grTripletESC16, 20, kBlue + 2);

  auto grTotalNSC97f = new TGraph();
  grTotalNSC97f->SetTitle(";#it{k}* (MeV/#it{c}); C(#it{k}*)");
  grTotalNSC97f->SetLineWidth(2);
  DreamPlot::SetStyleGraph(grTotalNSC97f, 20, kBlack);
  auto grSingletNSC97f = new TGraph();
  grSingletNSC97f->SetLineWidth(2);
  grSingletNSC97f->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grSingletNSC97f, 20, kRed + 2);
  auto grTripletNSC97f= new TGraph();
  grTripletNSC97f->SetLineWidth(2);
  grTripletNSC97f->SetLineStyle(2);
  DreamPlot::SetStyleGraph(grTripletNSC97f, 20, kBlue + 2);

  double cf_singletLednicky, cf_tripletLednicky;
  for (unsigned int i = 0; i < Ck_Lednicky->GetNbins(); ++i) {
    const float mom = Ck_Lednicky->GetBinCenter(0, i);
    grTotalLednicky->SetPoint(i, mom, Ck_Lednicky->Eval(mom));
    Lednicky_gauss_Sigma0(mom, radius, cf_singletLednicky, cf_tripletLednicky);
    grSingletLednicky->SetPoint(i, mom, cf_singletLednicky);
    grTripletLednicky->SetPoint(i, mom, cf_tripletLednicky);
  }

  FillWaveGraph(Kitty, grTotalHaidenbauer);
  FillWaveGraph(ResonantKitty, grTotalResonantHaidenbauer);

  FillWaveGraph(Kitty, grSingletHaidenbauer, { 1, 0, 0, 0, 0, 0 });
  FillWaveGraph(Kitty, grTripletHaidenbauer, { 0, 1, 0, 0, 0, 0 });
  FillWaveGraph(Kitty, grpS0Haidenbauer, { 0.25, 0.75, 0, 0, 0, 0 });
  FillWaveGraph(Kitty, grpLambdaHaidenbauer, { 0, 0, 0.25, 0.75, 0, 0 });
  FillWaveGraph(Kitty, grnSPlusHaidenbauer, { 0, 0, 0, 0, 0.25, 0.75 });
  FillWaveGraph(Kitty, grpS0pLHaidenbauer, { 0.25, 0.75, 0.25, 0.75, 0, 0 });
  FillWaveGraph(Kitty, grpS0nSplusHaidenbauer, { 0.25, 0.75, 0.25, 0.75, 0.25, 0.75 });

  FillWaveGraph(KittyESC16, grTotalESC16);
  FillWaveGraph(KittyESC16, grSingletESC16, { 1, 0 });
  FillWaveGraph(KittyESC16, grTripletESC16, { 0, 1 });

  FillWaveGraph(KittyNSC97f, grTotalNSC97f);
  FillWaveGraph(KittyNSC97f, grSingletNSC97f, { 1, 0 });
  FillWaveGraph(KittyNSC97f, grTripletNSC97f, { 0, 1 });


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Plotting

  const float ylow = 0.25;
  const float yup = 3;

  auto c = new TCanvas();
  grTotalLednicky->Draw("AL");
  grTotalLednicky->GetXaxis()->SetRangeUser(0, kmax);
  grTotalLednicky->GetYaxis()->SetRangeUser(ylow, yup);
  grSingletLednicky->Draw("lsame");
  grTripletLednicky->Draw("lsame");
  auto leg = new TLegend(0.6, 0.6, 0.89, 0.84);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg->AddEntry(grTotalLednicky, "Total", "l");
  leg->AddEntry(grSingletLednicky, "Singlet", "l");
  leg->AddEntry(grTripletLednicky, "Triplet", "l");
  leg->Draw("same");
  c->Print("CF_Lednicky_coupled.pdf");

  auto d = new TCanvas();
  grTotalHaidenbauer->Draw("AL");
  grTotalHaidenbauer->GetXaxis()->SetRangeUser(0, kmax);
  grTotalHaidenbauer->GetYaxis()->SetRangeUser(ylow, yup);
  grSingletHaidenbauer->Draw("lsame");
  grTripletHaidenbauer->Draw("lsame");
  leg->Draw("same");
  d->Print("CF_Haidenbauer.pdf");

  auto e = new TCanvas();
  grTotalLednicky->Draw("AL");
  grTotalLednicky->GetXaxis()->SetRangeUser(0, kmax);
  grTotalLednicky->GetYaxis()->SetRangeUser(ylow, yup);
  grTotalLednickySmeared->Draw("lsame");
  grTotalLednickyLambda->Draw("lsame");
  auto leg2 = new TLegend(0.5, 0.6, 0.89, 0.84);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(42);
  leg2->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg2->AddEntry(grTotalLednicky, "Model", "l");
  leg2->AddEntry(grTotalLednickySmeared, "Smeared", "l");
  leg2->AddEntry(grTotalLednickyLambda, "Lambda param", "l");
  leg2->Draw("same");
  e->Print("CF_Lednicky_coupled_smeared.pdf");

  auto f = new TCanvas();
  grTotalHaidenbauer->Draw("AL");
  grTotalHaidenbauer->GetXaxis()->SetRangeUser(0, kmax);
  grTotalHaidenbauer->GetYaxis()->SetRangeUser(ylow, yup);
  grTotalHaidenbauerSmeared->Draw("lsame");
  grTotalHaidenbauerLambda->Draw("lsame");
  leg2->Draw("same");
  f->Print("CF_Haidenbauer_smeared.pdf");

  auto g = new TCanvas();
  grTotalHaidenbauer->Draw("AL");
  grTotalResonantHaidenbauer->Draw("Lsame");
  auto leg3 = new TLegend(0.45, 0.72, 0.84, 0.84);
  leg3->SetBorderSize(0);
  leg3->SetTextFont(42);
  leg3->AddEntry(grTotalHaidenbauer,
                 Form("Gaussian source #it{r}_{0} = %.3f fm", radius[0]), "l");
  leg3->AddEntry(grTotalResonantHaidenbauer,
                 "p#minus#Sigma^{0} #LT #it{m}_{T} #GT = 2.07 GeV/#it{c}^{2}",
                 "l");
  leg3->Draw("same");
  g->Print("CF_Haidenbauer_resonances.pdf");

  auto h = new TCanvas();
  grTotalESC16->Draw("AL");
  grTotalESC16->GetXaxis()->SetRangeUser(0, kmax);
  grTotalESC16->GetYaxis()->SetRangeUser(ylow, yup);
  grSingletESC16->Draw("lsame");
  grTripletESC16->Draw("lsame");
  leg->Draw("same");
  h->Print("CF_ESC16.pdf");

  auto i = new TCanvas();
  grTotalNSC97f->Draw("AL");
  grTotalNSC97f->GetXaxis()->SetRangeUser(0, kmax);
  grTotalNSC97f->GetYaxis()->SetRangeUser(ylow, yup);
  grSingletNSC97f->Draw("lsame");
  grTripletNSC97f->Draw("lsame");
  leg->Draw("same");
  i->Print("CF_NSC97f.pdf");

  auto compareCF = new TCanvas();
  grTotalESC16->SetLineColor(kGreen + 2);
  grTotalHaidenbauer->SetLineColor(kAzure);
  grTotalLednicky->SetLineColor(kRed + 1 );
  grTotalNSC97f->SetLineColor(kOrange + 2 );
  grTotalESC16->Draw("AL");
  grTotalESC16->GetXaxis()->SetRangeUser(0, kmax);
  grTotalESC16->GetYaxis()->SetRangeUser(0.6, 2.25);
  grTotalHaidenbauer->Draw("lsame");
  grTotalLednicky->Draw("lsame");
  grTotalNSC97f->Draw("lsame");
  auto leg4 = new TLegend(0.65, 0.52, 0.5, 0.84);
  leg4->SetBorderSize(0);
  leg4->SetTextFont(42);
  leg4->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg4->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg4->AddEntry(grTotalLednicky, "fss2 (Lednicky)", "l");
  leg4->AddEntry(grTotalHaidenbauer, "#chiEFT (NLO)", "l");
  leg4->AddEntry(grTotalESC16, "ESC16", "l");
  leg4->AddEntry(grTotalNSC97f, "NSC97f", "l");
  leg4->Draw("same");
  compareCF->Print("allCF.pdf");

  auto HaidenbauerDecomp = new TCanvas();
  grpS0Haidenbauer->Draw("AL");
  grpS0Haidenbauer->GetXaxis()->SetRangeUser(0, 350);
  grpS0Haidenbauer->GetYaxis()->SetRangeUser(0, yup);
  grpLambdaHaidenbauer->Draw("Lsame");
  grnSPlusHaidenbauer->Draw("Lsame");
  grpS0pLHaidenbauer->Draw("Lsame");
  grTotalHaidenbauer->SetLineColor(kRed + 1 );
  grTotalHaidenbauer->Draw("lsame");
  grpS0nSplusHaidenbauer->Draw("Lsame");
  auto leg5 = new TLegend(0.65, 0.44, 0.5, 0.84);
  leg5->SetBorderSize(0);
  leg5->SetTextFont(42);
  leg5->SetTextSize(gStyle->GetTextSize() * 0.9);
  leg5->SetHeader(Form("#it{r}_{0} = %.3f fm", radius[0]));
  leg5->AddEntry(grpS0Haidenbauer, "p#minus#Sigma^{0} (no CC)", "l");
  leg5->AddEntry(grpLambdaHaidenbauer, "p#minus#Lambda (CC)", "l");
  leg5->AddEntry(grnSPlusHaidenbauer, "n#minus#Sigma^{+} (CC)", "l");
  leg5->AddEntry(grpS0pLHaidenbauer, "p#minus#Sigma^{0} + p#minus#Lambda (CC)", "l");
  leg5->AddEntry(grpS0nSplusHaidenbauer, "p#minus#Sigma^{0} +n#minus#Sigma^{+} (CC)", "l");
  leg5->AddEntry(grTotalHaidenbauer, "p#minus#Sigma^{0} (with CC)", "l");
  leg5->Draw("same");
  HaidenbauerDecomp->Print("CF_HaidenbauerDecomp.pdf");

}
