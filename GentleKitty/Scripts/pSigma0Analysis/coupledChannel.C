#include "DLM_CkModels.h"
#include "DLM_CkDecomposition.h"
#include "TApplication.h"
#include "TMath.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "DreamPlot.h"
#include <complex>
#include <cmath>
#include "gsl_sf_dawson.h"

/// =====================================================================================
void FillCkGraph(DLM_Ck *ck, TGraph *gr) {
  for (unsigned int i = 0; i < ck->GetNbins(); ++i) {
    const float mom = ck->GetBinCenter(0, i);
    gr->SetPoint(i, mom, ck->Eval(mom));
  }
}

/// =====================================================================================
double CoupledLednicky(const double &Momentum, const double* SourcePar,
                       const double* PotPar) {
  //This model tries to calculate the Sigma0 correlation function, one can simply say trying because nothing is known about
  //this kind of interaction

  //Parameter definitions:
  //*****************************************
  //*****************************************
  //             isospin 1/2 (==0)
  //________________________________________
  const TComplex scatt_length_1s0_0(PotPar[0], 0);
  const TComplex effrange_1s0_0(PotPar[1], 0);

  const double hbarc = 0.197326971844;
  //Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
  const double MomentumGeV = Momentum / 1000.;

  //
  //             isospin 3/2 (==1)
  //________________________________________
  const TComplex scatt_length_1s0_1(PotPar[0], 0.);
  const TComplex effrange_1s0_1(PotPar[1], 0.);

  const Double_t mass_proton = 0.938272;
  const Double_t mass_neutron = 0.939565;
  const Double_t mass_sigmaplus = 1.189377;
  const Double_t mass_sigma0 = 1.192642;

  Double_t mu1 = mass_proton * mass_sigma0 / (mass_proton + mass_sigma0);
  Double_t mu2 = mass_neutron * mass_sigmaplus
      / (mass_neutron + mass_sigmaplus);
  //*****************************************
  //*****************************************

  Double_t k1 = MomentumGeV / hbarc;      //momentum in the elastic channel
  Double_t k2 = TMath::Sqrt(
      2 * mu2
          * (MomentumGeV * MomentumGeV / (2 * mu1) + mass_proton + mass_sigma0
              - mass_neutron - mass_sigmaplus));
  k2 /= hbarc;  //momentum in the inelastic channel
  //(k2>k1)

  Double_t RG = SourcePar[0];  //Gaussian radius

  //Define inverse K-matrices for different spin configurations in isopin basis:

  //spin singlet:
  TComplex Kmin_1s0_0 = 1. / scatt_length_1s0_0
      + 0.5 * effrange_1s0_0 * k1 * k1;
  TComplex Kmin_1s0_1 = 1. / scatt_length_1s0_1
      + 0.5 * effrange_1s0_1 * k1 * k1;

  //transform to transition basis with help of Clebsch-Gordan coefficients:
  TComplex Kmin_11_1s0 = 2. / 3. * Kmin_1s0_1 + 1. / 3. * Kmin_1s0_0;

  TComplex Kmin_22_1s0 = 1. / 3. * Kmin_1s0_1 + 2. / 3. * Kmin_1s0_0;

  TComplex Kmin_12_1s0 = TMath::Sqrt(2.) / 3. * (Kmin_1s0_1 - Kmin_1s0_0);

  //Determinant of scattering amplitude matrix:
  TComplex D_1s0 = (Kmin_11_1s0 + TComplex(0., -k1))
      * (Kmin_22_1s0 + TComplex(0., -k2)) - TComplex::Power(Kmin_12_1s0, 2.);

  //What is really needed are the scattering amplitudes:

  TComplex f11_1s0 = (Kmin_22_1s0 + TComplex(0., -k2)) / D_1s0;

  TComplex f12_1s0 = -(Kmin_12_1s0) / D_1s0;

  Double_t fF1 = gsl_sf_dawson(2 * k1 * RG) / (2 * k1 * RG);
  Double_t fF2 = (1. - TMath::Exp(-pow(2 * k1 * RG, 2))) / (2 * k1 * RG);

  //elastic part of the correlation function(1 -> 1):
  //*********************************************
  Double_t corr_1s0 = 0.5 * f11_1s0.Rho2() / (RG * RG)
      + 2. * f11_1s0.Re() / (TMath::Sqrt(TMath::Pi()) * RG) * fF1
      - f11_1s0.Im() / RG * fF2;

  //*********************************************

  //inelastic part of the correlation function(2 -> 1):
  //*********************************************
  Double_t corr_1s0_inel = 0.5 * mu2 / mu1 * fabs(f12_1s0.Rho2() / (RG * RG));

  //*********************************************

  //correction term to non spheric distortions on the short range scale
  //Define matrices d_ij = 2Re d(K_ij)/dk²

  TComplex dKmin_11_1s0_dk = 2. / 3. * 0.5 * effrange_1s0_1
      + 1. / 3. * 0.5 * effrange_1s0_0;

  TComplex dKmin_22_1s0_dk = 1. / 3. * 0.5 * effrange_1s0_1
      + 2. / 3. * 0.5 * effrange_1s0_0;

  TComplex dKmin_12_1s0_dk = TMath::Sqrt(2.) / 3. * 0.5
      * (effrange_1s0_1 - effrange_1s0_0);

  TComplex factor1_1s0 = f11_1s0 * TComplex::Conjugate(f12_1s0);

  Double_t corr_1s0_distort = f11_1s0.Rho2() * 2. * dKmin_11_1s0_dk.Re()
      + f12_1s0.Rho2() * 2. * dKmin_22_1s0_dk.Re()
      + 2. * factor1_1s0.Re() * 2. * dKmin_12_1s0_dk.Re();
  corr_1s0_distort *= -1.f
      / (4 * TMath::Sqrt(TMath::Pi()) * TMath::Power(RG, 3.));

  std::cout << corr_1s0_inel << "\n";

  Double_t corr_fin = 1. + corr_1s0 + corr_1s0_inel + corr_1s0_distort;
  return corr_fin;

}

/// =====================================================================================
int main(int argc, char *argv[]) {
  TApplication *app = new TApplication("app", 0, 0);

  DreamPlot::SetStyle();

  const int momBins = 100;
  const double kmin = 0;
  const double kmax = 400;
  const double radius = 1.25;
  const int lineWidth = 3;

  const double scatLen = atof(argv[1]);
  const double effRange = atof(argv[2]);

  auto grnormalCk = new TGraph();
  DreamPlot::SetStyleGraph(grnormalCk, 20, kBlue + 3, 0.8);
  grnormalCk->SetLineWidth(lineWidth);

  auto normalCk = new DLM_Ck(1, 2, momBins, kmin, kmax,
                             Lednicky_Singlet);
  normalCk->SetPotPar(0, scatLen);
  normalCk->SetPotPar(1, effRange);
  normalCk->SetSourcePar(0, radius);
  normalCk->Update();
  FillCkGraph(normalCk, grnormalCk);

  auto grcoupledCk = new TGraph();
  DreamPlot::SetStyleGraph(grcoupledCk, 20, kGreen + 3, 0.8);
  grcoupledCk->SetLineWidth(lineWidth);

  auto coupledCk = new DLM_Ck(1, 2, momBins, kmin, kmax, CoupledLednicky);
  coupledCk->SetPotPar(0, scatLen);
  coupledCk->SetPotPar(1, effRange);
  coupledCk->SetSourcePar(0, radius);
  coupledCk->Update();
  FillCkGraph(coupledCk, grcoupledCk);

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            kmax);
  DreamPlot::SetStyleHisto(dummyHist, 20, kWhite);

  auto c = new TCanvas("CFComp", "CFComp", 0, 0, 650, 550);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.025);
  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0., 5);
  dummyHist->GetXaxis()->SetNdivisions(504);
  grnormalCk->Draw("L3");
  grcoupledCk->Draw("L3");

  app->Run();

}
