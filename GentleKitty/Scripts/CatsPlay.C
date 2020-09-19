#include "PlayWithCats.h"
#include "TidyCats.cxx"
#include "DLM_CkModels.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
int main(int argc, char *argv[]) {
  TFile* out = TFile::Open("pXiSource.root", "recreate"); 
  // TidyCats* tidy = new TidyCats();

  // CATS pS0_ESC16;
  // CATS pS0_NSC97;
  // CATS pS0_Haidenbro;
  // CATS pXi;
  
  // double Radius = 1.2;

  // int NumMomBins = 30;
  // double kMin = 0.5;
  // double kMax = 300.5; 

  
  // tidy->GetCatsProtonSigma0(&pS0_ESC16, NumMomBins, kMin, kMax, TidyCats::sGaussian, TidyCats::pSigma0ESC16);
  // pS0_ESC16.SetAnaSource(0,1.2); 
  // pS0_ESC16.KillTheCat();

  // TGraph* pS0_ESC16_1_2 = new TGraph();
  // pS0_ESC16_1_2->SetName("pS0_ESC16_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pS0_ESC16.GetMomentum(it);
  //   double mean = pS0_ESC16.GetCorrFun(it);
  //   pS0_ESC16_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pS0_ESC16_1_2->Write("pS0_ESC16_1_2");
  
  // tidy->GetCatsProtonSigma0(&pS0_NSC97, NumMomBins, kMin, kMax, TidyCats::sGaussian, TidyCats::pSigma0NSC97f);
  // pS0_NSC97.SetAnaSource(0,1.2); 
  // pS0_NSC97.KillTheCat();

  // TGraph* pS0_NSC97f_1_2 = new TGraph();
  // pS0_NSC97f_1_2->SetName("pS0_NSC97f_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pS0_NSC97.GetMomentum(it);
  //   double mean = pS0_NSC97.GetCorrFun(it);
  //   pS0_NSC97f_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pS0_NSC97f_1_2->Write("pS0_NSC97f_1_2");
  
  // tidy->GetCatsProtonSigma0(&pS0_Haidenbro, NumMomBins, kMin, kMax, TidyCats::sGaussian, TidyCats::pSigma0Haidenbauer);
  // pS0_Haidenbro.SetAnaSource(0,1.2); 
  // pS0_Haidenbro.KillTheCat();

  // TGraph* pS0_Haidenbro_1_2 = new TGraph();
  // pS0_Haidenbro_1_2->SetName("pS0_Haidenbro_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pS0_Haidenbro.GetMomentum(it);
  //   double mean = pS0_Haidenbro.GetCorrFun(it);
  //   pS0_Haidenbro_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pS0_Haidenbro_1_2->Write("pS0_Haidenbro_1_2");

  // TGraph* pS0_Ledni_1_2 = new TGraph();
  // pS0_Ledni_1_2->SetName("pS0_Ledni_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pS0_Haidenbro.GetMomentum(it);
  //   double mean = Lednicky_gauss_Sigma0(kStar,&Radius, nullptr); 
  //   pS0_Ledni_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pS0_Ledni_1_2->Write("pS0_Ledni_1_2");
  
  // tidy->GetCatsProtonXiMinus(&pXi, NumMomBins, kMin, kMax, TidyCats::sGaussian, TidyCats::pHALQCD, 24);
  // pXi.SetAnaSource(0, 1.2);
  // pXi.KillTheCat();

  // TGraph* pXi_1_2 = new TGraph();
  // pXi_1_2->SetName("pXi_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pXi.GetMomentum(it);
  //   double mean = pXi.GetCorrFun(it);
  //   pXi_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pXi_1_2->Write("pXi_1_2");
  
  // out->Close(); 


  // out->cd();
 
  // tidy->GetCatsProtonLambda(&pL_LO, NumMomBins, kMin, kMax, TheSource,
  // 			      TidyCats::pLOWF);
  // pL_LO.SetAnaSource(0, 1.2);
  // pL_LO.KillTheCat();

  // TGraph* pL_LO_1_2 = new TGraph();
  // pL_LO_1_2->SetName("pL_LO_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pL_LO.GetMomentum(it);
  //   double mean = pL_LO.GetCorrFun(it);
  //   pL_LO_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pL_LO_1_2->Write("pL_LO_1_2");
  
  // tidy->GetCatsProtonLambda(&pL_NLO, NumMomBins, kMin, kMax, TheSource,
  // 			    TidyCats::pNLOWF);
  // pL_NLO.SetAnaSource(0, 1.2);
  // pL_NLO.KillTheCat();

  // TGraph* pL_NLO_1_2 = new TGraph();
  // pL_NLO_1_2->SetName("pL_NLO_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pL_NLO.GetMomentum(it);
  //   double mean = pL_NLO.GetCorrFun(it);
  //   pL_NLO_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pL_NLO_1_2->Write("pL_NLO_1_2");

  // tidy->GetCatsProtonLambda(&pL_Usmani, NumMomBins, kMin, kMax, TheSource,
  // 			    TidyCats::pUsmani);
  // pL_Usmani.SetAnaSource(0, 1.2);
  // pL_Usmani.KillTheCat();

  // TGraph* pL_Usmani_1_2 = new TGraph();
  // pL_Usmani_1_2->SetName("pL_Usmani_1_2");
  // for (auto it = 0; it < NumMomBins; ++it) {
  //   double kStar = pL_Usmani.GetMomentum(it);
  //   double mean = pL_Usmani.GetCorrFun(it);
  //   pL_Usmani_1_2->SetPoint(it, kStar, mean);
  // }
  // out->cd();
  // pL_Usmani_1_2->Write("pL_Usmani_1_2");

  
  //TApplication theApp("App",&argc, argv); 
  PlayWithCats *catsPlay = new PlayWithCats();
  // TFile* out_100 = TFile::Open("OutSource_100.root", "RECREATE");
  // TFile* out_200 = TFile::Open("OutSource_200.root", "RECREATE");
  // TFile* out_300 = TFile::Open("OutSource_300.root", "RECREATE");
  // TFile* out_400 = TFile::Open("OutSource_400.root", "RECREATE"); 
  
  // auto gauss_1_2 = new TF1(
  // 			   "gaus12",
  // 			   [&](double *x, double *p) {return 4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) * std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
  // 			   },
  // 			   0, 12, 1);
  // gauss_1_2->SetParameter(0, 1.2);
  // gauss_1_2->SetNpx(1000);
  // gauss_1_2->SetLineColor(kBlue + 7);
  // gauss_1_2->SetLineWidth(3);
  // gauss_1_2->SetLineStyle(1);

  
  // auto gauss_4 =
  //   new TF1(
  // 	    "gaus4",
  // 	    [&](double *x, double *p) {return
  // 				       4. * TMath::Pi() * x[0] * x[0] / std::pow((4. * TMath::Pi() * p[0] * p[0]), 1.5) *
  // 				       std::exp(-x[0] * x[0] / (4. * p[0] * p[0]));
  // 	    },
  // 	    0, 12, 1);
  // gauss_4->SetParameter(0, 4.0);
  // gauss_4->SetNpx(1000);
  // gauss_4->SetLineColor(kPink + 7);
  // gauss_4->SetLineWidth(3);
  // gauss_4->SetLineStyle(1);
  // out->cd();
  // gauss_1_2->Write("gauss_1_2");
  // gauss_4->Write("gauss_4"); 
    
  catsPlay->GenerateSourceDistpxi(out);
  out->Close();
// catsPlay->GenerateSourceDistpp(out_100, 100);
  // catsPlay->GenerateSourceDistpp(out_200, 200);
  // catsPlay->GenerateSourceDistpp(out_300, 300);
  // catsPlay->GenerateSourceDistpp(out_400, 400);
  // catsPlay->GenerateSourceDistpL(out);
  /*
  const double Mass_p = TDatabasePDG::Instance()->GetParticle(2212)->Mass()
      * 1000;
  CATSparameters* cPars = nullptr;
  
  cPars = new CATSparameters(CATSparameters::tSource, 1, true);
  cPars->SetParameter(0, 1.2);
  
  CATS AB_pp_QSOnly; 
  AB_pp_QSOnly.SetAnaSource(GaussSource, *cPars);
  AB_pp_QSOnly.SetUseAnalyticSource(true);
  AB_pp_QSOnly.SetThetaDependentSource(false);
  AB_pp_QSOnly.SetExcludeFailedBins(false);
  AB_pp_QSOnly.SetMomBins(300, 0.5, 300.5);
  // AB_pp_QSOnly.SetNumChannels(2);
  // AB_pp_QSOnly.SetNumPW(0, 0);
  // AB_pp_QSOnly.SetSpin(0, 0);
  // AB_pp_QSOnly.SetChannelWeight(0, 1./4.);
  // AB_pp_QSOnly.SetNumPW(1, 0);
  // AB_pp_QSOnly.SetSpin(1, 1);
  // AB_pp_QSOnly.SetChannelWeight(1, 3./4.);
  AB_pp_QSOnly.SetNumChannels(4);
  AB_pp_QSOnly.SetNumPW(0, 3);  // the maximum number of partial waves flying around = the number of l (s,p,d -> l = 3)
  AB_pp_QSOnly.SetNumPW(1, 3);
  AB_pp_QSOnly.SetNumPW(2, 3);
  AB_pp_QSOnly.SetNumPW(3, 3);
  AB_pp_QSOnly.SetSpin(0, 0);
  AB_pp_QSOnly.SetSpin(1, 1);
  AB_pp_QSOnly.SetSpin(2, 1);
  AB_pp_QSOnly.SetSpin(3, 1);
  AB_pp_QSOnly.SetChannelWeight(0, 3./12.);
  AB_pp_QSOnly.SetChannelWeight(1, 1./12.);
  AB_pp_QSOnly.SetChannelWeight(2, 3./12.);
  AB_pp_QSOnly.SetChannelWeight(3, 5./12.);

  AB_pp_QSOnly.SetQ1Q2(0);
  AB_pp_QSOnly.SetPdgId(2212, 2212);
  AB_pp_QSOnly.SetRedMass((Mass_p * Mass_p) / (Mass_p + Mass_p));
  
  AB_pp_QSOnly.SetAnaSource(0, 1.2);
  AB_pp_QSOnly.KillTheCat();

  TGraph* pp_QS_1_2= new TGraph();
  pp_QS_1_2->SetName("pp_QS_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_QSOnly.GetMomentum(it);
    double mean = AB_pp_QSOnly.GetCorrFun(it);
    pp_QS_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_QS_1_2->Write("pp_QS_1_2");
 
   
  AB_pp_QSOnly.SetAnaSource(0, 4);
  AB_pp_QSOnly.KillTheCat();

  TGraph* pp_QS_4 = new TGraph();
  pp_QS_4->SetName("pp_QS_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_QSOnly.GetMomentum(it);
    double mean = AB_pp_QSOnly.GetCorrFun(it);
    pp_QS_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_QS_4->Write("pp_QS_4");

  AB_pp_QSOnly.SetAnaSource(0, 1.2);
  AB_pp_QSOnly.KillTheCat();
  
  AB_pp_QSOnly.SetChannelWeight(0, 1);
  AB_pp_QSOnly.SetChannelWeight(1, 0);
  AB_pp_QSOnly.SetChannelWeight(2, 0);
  AB_pp_QSOnly.SetChannelWeight(3, 0);
  AB_pp_QSOnly.KillTheCat();
  
  TGraph* pp_QSOnly_1_2_1S0= new TGraph();
  pp_QSOnly_1_2_1S0->SetName("pp_QSOnly_1_2_1S0");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_QSOnly.GetMomentum(it);
    double mean = AB_pp_QSOnly.GetCorrFun(it);
    pp_QSOnly_1_2_1S0->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_QSOnly_1_2_1S0->Write("pp_QSOnly_1_2_1S0");

  AB_pp_QSOnly.SetChannelWeight(0, 0);
  AB_pp_QSOnly.SetChannelWeight(1, 0);
  AB_pp_QSOnly.SetChannelWeight(2, 1);
  AB_pp_QSOnly.SetChannelWeight(3, 0);
  AB_pp_QSOnly.KillTheCat();
  
  TGraph* pp_QSOnly_1_2_3P1= new TGraph();
  pp_QSOnly_1_2_3P1->SetName("pp_QSOnly_1_2_3P1");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_QSOnly.GetMomentum(it);
    double mean = AB_pp_QSOnly.GetCorrFun(it);
    pp_QSOnly_1_2_3P1->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_QSOnly_1_2_3P1->Write("pp_QSOnly_1_2_3P1");

  
  CATS AB_pp_CoulombOnly; 
  AB_pp_CoulombOnly.SetAnaSource(GaussSource, *cPars);
  AB_pp_CoulombOnly.SetUseAnalyticSource(true);
  AB_pp_CoulombOnly.SetThetaDependentSource(false);
  AB_pp_CoulombOnly.SetExcludeFailedBins(false);
  AB_pp_CoulombOnly.SetMomBins(300, 0.5, 300.5);
  AB_pp_CoulombOnly.SetNumChannels(1);
  AB_pp_CoulombOnly.SetNumPW(0, 1);
  AB_pp_CoulombOnly.SetSpin(0, 0);
  AB_pp_CoulombOnly.SetChannelWeight(0, 1.);
  AB_pp_CoulombOnly.SetQ1Q2(1);
  AB_pp_CoulombOnly.SetPdgId(2213, 2212);
  AB_pp_CoulombOnly.SetRedMass((Mass_p * Mass_p) / (Mass_p + Mass_p));
  
  AB_pp_CoulombOnly.SetAnaSource(0, 1.2);
  AB_pp_CoulombOnly.KillTheCat();

  TGraph* pp_Coulomb_1_2= new TGraph();
  pp_Coulomb_1_2->SetName("pp_Coulomb_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_CoulombOnly.GetMomentum(it);
    double mean = AB_pp_CoulombOnly.GetCorrFun(it);
    pp_Coulomb_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Coulomb_1_2->Write("pp_Coulomb_1_2");
 
   
  AB_pp_CoulombOnly.SetAnaSource(0, 4);
  AB_pp_CoulombOnly.KillTheCat();

  TGraph* pp_Coulomb_4 = new TGraph();
  pp_Coulomb_4->SetName("pp_Coulomb_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_CoulombOnly.GetMomentum(it);
    double mean = AB_pp_CoulombOnly.GetCorrFun(it);
    pp_Coulomb_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Coulomb_4->Write("pp_Coulomb_4");
 

  CATS AB_pp_CoulombQS; 
  AB_pp_CoulombQS.SetAnaSource(GaussSource, *cPars);
  AB_pp_CoulombQS.SetUseAnalyticSource(true);
  AB_pp_CoulombQS.SetThetaDependentSource(false);
  AB_pp_CoulombQS.SetExcludeFailedBins(false);
  AB_pp_CoulombQS.SetMomBins(300, 0.5, 300.5);
  AB_pp_CoulombQS.SetNumChannels(2);
  AB_pp_CoulombQS.SetNumPW(0, 0);
  AB_pp_CoulombQS.SetSpin(0, 0);
  AB_pp_CoulombQS.SetChannelWeight(0, 1./4.);
  AB_pp_CoulombQS.SetNumPW(1, 0);
  AB_pp_CoulombQS.SetSpin(1, 1);
  AB_pp_CoulombQS.SetChannelWeight(1, 3./4.);
  AB_pp_CoulombQS.SetQ1Q2(1);
  AB_pp_CoulombQS.SetPdgId(2212, 2212);
  AB_pp_CoulombQS.SetRedMass((Mass_p * Mass_p) / (Mass_p + Mass_p));
  
  AB_pp_CoulombQS.SetAnaSource(0, 1.2);
  AB_pp_CoulombQS.KillTheCat();

  TGraph* pp_Coulomb_QS_1_2= new TGraph();
  pp_Coulomb_QS_1_2->SetName("pp_Coulomb_QS_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_CoulombQS.GetMomentum(it);
    double mean = AB_pp_CoulombQS.GetCorrFun(it);
    pp_Coulomb_QS_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Coulomb_QS_1_2->Write("pp_Coulomb_QS_1_2");
 
   
  AB_pp_CoulombQS.SetAnaSource(0, 4);
  AB_pp_CoulombQS.KillTheCat();

  TGraph* pp_Coulomb_QS_4 = new TGraph();
  pp_Coulomb_QS_4->SetName("pp_Coulomb_QS_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_CoulombQS.GetMomentum(it);
    double mean = AB_pp_CoulombQS.GetCorrFun(it);
    pp_Coulomb_QS_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Coulomb_QS_4->Write("pp_Coulomb_QS_4");
 

  CATS AB_pp_Full; 
  TidyCats* tidy = new TidyCats();
  tidy->GetCatsProtonProton(&AB_pp_Full, 300, 0.5, 300.5, TidyCats::sGaussian);

  AB_pp_Full.SetAnaSource(0, 1.2);
  AB_pp_Full.KillTheCat();

  TGraph* pp_Full_1_2= new TGraph();
  pp_Full_1_2->SetName("pp_Full_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Full.GetMomentum(it);
    double mean = AB_pp_Full.GetCorrFun(it);
    pp_Full_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Full_1_2->Write("pp_Full_1_2");
 
   
  AB_pp_Full.SetAnaSource(0, 4);
  AB_pp_Full.KillTheCat();

  TGraph* pp_Full_4 = new TGraph();
  pp_Full_4->SetName("pp_Full_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Full.GetMomentum(it);
    double mean = AB_pp_Full.GetCorrFun(it);
    pp_Full_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Full_4->Write("pp_Full_4");
  
  CATS AB_pp_Strong; 
  tidy->GetCatsProtonProton(&AB_pp_Strong, 300, 0.5, 300.5, TidyCats::sGaussian);
  AB_pp_Strong.SetQ1Q2(0);
  AB_pp_Strong.SetAnaSource(0, 1.2);
  AB_pp_Strong.SetPdgId(2212, -2212);
  AB_pp_Strong.KillTheCat();

  TGraph* pp_Strong_1_2= new TGraph();
  pp_Strong_1_2->SetName("pp_Strong_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_1_2->Write("pp_Strong_1_2");
 
   
  AB_pp_Strong.SetAnaSource(0, 4);
  AB_pp_Strong.KillTheCat();

  TGraph* pp_Strong_4 = new TGraph();
  pp_Strong_4->SetName("pp_Strong_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_4->Write("pp_Strong_4");

  AB_pp_Strong.SetAnaSource(0,1.2);

  AB_pp_Strong.SetChannelWeight(0, 1);
  AB_pp_Strong.SetChannelWeight(1, 0);
  AB_pp_Strong.SetChannelWeight(2, 0);
  AB_pp_Strong.SetChannelWeight(3, 0);
  AB_pp_Strong.KillTheCat();
  
  TGraph* pp_Strong_1_2_1S0= new TGraph();
  pp_Strong_1_2_1S0->SetName("pp_Strong_1_2_1S0");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_1_2_1S0->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_1_2_1S0->Write("pp_Strong_1_2_1S0");
 
  
  AB_pp_Strong.SetChannelWeight(0, 0);
  AB_pp_Strong.SetChannelWeight(1, 1);
  AB_pp_Strong.SetChannelWeight(2, 0);
  AB_pp_Strong.SetChannelWeight(3, 0);
  AB_pp_Strong.KillTheCat();
  
  TGraph* pp_Strong_1_2_3P0= new TGraph();
  pp_Strong_1_2_3P0->SetName("pp_Strong_1_2_3P0");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_1_2_3P0->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_1_2_3P0->Write("pp_Strong_1_2_3P0");

  
  AB_pp_Strong.SetChannelWeight(0, 0);
  AB_pp_Strong.SetChannelWeight(1, 0);
  AB_pp_Strong.SetChannelWeight(2, 1);
  AB_pp_Strong.SetChannelWeight(3, 0);
  AB_pp_Strong.KillTheCat();
  
  TGraph* pp_Strong_1_2_3P1= new TGraph();
  pp_Strong_1_2_3P1->SetName("pp_Strong_1_2_3P1");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_1_2_3P1->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_1_2_3P1->Write("pp_Strong_1_2_3P1");

  
  AB_pp_Strong.SetChannelWeight(0, 0);
  AB_pp_Strong.SetChannelWeight(1, 0);
  AB_pp_Strong.SetChannelWeight(2, 0);
  AB_pp_Strong.SetChannelWeight(3, 1);
  AB_pp_Strong.KillTheCat();
  
  TGraph* pp_Strong_1_2_3P2= new TGraph();
  pp_Strong_1_2_3P2->SetName("pp_Strong_1_2_3P2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_Strong.GetMomentum(it);
    double mean = AB_pp_Strong.GetCorrFun(it);
    pp_Strong_1_2_3P2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_Strong_1_2_3P2->Write("pp_Strong_1_2_3P2");
   
  CATS AB_pp_StrongQS; 
  tidy->GetCatsProtonProton(&AB_pp_StrongQS, 300, 0.5, 300.5, TidyCats::sGaussian);
  AB_pp_StrongQS.SetQ1Q2(0);
  AB_pp_StrongQS.SetAnaSource(0, 1.2);
  //AB_pp_StrongQS.SetPdgId(2212, -2212);
  AB_pp_StrongQS.KillTheCat();

  TGraph* pp_StrongQS_1_2= new TGraph();
  pp_StrongQS_1_2->SetName("pp_StrongQS_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_1_2->Write("pp_StrongQS_1_2");
 
   
  AB_pp_StrongQS.SetAnaSource(0, 4);
  AB_pp_StrongQS.KillTheCat();

  TGraph* pp_StrongQS_4 = new TGraph();
  pp_StrongQS_4->SetName("pp_StrongQS_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_4->Write("pp_StrongQS_4");

  AB_pp_StrongQS.SetAnaSource(0,1.2);

  AB_pp_StrongQS.SetChannelWeight(0, 1);
  AB_pp_StrongQS.SetChannelWeight(1, 0);
  AB_pp_StrongQS.SetChannelWeight(2, 0);
  AB_pp_StrongQS.SetChannelWeight(3, 0);
  AB_pp_StrongQS.KillTheCat();
  
  TGraph* pp_StrongQS_1_2_1S0= new TGraph();
  pp_StrongQS_1_2_1S0->SetName("pp_StrongQS_1_2_1S0");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_1_2_1S0->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_1_2_1S0->Write("pp_StrongQS_1_2_1S0");
 
  
  AB_pp_StrongQS.SetChannelWeight(0, 0);
  AB_pp_StrongQS.SetChannelWeight(1, 1);
  AB_pp_StrongQS.SetChannelWeight(2, 0);
  AB_pp_StrongQS.SetChannelWeight(3, 0);
  AB_pp_StrongQS.KillTheCat();
  
  TGraph* pp_StrongQS_1_2_3P0= new TGraph();
  pp_StrongQS_1_2_3P0->SetName("pp_StrongQS_1_2_3P0");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_1_2_3P0->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_1_2_3P0->Write("pp_StrongQS_1_2_3P0");

  
  AB_pp_StrongQS.SetChannelWeight(0, 0);
  AB_pp_StrongQS.SetChannelWeight(1, 0);
  AB_pp_StrongQS.SetChannelWeight(2, 1);
  AB_pp_StrongQS.SetChannelWeight(3, 0);
  AB_pp_StrongQS.KillTheCat();
  
  TGraph* pp_StrongQS_1_2_3P1= new TGraph();
  pp_StrongQS_1_2_3P1->SetName("pp_StrongQS_1_2_3P1");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_1_2_3P1->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_1_2_3P1->Write("pp_StrongQS_1_2_3P1");

  
  AB_pp_StrongQS.SetChannelWeight(0, 0);
  AB_pp_StrongQS.SetChannelWeight(1, 0);
  AB_pp_StrongQS.SetChannelWeight(2, 0);
  AB_pp_StrongQS.SetChannelWeight(3, 1);
  AB_pp_StrongQS.KillTheCat();
  
  TGraph* pp_StrongQS_1_2_3P2= new TGraph();
  pp_StrongQS_1_2_3P2->SetName("pp_StrongQS_1_2_3P2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongQS.GetMomentum(it);
    double mean = AB_pp_StrongQS.GetCorrFun(it);
    pp_StrongQS_1_2_3P2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongQS_1_2_3P2->Write("pp_StrongQS_1_2_3P2");
 
  
  CATS AB_pp_StrongCoulomb; 
  tidy->GetCatsProtonProton(&AB_pp_StrongCoulomb, 300, 0.5, 300.5, TidyCats::sGaussian);
  //AB_pp_StrongCoulomb.SetQ1Q2(0);
  AB_pp_StrongCoulomb.SetAnaSource(0, 1.2);
  AB_pp_StrongCoulomb.SetPdgId(2212, -2212);
  AB_pp_StrongCoulomb.KillTheCat();

  TGraph* pp_StrongCoulomb_1_2= new TGraph();
  pp_StrongCoulomb_1_2->SetName("pp_StrongCoulomb_1_2");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongCoulomb.GetMomentum(it);
    double mean = AB_pp_StrongCoulomb.GetCorrFun(it);
    pp_StrongCoulomb_1_2->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongCoulomb_1_2->Write("pp_StrongCoulomb_1_2");
 
   
  AB_pp_StrongCoulomb.SetAnaSource(0, 4);
  AB_pp_StrongCoulomb.KillTheCat();

  TGraph* pp_StrongCoulomb_4 = new TGraph();
  pp_StrongCoulomb_4->SetName("pp_StrongCoulomb_4");
  for (auto it = 0; it < 300; ++it) {
    double kStar = AB_pp_StrongCoulomb.GetMomentum(it);
    double mean = AB_pp_StrongCoulomb.GetCorrFun(it);
    pp_StrongCoulomb_4->SetPoint(it, kStar, mean);
  }
  out->cd();
  pp_StrongCoulomb_4->Write("pp_StrongCoulomb_4");

  int nRadBins = 2400;
  double dRad = 0.005;
  double radMax = dRad * nRadBins;

  TString Name1S0 = Form("1S0");
  TString Name3P0 = Form("3P0");
  TString Name3P1 = Form("3P1");
  TString Name3P2 = Form("3P2");

  TGraph* gr1S0 = new TGraph(nRadBins-1);
  TGraph* gr3P0 = new TGraph(nRadBins-1);
  TGraph* gr3P1 = new TGraph(nRadBins-1);
  TGraph* gr3P2 = new TGraph(nRadBins-1); 

  gr1S0->SetName("1S0");
  gr3P0->SetName("3P0");
  gr3P1->SetName("3P1");
  gr3P2->SetName("3P2");

  for (int iRad = 1; iRad < nRadBins; ++iRad) {
    double PotPars1S0[10] = { iRad * dRad, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
    double PotPars1D2[10] = { iRad * dRad, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 2, 2 };  //the last 3 digits are s,l,j
    double PotPars3P0[10] = { iRad * dRad, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
    double PotPars3P1[10] = { iRad * dRad, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
    double PotPars3P2[10] = { iRad * dRad, 0, NN_AV18, v18_SingleChannelMagic, 1, 1, 1, 1, 1, 2 };  // magically accounts for the coupling to F something
    gr1S0->SetPoint(iRad-1,iRad*dRad,fDlmPot(PotPars1S0));
    gr3P0->SetPoint(iRad-1,iRad*dRad,fDlmPot(PotPars3P0));
    gr3P1->SetPoint(iRad-1,iRad*dRad,fDlmPot(PotPars3P1));
    gr3P2->SetPoint(iRad-1,iRad*dRad,fDlmPot(PotPars3P2));
  }
  out->cd();
  gr1S0->Write("1S0");
  gr3P0->Write("3P0");
  gr3P1->Write("3P1");
  gr3P2->Write("3P2");
    
  */  

  // TFile* out = TFile::Open("Yuki.root", "recreate");
  // out->cd();
  // catsPlay->GenerateYukiCurves_200515(out);
  // out->Close(); 
  //  catsPlay->CloseFile();

  //const char* Data = (argv[1]) ? argv[1] : "";
  //const char* Fit = (argv[2]) ? argv[2] : "";
  //if (Fit != "")catsPlay->ExtractUncertaintyFit(Fit);
  //if (Data != "")catsPlay->ExtractUncertaintyData(Data);
  //catsPlay->GenerateDefault();
//  catsPlay->ShiftBinning();
  //catsPlay->GenerateCoulombOnly();
  //catsPlay->PlotPotentials();
  //catsPlay->PlotPotentialSum();

//   TidyCats* tidy = new TidyCats();
//
//  TFile* out = TFile::Open("out.root","recreate");
//  for (int irad = 0 ; irad < 7; ++irad) {
//
//    CATS AB_ppup;
//    CATS AB_pplow;
//
//    tidy->GetCatsProtonProton(&AB_ppup, 125, 0, 250, TidyCats::sGaussian);
//    tidy->GetCatsProtonProton(&AB_pplow, 125, 0, 250, TidyCats::sGaussian);
//
//
//    TGraphErrors* graph = new TGraphErrors();
//    AB_ppup.SetAnaSource(0,1.48-0.1*irad);
//    AB_ppup.KillTheCat();
//    AB_pplow.SetAnaSource(0,1.52-0.1*irad);
//    AB_pplow.KillTheCat();
//
//    for (int ikStar = 0; ikStar < 125; ++ikStar ) {
//      double mean  = (AB_ppup.GetCorrFun(ikStar)+AB_pplow.GetCorrFun(ikStar))/2.;
//      graph->SetPoint(ikStar,AB_ppup.GetMomentum(ikStar),mean);
//      double Err = TMath::Abs(AB_ppup.GetCorrFun(ikStar)-mean);
//      graph->SetPointError(ikStar,0,Err);
//    }
//    graph->SetName(TString::Format("mTBin_%u",irad).Data());
//    out->cd();
//    graph->Write();
//  }
//  std::cout << "Done looping \n";
  // out_100->Close();
  // out_200->Close();
  // out_300->Close();
  // out_400->Close();
  return 0;
}
