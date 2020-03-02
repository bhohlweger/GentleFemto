#include "CATSInputPhi.h"
#include "ForgivingReader.h"
#include "DecayQA.h"
#include "TObject.h"
#include <iostream>
#include <TDatabasePDG.h>

CATSInputPhi::CATSInputPhi()
    : fCF_pPhi(nullptr),
      fCF_SidebandUp(nullptr),
      fCF_SidebandLow(nullptr),
      fnProtons(0),
      fnAntiProtons(0),
      fnPhi(0),
      fPurityPhi(0),
      fPurityPhipt(0) {
}

CATSInputPhi::~CATSInputPhi() {
  delete fCF_pPhi;
  delete fCF_SidebandUp;
  delete fCF_SidebandLow;
}

void CATSInputPhi::ReadPhiCorrelationFile(const char* path,
                                                const char* trigger,
                                                const char* suffixChar) {
  TString filename = Form("%s/AnalysisResults.root", path);
  fDreamFile = new ReadDreamFile(3, 3);
  fDreamFile->SetAnalysisFile(filename.Data(), "Results", trigger, suffixChar);
  return;
}

void CATSInputPhi::ObtainCFs(int rebin, float normleft, float normright,
                                int rebinSyst, bool isAllCF) {
//normleft & right in MeV!
  normleft /= 1000.;
  normright /= 1000.;

  if (!fDreamFile) {
    std::cerr
        << "ERROR CATSInputPhi: No File was set via ReadPhiCorrelationFile\n";
    return;
  }
  if (fnormalizationLeft == normleft & fnormalizationRight == normright) {
    std::cerr << "ERROR CATSInputPhi: Already existing normalization \n";
    return;
  }

  fCF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", normleft, normright);
  DreamPair* ApAp = new DreamPair("AntiPart", normleft, normright);

  fCF_pPhi = new DreamCF();
  DreamPair* pPhi = new DreamPair("Part", normleft, normright);
  DreamPair* ApPhi = new DreamPair("AntiPart", normleft, normright);

  fCF_SidebandUp = new DreamCF();
  DreamPair* pPhiSBup = new DreamPair("Part_SB_up", normleft, normright);
  DreamPair* ApPhiSBup = new DreamPair("AntiPart_SB_up", normleft, normright);

  fCF_SidebandLow = new DreamCF();
  DreamPair* pPhiSBlow = new DreamPair("Part_SB_low", normleft, normright);
  DreamPair* ApPhiSBlow = new DreamPair("AntiPart_SB_low", normleft, normright);

  pp->SetPair(fDreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(fDreamFile->GetPairDistributions(1, 1, ""));

  pPhi->SetPair(fDreamFile->GetPairDistributions(0, 2, ""));
  ApPhi->SetPair(fDreamFile->GetPairDistributions(1, 2, ""));

  if (isAllCF) {
    pPhiSBup->SetPair(fDreamFile->GetPairDistributions(0, 2, ""));
    ApPhiSBup->SetPair(fDreamFile->GetPairDistributions(1, 2, ""));

    pPhiSBlow->SetPair(fDreamFile->GetPairDistributions(0, 2, ""));
    ApPhiSBlow->SetPair(fDreamFile->GetPairDistributions(1, 2, ""));
  }

  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pPhi->ShiftForEmpty(pPhi->GetPair());
  ApPhi->ShiftForEmpty(ApPhi->GetPair());

  if (isAllCF) {
    pPhiSBup->ShiftForEmpty(pPhiSBup->GetPair());
    ApPhiSBup->ShiftForEmpty(ApPhiSBup->GetPair());

    pPhiSBlow->ShiftForEmpty(pPhiSBlow->GetPair());
    ApPhiSBlow->ShiftForEmpty(ApPhiSBlow->GetPair());
  }

  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());

  pPhi->FixShift(pPhi->GetPairShiftedEmpty(0),
                   ApPhi->GetPairShiftedEmpty(0), ApPhi->GetFirstBin());
  ApPhi->FixShift(ApPhi->GetPairShiftedEmpty(0),
                     pPhi->GetPairShiftedEmpty(0), pPhi->GetFirstBin());

  if (isAllCF) {
    pPhiSBup->FixShift(pPhiSBup->GetPairShiftedEmpty(0),
                      pPhi->GetPairShiftedEmpty(0),
                      ApPhi->GetPairShiftedEmpty(0), pPhi->GetFirstBin(),
                      ApPhi->GetFirstBin());
    ApPhiSBup->FixShift(ApPhiSBup->GetPairShiftedEmpty(0),
                        pPhi->GetPairShiftedEmpty(0),
                        ApPhi->GetPairShiftedEmpty(0), pPhi->GetFirstBin(),
                        ApPhi->GetFirstBin());

    pPhiSBlow->FixShift(pPhiSBlow->GetPairShiftedEmpty(0),
                       pPhi->GetPairShiftedEmpty(0),
                       ApPhi->GetPairShiftedEmpty(0), pPhi->GetFirstBin(),
                       ApPhi->GetFirstBin());
    ApPhiSBlow->FixShift(ApPhiSBlow->GetPairShiftedEmpty(0),
                         pPhi->GetPairShiftedEmpty(0),
                         ApPhi->GetPairShiftedEmpty(0),
                         pPhi->GetFirstBin(), ApPhi->GetFirstBin());
  }

  if (rebinSyst != 1) {
    pp->Rebin(pp->GetPairFixShifted(0), rebinSyst, true);
    ApAp->Rebin(ApAp->GetPairFixShifted(0), rebinSyst, true);
  }

  pPhi->Rebin(pPhi->GetPairFixShifted(0), rebin * rebinSyst, true);
  ApPhi->Rebin(ApPhi->GetPairFixShifted(0), rebin * rebinSyst, true);

  if (isAllCF) {
    pPhiSBup->Rebin(pPhiSBup->GetPairFixShifted(0), rebin * rebinSyst, true);
    ApPhiSBup->Rebin(ApPhiSBup->GetPairFixShifted(0), rebin * rebinSyst, true);

    pPhiSBlow->Rebin(pPhiSBlow->GetPairFixShifted(0), rebin * rebinSyst, true);
    ApPhiSBlow->Rebin(ApPhiSBlow->GetPairFixShifted(0), rebin * rebinSyst, true);
  }

  if (rebinSyst != 1) {
    pp->ReweightMixedEvent(pp->GetPairRebinned(0), 0.2, 0.9, pp->GetPair());
    ApAp->ReweightMixedEvent(ApAp->GetPairRebinned(0), 0.2, 0.9, ApAp->GetPair());
  } else {
    pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9, pp->GetPair());
    ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9, ApAp->GetPair());
  }

  pPhi->ReweightMixedEvent(pPhi->GetPairRebinned(0), 0.2, 0.9, pPhi->GetPair());
  ApPhi->ReweightMixedEvent(ApPhi->GetPairRebinned(0), 0.2, 0.9, ApPhi->GetPair());

  if (isAllCF) {
    pPhiSBup->ReweightMixedEvent(pPhiSBup->GetPairRebinned(0), 0.2, 0.9, pPhiSBup->GetPair());
    ApPhiSBup->ReweightMixedEvent(ApPhiSBup->GetPairRebinned(0), 0.2, 0.9, ApPhiSBup->GetPair());

    pPhiSBlow->ReweightMixedEvent(pPhiSBlow->GetPairRebinned(0), 0.2, 0.9, pPhiSBlow->GetPair());
    ApPhiSBlow->ReweightMixedEvent(ApPhiSBlow->GetPairRebinned(0), 0.2, 0.9, ApPhiSBlow->GetPair());
  }

  fCF_pp->SetPairs(pp, ApAp);
  fCF_pp->GetCorrelations("pp");

  fCF_pPhi->SetPairs(pPhi, ApPhi);
  fCF_pPhi->GetCorrelations("pPhi");

  if (isAllCF) {
    fCF_SidebandUp->SetPairs(pPhiSBup, ApPhiSBup);
    fCF_SidebandUp->GetCorrelations("pPhiSBUp");

    fCF_SidebandLow->SetPairs(pPhiSBlow, ApPhiSBlow);
    fCF_SidebandLow->GetCorrelations("pPhiSBLow");
  }
  fnormalizationLeft = normleft;
  fnormalizationRight = normright;
}

void CATSInputPhi::CountPairs(const char* path, const char* trigger,
                                 const char* suffixChar) {
  TString filename = Form("%s/AnalysisResults.root", path);
//  ForgivingReader* reader = new ForgivingReader(filename, trigger, suffixChar);

//  auto pTProton = (TH1F*) reader->Get1DHistInList(reader->GetTrackCuts(), {
//                                                      "pTDist_after" });
//  if (pTProton) {
//    fnProtons = pTProton->GetEntries();
//    delete pTProton;
//  } else {
//    pTProton = (TH1F*) reader->Get1DHistInList(
//        reader->GetListInList(reader->GetTrackCuts(), { { "after" } }), {
//            "pTDist_after" });
//    fnProtons = pTProton->GetEntries();
//    delete pTProton;
//  }

//  auto pTAntiProton = (TH1F*) reader->Get1DHistInList(
//      reader->GetAntiTrackCuts(), { "pTDist_after" });
//  if (pTAntiProton) {
//    fnAntiProtons = pTAntiProton->GetEntries();
//    delete pTAntiProton;
//  } else {
//    pTAntiProton = (TH1F*) reader->Get1DHistInList(
//        reader->GetListInList(reader->GetAntiTrackCuts(), { { "after" } }), {
//            "pTDist_after" });
//    fnAntiProtons = pTAntiProton->GetEntries();
//    delete pTAntiProton;
//  }

  auto file = TFile::Open(filename);
  TList *dirResults;
  TDirectoryFile *dir=(TDirectoryFile*)(file->FindObjectAny(Form("%sResults%s", trigger, suffixChar)));
  dir->GetObject(Form("%sResults%s", trigger, suffixChar), dirResults);

  //PROTON
  TList* plist = (TList*)dirResults->FindObject("Proton");
  if (suffixChar == "0")
  {
    TList* tracklistp =  (TList*)plist->FindObject("after");
    auto pTProton1 = (TH1F*)tracklistp->FindObject("pTDist_after");
    if (pTProton1) {
     fnProtons = pTProton1->GetEntries();
        delete pTProton1;
    } else {
        std::cout<<"no pT proton found"<<std::endl;
        delete pTProton1;
    }
  }
  else {
      auto pTProton2 = (TH1F*)plist->FindObject("pTDist_after");
      if (pTProton2) {
       fnProtons = pTProton2->GetEntries();
          delete pTProton2;
      } else {
          std::cout<<"no pT proton found"<<std::endl;
          delete pTProton2;
      }
  }

  //ANTIPROTON
  TList* aplist = (TList*)dirResults->FindObject("AntiProton");
  if (suffixChar == "0")
  {
    TList* tracklistap =  (TList*)aplist->FindObject("after");
    auto pTAntiProton1 = (TH1F*)tracklistap->FindObject("pTDist_after");
    if (pTAntiProton1) {
     fnAntiProtons = pTAntiProton1->GetEntries();
        delete pTAntiProton1;
    } else {
        std::cout<<"no pT anti proton found"<<std::endl;
        delete pTAntiProton1;
    }
  }
  else {
      auto pTAntiProton2 = (TH1F*)aplist->FindObject("pTDist_after");
      if (pTAntiProton2) {
       fnAntiProtons = pTAntiProton2->GetEntries();
          delete pTAntiProton2;
      } else {
          std::cout<<"no pT ani proton found"<<std::endl;
          delete pTAntiProton2;
      }
  }

  //PHI
//  TList* philist = (TList*)dirResults->FindObject("Phi");
//  DecayQA* PhiQA = new DecayQA("#Phi", "#Lambda#gamma");
//  ForgivingFitter* fFitter= new ForgivingFitter();

//  if (suffixChar == "0")
//  {
//    TList* listphi =  (TList*)philist->FindObject("v0Cuts");
//    auto pTPhi1 = (TH1F*)listphi->FindObject("InvMassPt");
//    PhiQA->SetRangesFitting(0.99, 1.06, 0.99, 1.06);
//    PhiQA->SetInvMasspTStartBin(1);
//    PhiQA->SetIMHistoScale(1.75, 0.8, 0.35);
//    fFitter->FitInvariantMassPhi(pTPhi1, 0.008);
//    fPurityPhi = PhiQA->GetPurity();
//    fnPhi = PhiQA->GetSignalCounts() + PhiQA->GetBackgroundCounts();

//    PhiQA->FitInvariantMassPhi(0.003, "Phi", "");
//    auto purityGraph = PhiQA->GetPurityGraph();
//    fPurityPhipt = PhiQA->GetPurity(averagePtLowKstarPhi);
//    Warning("CATSInputPhi", "pT averaged purity: %.2f, purity at the pT of the Phi at low k*: %.2f", fPurityPhi, fPurityPhipt);

//    delete PhiQA;
//    delete pTPhi1;
//  }
//  else {
//      TList* tracklistphi2 =  (TList*)philist->FindObject("MinimalBooking");
//      auto pTPhi2 = (TH1F*)tracklistphi2->FindObject("InvMassPt");
//      PhiQA->SetRangesFitting(0.99, 1.06, 0.99, 1.06);
//      PhiQA->SetInvMasspTStartBin(1);
//      PhiQA->SetIMHistoScale(1.75, 0.8, 0.35);
//      fFitter->FitInvariantMassPhi(pTPhi2, 0.008, "");
//      fPurityPhi = PhiQA->GetPurity();
//      fnPhi = PhiQA->GetSignalCounts() + PhiQA->GetBackgroundCounts();

//      PhiQA->FitInvariantMassPhi(0.003, "Phi", "");
//      auto purityGraph = PhiQA->GetPurityGraph();
//      fPurityPhipt = PhiQA->GetPurity(averagePtLowKstarPhi);
//      Warning("CATSInputPhi", "pT averaged purity: %.2f, purity at the pT of the Phi at low k*: %.2f", fPurityPhi, fPurityPhipt);

//      delete PhiQA;
//      delete pTPhi2;
//  }






//    TList* philist = (TList*)dirResults->FindObject("Phi");
//    if (suffixChar == "0")
//    {
//      TList* listphi =  (TList*)philist->FindObject("v0Cuts");
//      auto pTPhi1 = (TH1F*)listphi->FindObject("InvMassPt");

////      fPurityPhi = PhiQA->GetPurity();
////      fnPhi = PhiQA->GetSignalCounts() + PhiQA->GetBackgroundCounts();

//      pTPhi1->SetTitle("; M_{K^{-}K^{+}} (GeV/#it{c^{2}}); Entries");
//      pTPhi1->GetXaxis()->SetRangeUser(0.99, 1.06);

//    int nPars = 7; // number of parameters for the fit

//    auto test = new TF1("peak", [&](double *x, double *p) {

//             return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3])+(p[4]+p[5]*x[0]+p[6]*TMath::Power(x[0],2));

//            }, 0.99, 2, nPars);

//    test->SetParameter(0,pTPhi1->GetBinContent(pTPhi1->GetMaximumBin())*1.2);
//    test->SetParameter(1,1.02);//peak phi
//    test->SetParameter(2,0.003);//sigma gauss
//    test->FixParameter(3,0.00425);//gamma decay width
//    test->SetParLimits(1,1.0,1.04);
//    test->SetParLimits(2,0.0005,0.002);
//    //test->SetParameter(4,4);//r
//    //test->FixParameter(0);
//    pTPhi1->Fit(test);


//    int mPars = 3;
//    auto background = new TF1("background", [&](double *x, double *p) {

//             return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2);

//            }, 0.99, 2, mPars);

//    background->FixParameter(0,test->GetParameter(4));
//    background->FixParameter(1,test->GetParameter(5));
//    background->FixParameter(2,test->GetParameter(6));
//    background->SetLineColor(kRed+2);


//    int oPars = 4; // number of parameters for the fit

//    auto signal = new TF1("signal", [&](double *x, double *p) {

//             return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3]);

//            }, 0.99, 2, oPars);

//    signal->FixParameter(0,test->GetParameter(0));
//    signal->FixParameter(1,test->GetParameter(1));
//    signal->FixParameter(2,test->GetParameter(2));
//    signal->FixParameter(3,test->GetParameter(3));
//    signal->SetLineColor(kBlue+2);

//    double massPhi=TDatabasePDG::Instance()->GetParticle(333)->Mass();
//    double breite=0.008;
//    double intmin= massPhi-breite;
//    double intmax= massPhi+breite;
//    double a=test->Integral(intmin,intmax);
//    double b=background->Integral(intmin,intmax);
//    double s=signal->Integral(intmin,intmax);
//    fPurityPhi=(s/a);
//   // std::cout<<"purity: "<<fPurityPhi<<std::endl;
//    fnPhi=a;

//    }
//    else {
//        TList* tracklistphi2 =  (TList*)philist->FindObject("MinimalBooking");
//        auto pTPhi2 = (TH1F*)tracklistphi2->FindObject("InvMassPt");

//        pTPhi2->SetTitle("; M_{K^{-}K^{+}} (GeV/#it{c^{2}}); Entries");
//        pTPhi2->GetXaxis()->SetRangeUser(0.99, 1.06);

//      int nPars = 7; // number of parameters for the fit

//      auto test = new TF1("peak", [&](double *x, double *p) {

//               return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3])+(p[4]+p[5]*x[0]+p[6]*TMath::Power(x[0],2));

//              }, 0.99, 2, nPars);

//      test->SetParameter(0,pTPhi2->GetBinContent(pTPhi2->GetMaximumBin())*1.2);
//      test->SetParameter(1,1.02);//peak phi
//      test->SetParameter(2,0.003);//sigma gauss
//      test->FixParameter(3,0.00425);//gamma decay width
//      test->SetParLimits(1,1.0,1.04);
//      test->SetParLimits(2,0.0005,0.002);
//      //test->SetParameter(4,4);//r
//      //test->FixParameter(0);
//      pTPhi2->Fit(test);


//      int mPars = 3;
//      auto background = new TF1("background", [&](double *x, double *p) {

//               return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2);

//              }, 0.99, 2, mPars);

//      background->FixParameter(0,test->GetParameter(4));
//      background->FixParameter(1,test->GetParameter(5));
//      background->FixParameter(2,test->GetParameter(6));
//      background->SetLineColor(kRed+2);


//      int oPars = 4; // number of parameters for the fit

//      auto signal = new TF1("signal", [&](double *x, double *p) {

//               return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3]);

//              }, 0.99, 2, oPars);

//      signal->FixParameter(0,test->GetParameter(0));
//      signal->FixParameter(1,test->GetParameter(1));
//      signal->FixParameter(2,test->GetParameter(2));
//      signal->FixParameter(3,test->GetParameter(3));
//      signal->SetLineColor(kBlue+2);

//      double massPhi=TDatabasePDG::Instance()->GetParticle(333)->Mass();
//      double breite=0.008;
//      double intmin= massPhi-breite;
//      double intmax= massPhi+breite;
//      double a=test->Integral(intmin,intmax);
//      double b=background->Integral(intmin,intmax);
//      double s=signal->Integral(intmin,intmax);
//      fPurityPhi=(s/a);
//      //fPurityPhipt=(s/a);
//       std::cout<<"purity: "<<fPurityPhi<<std::endl;

//      fnPhi=a;
//    }




  TList* philist = (TList*)dirResults->FindObject("Phi");
  TH1F* pTPhi= nullptr;


  if (suffixChar == "0") {
   TList* listphi =  (TList*)philist->FindObject("v0Cuts");
   pTPhi =  (TH1F*)(((TH2F*)listphi->FindObject("InvMassPt"))->ProjectionY());
   if(pTPhi==nullptr){
        std::cout<<"no pTPhi!!"<<std::endl;
   }
  }
  else {
      TList* tracklistphi2 =  (TList*)philist->FindObject("MinimalBooking");
      pTPhi = (TH1F*)(((TH2F*)tracklistphi2->FindObject("InvMassPt"))->ProjectionY());
      if(pTPhi==nullptr){
           std::cout<<"no pTPhi!!"<<std::endl;
      }
  }

//      fPurityPhi = PhiQA->GetPurity();
//      fnPhi = PhiQA->GetSignalCounts() + PhiQA->GetBackgroundCounts();

    pTPhi->SetTitle("; M_{K^{-}K^{+}} (GeV/#it{c^{2}}); Entries");
    pTPhi->GetXaxis()->SetRangeUser(0.99, 1.06);

  int nPairs = 7; // number of parameters for the fit

  auto test = new TF1("peak", [&](double *x, double *p) {

           return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3])+(p[4]+p[5]*x[0]+p[6]*TMath::Power(x[0],2));

          }, 0.99, 2, nPairs);

  test->SetParameter(0,pTPhi->GetBinContent(pTPhi->GetMaximumBin())*1.2);
  test->SetParameter(1,1.02);//peak phi
  test->SetParameter(2,0.003);//sigma gauss
  test->FixParameter(3,0.00425);//gamma decay width
  test->SetParLimits(1,1.0,1.04);
  test->SetParLimits(2,0.0005,0.002);
  //test->SetParameter(4,4);//r
  //test->FixParameter(0);
  pTPhi->Fit(test);


  int mPairs = 3;
  auto background = new TF1("background", [&](double *x, double *p) {

           return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2);

          }, 0.99, 2, mPairs);

  background->FixParameter(0,test->GetParameter(4));
  background->FixParameter(1,test->GetParameter(5));
  background->FixParameter(2,test->GetParameter(6));
  background->SetLineColor(kRed+2);


  int oPars = 4; // number of parameters for the fit

  auto signal = new TF1("signal", [&](double *x, double *p) {

           return p[0]*TMath::Voigt(x[0]-p[1],p[2],p[3]);

          }, 0.99, 2, oPars);

  signal->FixParameter(0,test->GetParameter(0));
  signal->FixParameter(1,test->GetParameter(1));
  signal->FixParameter(2,test->GetParameter(2));
  signal->FixParameter(3,test->GetParameter(3));
  signal->SetLineColor(kBlue+2);

  double massPhi=TDatabasePDG::Instance()->GetParticle(333)->Mass();
  double breite=0.008;
  double intmin= massPhi-breite;
  double intmax= massPhi+breite;
  double a=test->Integral(intmin,intmax);
  double b=background->Integral(intmin,intmax);
  double s=signal->Integral(intmin,intmax);
  fPurityPhi=(s/a);
  std::cout<<"inmin: "<<intmin<<std::endl;
  std::cout<<"intmax: "<<intmax<<std::endl;
  std::cout<<"all: "<<a<<std::endl;
  std::cout<<"bg: "<<b<<std::endl;
  std::cout<<"signal: "<<s<<std::endl;
  std::cout<<"purity: "<<fPurityPhi<<std::endl;
  fnPhi=a;
  delete pTPhi;





//  TList* philist = (TList*)dirResults->FindObject("Phi");
//  if (suffixChar == "0")
//  {
//    TList* listphi =  (TList*)philist->FindObject("v0Cuts");
//    TList* tracklistphi1 =  (TList*)listphi->FindObject("after");
//    auto pTPhi1 = (TH1F*)tracklistphi1->FindObject("pTDist_after");
//    if (pTPhi1) {
//     fnPhi = pTPhi1->GetEntries();
//        delete pTPhi1;
//    } else {
//        std::cout<<"no pT phi found"<<std::endl;
//        delete pTPhi1;
//    }
//  }
//  else {
//      TList* tracklistphi2 =  (TList*)philist->FindObject("PosCuts");
//      auto pTPhi2 = (TH1F*)tracklistphi2->FindObject("pTDist_after");
//      if (pTPhi2) {
//       fnPhi = pTPhi2->GetEntries();
//          delete pTPhi2;
//      } else {
//          std::cout<<"no pT phi found"<<std::endl;
//          delete pTPhi2;
//      }
//  }

//----------------------------------------------------------------------------------

//  float averagePtLowKstarPhi = 0;
//  auto ptList = reader->GetOtherCuts("PairQA");
//  auto pPhiList = (TList*)ptList->FindObject("QA_Particle0_Particle2");
//  auto pbarPhiList = (TList*)ptList->FindObject("QA_Particle1_Particle2");
//  if (pPhiList && pbarPhiList) {
//    auto pPhiptHist = (TH2F*) pPhiList->FindObject("PtQA_Particle0_Particle2");
//    pPhiptHist->Add(
//        (TH2F*) pbarPhiList->FindObject("PtQA_Particle1_Particle2"));
//    averagePtLowKstarPhi = pPhiptHist->GetMean(2);
//    Warning("CATSInputPhi",
//            "Average pT at low kstar - p: %.2f - Phi : %.2f",
//            pPhiptHist->GetMean(1), pPhiptHist->GetMean(2));
//  }

//  DecayQA* PhiQA = new DecayQA("#Sigma", "#Lambda#gamma");
//  PhiQA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
//  PhiQA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
//  PhiQA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
//  PhiQA->SetInvMasspTStartBin(3);
//  PhiQA->SetIMHistoScale(1.75, 0.8, 0.35);
//  PhiQA->GetPeriodQAPhi(0.003, "Sigma0");
//  fPuritySigma0 = PhiQA->GetPurity();
//  fnPhi = PhiQA->GetSignalCounts() + PhiQA->GetBackgroundCounts();

//  PhiQA->InvariantMassPhi(0.003, "Sigma0", false);
//  auto purityGraph = PhiQA->GetPurityGraph();
//  fPurityPhipt = PhiQA->GetPurity(averagePtLowKstarSigma0);
//  Warning("CATSInputSigma0", "pT averaged purity: %.2f, purity at the pT of the Sigma0 at low k*: %.2f", fPuritySigma0, fPuritySigma0pt);

//  delete PhiQA;
//  delete reader;




}

TH1F* CATSInputPhi::GetCF(TString pair, TString hist) {
  TH1F* output = nullptr;
  if (pair == TString("pp")) {
    for (const auto &it : fCF_pp->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhi")) {
    for (const auto &it : fCF_pPhi->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhiSBUp")) {
    for (const auto &it : fCF_SidebandUp->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhiSBLow")) {
    for (const auto &it : fCF_SidebandLow->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else {
    std::cout << pair << " does not exist\n";
  }
  if (!output) {
    std::cout << "Danger! Histogram not set, maybe histname " << hist
              << " does not exist? \n";
  }
  return (TH1F*) output->Clone(Form("%sCloned", output->GetName()));
}

TGraphAsymmErrors* CATSInputPhi::GetCFGr(TString pair, TString hist) {
  TGraphAsymmErrors* output = nullptr;
  if (pair == TString("pp")) {
    for (const auto &it : fCF_pp->GetCorrelationFunctionGraphs()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhi")) {
    for (const auto &it : fCF_pPhi->GetCorrelationFunctionGraphs()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhiSBUp")) {
    for (const auto &it : fCF_SidebandUp->GetCorrelationFunctionGraphs()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pPhiSBLow")) {
    for (const auto &it : fCF_SidebandLow->GetCorrelationFunctionGraphs()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else {
    std::cout << pair << " does not exist\n";
  }
  if (!output) {
    std::cout << "Danger! Histogram not set, maybe histname " << hist
              << " does not exist? \n";
  }
  return (TGraphAsymmErrors*) output->Clone(Form("%sCloned", output->GetName()));
}

unsigned int CATSInputPhi::GetFemtoPairs(float kMin, float kMax,
                                            TString pair) {
  if (pair == TString("pp")) {
    return fCF_pp->GetFemtoPairs(kMin, kMax);
  } else if (pair == TString("pPhi")) {
    return fCF_pPhi->GetFemtoPairs(kMin, kMax);
  } else if (pair == TString("pPhiSBUp")) {
    return fCF_SidebandUp->GetFemtoPairs(kMin, kMax);
  } else if (pair == TString("pPhiSBLow")) {
    return fCF_SidebandLow->GetFemtoPairs(kMin, kMax);
  } else {
    std::cout << pair << " does not exist\n";
  }
  return 0;
}
