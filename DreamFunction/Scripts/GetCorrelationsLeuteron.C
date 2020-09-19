/*
 * GetCorrelationsLeuteron.C
 *
 *  Created on: 27 January 2020
 *      Author: Michael Jung
 */

#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "ReadDreamFile.h"
#include <iostream>

void GetCorrelationsLeuteron(const char* file, const char* prefix = "", const char* addon = "", bool BruteForceDebugging = true){

  ReadDreamFile* DreamFile = new ReadDreamFile(6,6);	  // number of particle species (Proton, Antiproton, Deuteron, Antideuteron, Lambda, Antilambda)

    // Numbering scheme of the particles:
    // 0: Proton
    // 1: Antiproton
    // 2: Deuteron
    // 3: Antideuteron
    // 4: Lambda
    // 5: Antilambda

  DreamFile->SetAnalysisFile(file,prefix,addon);
    // SetAnalysisFile(1,2,3)
    // 1. argument (char) name ond path of the root file
    // 2. argument (char) if the final results folder is called "ppResults0", "pp" is the prefix 
    // 3. argument (char) if the final results folder is called "ppResults0", "0" is the addon

  DreamCF* CF_ProtonProton = new DreamCF();
  DreamCF* CF_ProtonAntiproton = new DreamCF();
  DreamCF* CF_ProtonDeuteron = new DreamCF();
  DreamCF* CF_ProtonAntideuteron = new DreamCF();
  DreamCF* CF_ProtonLambda = new DreamCF();
  DreamCF* CF_ProtonAntilambda = new DreamCF();

  DreamCF* CF_DeuteronDeuteron = new DreamCF();
  DreamCF* CF_DeuteronAntideuteron = new DreamCF();
  DreamCF* CF_DeuteronLambda = new DreamCF();
  DreamCF* CF_DeuteronAntilambda = new DreamCF();

  DreamCF* CF_LambdaLambda = new DreamCF();
  DreamCF* CF_LambdaAntilambda = new DreamCF();


  DreamPair* ProtonProton = new DreamPair("ProtonProton",0.24,0.34);
  DreamPair* ProtonAntiproton = new DreamPair("ProtonAntiproton",0.24,0.34);
  DreamPair* ProtonDeuteron = new DreamPair("ProtonDeuteron",0.24,0.34);
  DreamPair* ProtonAntideuteron = new DreamPair("ProtonAntideuteron",0.24,0.34);
  DreamPair* ProtonLambda = new DreamPair("ProtonLambda",0.24,0.34);
  DreamPair* ProtonAntilambda = new DreamPair("ProtonAntilambda",0.24,0.34);

  DreamPair* AntiprotonAntiproton = new DreamPair("AntiprotonAntiproton",0.24,0.34);
  DreamPair* AntiprotonDeuteron = new DreamPair("AntiprotonDeuteron",0.24,0.34);
  DreamPair* AntiprotonAntideuteron = new DreamPair("AntiprotonAntideuteron",0.24,0.34);
  DreamPair* AntiprotonLambda = new DreamPair("AntiprotonLambda",0.24,0.34);
  DreamPair* AntiprotonAntilambda = new DreamPair("AntiprotonAntilambda",0.24,0.34);

  DreamPair* DeuteronDeuteron = new DreamPair("DeuteronDeuteron",0.24,0.34);
  DreamPair* DeuteronAntideuteron = new DreamPair("DeuteronAntideuteron",0.24,0.34);
  DreamPair* DeuteronLambda = new DreamPair("DeuteronLambda",0.24,0.34);
  DreamPair* DeuteronAntilambda = new DreamPair("DeuteronAntilambda",0.24,0.34);

  DreamPair* AntideuteronAntideuteron = new DreamPair("AntideuteronAntideuteron",0.24,0.34);
  DreamPair* AntideuteronLambda = new DreamPair("AntideuteronLambda",0.24,0.34);
  DreamPair* AntideuteronAntilambda = new DreamPair("AntideuteronAntilambda",0.24,0.34);

  DreamPair* LambdaLambda = new DreamPair("LambdaLambda",0.24,0.34);
  DreamPair* LambdaAntilambda = new DreamPair("LambdaAntilambda",0.24,0.34);

  DreamPair* AntilambdaAntilambda = new DreamPair("AntilambdaAntilambda",0.24,0.34);

    // Dreampair(1,2,3)
    // 1. argument (char) name of the object
    // 2. argument (float) lower limit for the determination of the normalization parameter (N) for mixed and same event yields
    // 3. argument (float) upper limit for the determination of the normalization parameter (N) for mixed and same event yields
    // -> correlation function should be flat within these two limits!


  std::cout << ""			   << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << ""			   << std::endl;

    // GetPairDistributions(1,2,3)
    // 1. argument (integer) particle 1
    // 2. argument (integer) particle 2
    // 3. argument (char) name
  
    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonProton pairs" << std::endl;}
  ProtonProton->SetPair(DreamFile->GetPairDistributions(0,0,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonProton pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonAntiproton pairs" << std::endl;}
  ProtonAntiproton->SetPair(DreamFile->GetPairDistributions(0,1,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonAntiproton pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonDeuteron pairs" << std::endl;}
  ProtonDeuteron->SetPair(DreamFile->GetPairDistributions(0,2,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonDeuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonAntideuteron pairs" << std::endl;}
  ProtonAntideuteron->SetPair(DreamFile->GetPairDistributions(0,3,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonAntideuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonLambda pairs" << std::endl;}
  ProtonLambda->SetPair(DreamFile->GetPairDistributions(0,4,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonLambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting ProtonAntilambda pairs" << std::endl;}
  ProtonAntilambda->SetPair(DreamFile->GetPairDistributions(0,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonAntilambda pairs recieved" << std::endl;}


    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntiprotonAntiproton pairs" << std::endl;}
  AntiprotonAntiproton->SetPair(DreamFile->GetPairDistributions(1,1,""));
    if(BruteForceDebugging){std::cout << "x-x-> ProtonAntilambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntiprotonDeuteron pairs" << std::endl;}
  AntiprotonDeuteron->SetPair(DreamFile->GetPairDistributions(1,2,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntiprotonDeuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntiprotonAntideuteron pairs" << std::endl;}
  AntiprotonAntideuteron->SetPair(DreamFile->GetPairDistributions(1,3,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntiprotonAntideuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntiprotonLambda pairs" << std::endl;}
  AntiprotonLambda->SetPair(DreamFile->GetPairDistributions(1,4,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntiprotonLambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntiprotonAntilambda pairs" << std::endl;}
  AntiprotonAntilambda->SetPair(DreamFile->GetPairDistributions(1,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntiprotonAntilambda pairs recieved" << std::endl;}


    if(BruteForceDebugging){std::cout << "x-x-> Start getting DeuteronDeuteron pairs" << std::endl;}
  DeuteronDeuteron->SetPair(DreamFile->GetPairDistributions(2,2,""));
    if(BruteForceDebugging){std::cout << "x-x-> DeuteronDeuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting DeuteronAntideuteron pairs" << std::endl;}
  DeuteronAntideuteron->SetPair(DreamFile->GetPairDistributions(2,3,""));
    if(BruteForceDebugging){std::cout << "x-x-> DeuteronAntideuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting DeuteronLambda pairs" << std::endl;}
  DeuteronLambda->SetPair(DreamFile->GetPairDistributions(2,4,""));
    if(BruteForceDebugging){std::cout << "x-x-> DeuteronLambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting DeuteronAntilambda pairs" << std::endl;}
  DeuteronAntilambda->SetPair(DreamFile->GetPairDistributions(2,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntideuteronAntilambda pairs recieved" << std::endl;}


    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntideuteronAntideuteron pairs" << std::endl;}
  AntideuteronAntideuteron->SetPair(DreamFile->GetPairDistributions(3,3,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntideuteronAntideuteron pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntideuteronLambda pairs" << std::endl;}
  AntideuteronLambda->SetPair(DreamFile->GetPairDistributions(3,4,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntideuteronLambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntideuteronAntilambda pairs" << std::endl;}
  AntideuteronAntilambda->SetPair(DreamFile->GetPairDistributions(3,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntideuteronAntilambda pairs recieved" << std::endl;}


    if(BruteForceDebugging){std::cout << "x-x-> Start getting LambdaLambda pairs" << std::endl;}
  LambdaLambda->SetPair(DreamFile->GetPairDistributions(4,4,""));
    if(BruteForceDebugging){std::cout << "x-x-> LambdaLambda pairs recieved" << std::endl;}

    if(BruteForceDebugging){std::cout << "x-x-> Start getting LambdaAntilambda pairs" << std::endl;}
  LambdaAntilambda->SetPair(DreamFile->GetPairDistributions(4,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> LambdaAntilambda pairs recieved" << std::endl;}


    if(BruteForceDebugging){std::cout << "x-x-> Start getting AntilambdaAntilambda pairs" << std::endl;}
  AntilambdaAntilambda->SetPair(DreamFile->GetPairDistributions(5,5,""));
    if(BruteForceDebugging){std::cout << "x-x-> AntilambdaAntilambda pairs recieved" << std::endl;}
 

  std::cout << ""			   << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << ""			   << std::endl;
  
  DeuteronLambda->ShiftForEmpty(DeuteronLambda->GetPair());
  AntideuteronAntilambda->ShiftForEmpty(AntideuteronAntilambda->GetPair());

  std::cout << ""			   << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << ""			   << std::endl;

  DeuteronLambda->FixShift(DeuteronLambda->GetPairShiftedEmpty(0), AntideuteronAntilambda->GetPairShiftedEmpty(0),AntideuteronAntilambda->GetFirstBin());
  AntideuteronAntilambda->FixShift(AntideuteronAntilambda->GetPairShiftedEmpty(0), DeuteronLambda->GetPairShiftedEmpty(0),DeuteronLambda->GetFirstBin());

  DeuteronLambda->FixShift(DeuteronLambda->GetPair(), AntideuteronAntilambda->GetPair(), 0.004, true);
  AntideuteronAntilambda->FixShift(AntideuteronAntilambda->GetPair(), DeuteronLambda->GetPair(), 0.004, true);


  std::cout << ""			   << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << ""			   << std::endl;

  DeuteronLambda->ReweightMixedEvent(DeuteronLambda->GetPairFixShifted(0), 0.2, 0.9);
  AntideuteronAntilambda->ReweightMixedEvent(AntideuteronAntilambda->GetPairFixShifted(0), 0.2, 0.9);

  DeuteronLambda->ReweightMixedEvent(DeuteronLambda->GetPairFixShifted(1), 0.2, 0.9);
  AntideuteronAntilambda->ReweightMixedEvent(AntideuteronAntilambda->GetPairFixShifted(1), 0.2, 0.9);

  std::cout << ""			   << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << ""			   << std::endl;

  TString foldername = file;
  foldername.ReplaceAll("AnalysisResults.root", "");

  CF_ProtonProton->SetPairs(ProtonProton, AntiprotonAntiproton);
  CF_ProtonDeuteron->SetPairs(ProtonDeuteron, AntiprotonAntideuteron);
  CF_ProtonLambda->SetPairs(ProtonLambda, AntiprotonAntilambda);
  CF_DeuteronDeuteron->SetPairs(DeuteronDeuteron, AntideuteronAntideuteron);
  CF_DeuteronLambda->SetPairs(DeuteronLambda, AntideuteronAntilambda);
  CF_LambdaLambda->SetPairs(LambdaLambda, AntilambdaAntilambda);
  
  CF_ProtonProton->GetCorrelations();
  CF_ProtonDeuteron->GetCorrelations();
  CF_ProtonLambda->GetCorrelations();
  CF_DeuteronDeuteron->GetCorrelations();
  CF_DeuteronLambda->GetCorrelations();
  CF_LambdaLambda->GetCorrelations();

  CF_ProtonProton->WriteOutput(Form("%s1_ProtonProton.root", foldername.Data()));
  CF_ProtonDeuteron->WriteOutput(Form("%s2_ProtonDeuteron.root", foldername.Data()));
  CF_ProtonLambda->WriteOutput(Form("%s3_ProtonLambda.root", foldername.Data()));
  CF_DeuteronDeuteron->WriteOutput(Form("%s4_DeuteronDeuteron.root", foldername.Data()));
  CF_DeuteronLambda->WriteOutput(Form("%s5_DeuteronLambda.root", foldername.Data()));
  CF_LambdaLambda->WriteOutput(Form("%s6_LambdaLambda.root", foldername.Data()));

}


