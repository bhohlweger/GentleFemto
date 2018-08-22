/*
 * DreamCF.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamCF.h"
#include "TCanvas.h"
#include <iostream>
DreamCF::DreamCF()
:fPair(nullptr)
,fPairShifted(nullptr)
,fPairLast(nullptr)
{
}

DreamCF::~DreamCF()
{
}

void DreamCF::ShiftForEmpty(DreamPair* pair) {
  fPairShifted=new DreamPair("fted");
  if (!fPairLast)
  {
    fPairLast=fPairShifted;
  }
  TH1F* SE=fPair->GetSEDist();
  TH1F* ME=fPair->GetMEDist();

  TH2F* SEMult=fPair->GetSEMultDist();
  TH2F* MEMult=fPair->GetMEMultDist();
  const int SEkMinBin=SE->FindFirstBinAbove(0);
  //+1 since we start counting bins at 1
  const int nBinsEffSE=SE->GetNbinsX()-SEkMinBin+1;
  const float SEkMin=SE->GetXaxis()->GetBinLowEdge(SEkMinBin);
  const float SEkMax=SE->GetXaxis()->GetBinUpEdge(SE->GetNbinsX());
  const int multBins=SEMult->GetNbinsY();
  const int multMax=SEMult->GetYaxis()->GetBinUpEdge(multBins);

  const char* SEHistName=Form("%sShi",SE->GetName());
  TH1F* SEShifted=new TH1F(SEHistName,SEHistName,nBinsEffSE,SEkMin,SEkMax);
  SEShifted->Sumw2();
  const char* MEHistName=Form("%sShi",ME->GetName());
  TH1F* MEShifted=new TH1F(MEHistName,MEHistName,nBinsEffSE,SEkMin,SEkMax);
  MEShifted->Sumw2();
  const char* SEMultHistName=Form("%sShi",SEMult->GetName());
  TH2F* SEMultShifted=new TH2F(SEMultHistName,SEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,1,multMax);
  SEMultShifted->Sumw2();
  const char* MEMultHistName=Form("%sShi",MEMult->GetName());
  TH2F* MEMultShifted=new TH2F(MEMultHistName,MEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,1,multMax);
  MEMultShifted->Sumw2();
  int ckBin=1;
  for (int ikBin=SEkMinBin;ikBin<=SE->GetNbinsX();++ikBin) {
    SEShifted->SetBinContent(ckBin,SE->GetBinContent(ikBin));
    SEShifted->SetBinError(ckBin,SE->GetBinError(ikBin));

    MEShifted->SetBinContent(ckBin,ME->GetBinContent(ikBin));
    MEShifted->SetBinError(ckBin,ME->GetBinError(ikBin));
    for (int iMult=1;iMult<=multBins;++iMult)
    {
      SEMultShifted->SetBinContent(ckBin,iMult,SEMult->GetBinContent(ikBin,iMult));
      SEMultShifted->SetBinError(ckBin,iMult,SEMult->GetBinError(ikBin,iMult));

      MEMultShifted->SetBinContent(ckBin,iMult,MEMult->GetBinContent(ikBin,iMult));
      MEMultShifted->SetBinError(ckBin,iMult,MEMult->GetBinError(ikBin,iMult));
    }
    ckBin++;
  }
  fPairShifted->SetSEDist(SEShifted);
  fPairShifted->SetSEMultDist(SEMultShifted);
  fPairShifted->SetMEDist(MEShifted);
  fPairShifted->SetMEMultDist(MEMultShifted);
//  TH1F* SEProjMult[multBins];
//  TH1F* SEProjMultShifted[multBins];
//  TH1F* MEProjMult[multBins];
//  TH1F* MEProjMultShifted[multBins];
//  const char* can1Name=Form("c1%s",SE->GetName());
//  TCanvas *c1 = new TCanvas(can1Name,can1Name,2000,1000);
//  c1->Divide(2,2);
//  c1->cd(1);
//  SEMult->Draw();
//  c1->cd(2);
//  SEMultShifted->Draw();
//  c1->cd(3);
//  MEMult->Draw();
//  c1->cd(4);
//  MEMultShifted->Draw();
//
//  const char* can2Name=Form("c2%s",SE->GetName());
//  TCanvas *c2 = new TCanvas(can2Name,can2Name,2000,1000);
//  c2->Divide(7,4);

//  for (int iMult=1;iMult<=multBins;++iMult)
//  {
//    TString projSEMultName=Form("%sProj%i",SEMult->GetName(),iMult);
//    SEProjMult[iMult-1]=(TH1F*)SEMult->ProjectionX(projSEMultName.Data(),iMult,iMult);
////    TString projSEMultShiftedName=Form("%sProj%i",SEMultShifted->GetName(),iMult);
////    SEProjMultShifted[iMult-1]=(TH1F*)SEMultShifted->ProjectionX("",iMult,iMult);
//    SEProjMult[iMult-1]->Divide((TH1F*)SEMultShifted->ProjectionX("",iMult,iMult));
//    c1->cd(iMult);
//    SEProjMult[iMult-1]->Draw("hist");
//
//    TString projMEMultName=Form("%sProj%i",MEMult->GetName(),iMult);
//    MEProjMult[iMult-1]=(TH1F*)MEMult->ProjectionX(projMEMultName.Data(),iMult,iMult);
////    TString projMEMultShiftedName=Form("%sProj%i",MEMultShifted->GetName(),iMult);
////    MEProjMultShifted[iMult-1]=(TH1F*)MEMultShifted->ProjectionX("",iMult,iMult);
//    MEProjMult[iMult-1]->Divide((TH1F*)MEMultShifted->ProjectionX("",iMult,iMult));
//    c2->cd(iMult);
//    MEProjMult[iMult-1]->Draw("hist");
//  }
}
