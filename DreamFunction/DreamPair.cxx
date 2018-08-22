/*
 * DreamCF.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamPair.h"

#include "TCanvas.h"
#include <iostream>
DreamPair::DreamPair()
:fPair(nullptr)
,fPairShifted()
,fPairRebinned()
,fPairReweighted()
,fFirstBin(-99)
{
}

DreamPair::~DreamPair()
{
}

void DreamPair::ShiftForEmpty(DreamDist* pair)
{
  DreamDist* PairShifted=new DreamDist();

  TH1F* SE=pair->GetSEDist();
  TH1F* ME=pair->GetMEDist();

  TH2F* SEMult=pair->GetSEMultDist();
  TH2F* MEMult=pair->GetMEMultDist();
  fFirstBin=SE->FindFirstBinAbove(0);
  //+1 since we start counting bins at 1
  const int nBinsEffSE=SE->GetNbinsX()-fFirstBin+1;
  const float SEkMin=SE->GetXaxis()->GetBinLowEdge(fFirstBin);
  const float SEkMax=SE->GetXaxis()->GetBinUpEdge(SE->GetNbinsX());
  const int multBins=SEMult->GetNbinsY();
  const int multMax=SEMult->GetYaxis()->GetBinUpEdge(multBins);

  const char* SEHistName=Form("%s_",SE->GetName());
  TH1F* SEShifted=new TH1F(SEHistName,SEHistName,nBinsEffSE,SEkMin,SEkMax);
  SEShifted->Sumw2();
  const char* MEHistName=Form("%s_",ME->GetName());
  TH1F* MEShifted=new TH1F(MEHistName,MEHistName,nBinsEffSE,SEkMin,SEkMax);
  MEShifted->Sumw2();
  const char* SEMultHistName=Form("%s_",SEMult->GetName());
  TH2F* SEMultShifted=new TH2F(SEMultHistName,SEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,1,multMax);
  SEMultShifted->Sumw2();
  const char* MEMultHistName=Form("%s_",MEMult->GetName());
  TH2F* MEMultShifted=new TH2F(MEMultHistName,MEMultHistName,
                               nBinsEffSE,SEkMin,SEkMax,
                               multBins,1,multMax);
  MEMultShifted->Sumw2();
  int ckBin=1;
  for (int ikBin=fFirstBin;ikBin<=SE->GetNbinsX();++ikBin) {
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
  PairShifted->SetSEDist(SEShifted,"Shifted");
  PairShifted->SetSEMultDist(SEMultShifted,"Shifted");
  PairShifted->SetMEDist(MEShifted,"Shifted");
  PairShifted->SetMEMultDist(MEMultShifted,"Shifted");

  fPairShifted.push_back(PairShifted);
  delete SEShifted;
  delete SEMultShifted;
  delete MEShifted;
  delete MEMultShifted;
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

void DreamPair::Rebin(DreamDist* pair, int rebin)
{
  DreamDist* Rebinned=new DreamDist(pair,Form("_Rebinned_%i",rebin));
  Rebinned->GetSEDist()->Rebin(rebin);
  Rebinned->GetSEMultDist()->Rebin2D(rebin,1);
  Rebinned->GetMEDist()->Rebin(rebin);
  Rebinned->GetMEMultDist()->Rebin2D(rebin,1);

  fPairRebinned.push_back(Rebinned);
}

void DreamPair::ReweightMixedEvent(DreamDist* pair,float kSMin,float kSMax)
{
  DreamDist* PairReweighted=new DreamDist(pair,"_Reweighted");

  TH1F* SE=pair->GetSEDist();
  TH1F* ME=pair->GetMEDist();
  TH2F* SEMult=pair->GetSEMultDist();
  TH2F* MEMult=pair->GetMEMultDist();

  TH1F* MEReweighted=PairReweighted->GetMEDist();
  TH2F* MEMultReweighted=PairReweighted->GetMEMultDist();
  MEReweighted->Reset(); // this one we fill from scratch, we only want to keep the dimensions
  MEMultReweighted->Reset();

  int nKSbins=SEMult->GetXaxis()->GetNbins();
  double kSMaxVal=SEMult->GetXaxis()->GetBinUpEdge(nKSbins);
  //k* range of the normalization of ME and SE in Multiplicity
  int firstBin=SEMult->FindBin(kSMin);
  int lastBin=SEMult->FindBin(kSMax);
  TH1D* MultProjSE=SEMult->ProjectionY(Form("%skSProj",SEMult->GetName()),firstBin,lastBin);
  TH1D* MultProjME=MEMult->ProjectionY(Form("%skSProj",MEMult->GetName()),firstBin,lastBin);

  int nMultBins=MEMult->GetYaxis()->GetNbins();
  int multMax=MEMult->GetYaxis()->GetBinUpEdge(nMultBins);
  double weight=0;
  for (int iMult = 1;iMult<=nMultBins;++iMult) {
    if (MultProjME->GetBinContent(iMult)>0) {
      weight=MultProjSE->GetBinContent(iMult)/MultProjME->GetBinContent(iMult);
    } else {
      weight=0;
      std::cout << "Weight = 0 for Mult Bin iMult = "
          << iMult  << " in the case of " << SE->GetName() << std::endl;
    }
    TString MultBinName=Form("MEBin%i",iMult);
    MEReweighted->Add(MEMult->ProjectionX(MultBinName.Data(),iMult,iMult),weight);

    for (int ikStar=1;ikStar<=MEMult->GetNbinsX();++ikStar)
    {
      MEMultReweighted->SetBinContent(ikStar,iMult,MEMult->GetBinContent(ikStar,iMult)*weight);
      MEMultReweighted->SetBinError(ikStar,iMult,MEMult->GetBinError(ikStar,iMult)*weight);
    }
  }
  fPairReweighted.push_back(PairReweighted);
//  delete SEReweighted;
//  delete SEMultReweighted;
//  delete MEReweighted;
//  delete MEMultReweighted;
  delete MultProjSE;
  delete MultProjME;

//  const char* can1Name=Form("c1%s",SE->GetName());
//  TCanvas *c1 = new TCanvas(can1Name,can1Name,2000,1000);
//  c1->Divide(2,2);
//  c1->cd(1);
//  PairReweighted->GetSEMultDist()->Draw("COLZ");
//  c1->cd(2);
//  SEMult->Draw("COLZ");
//  c1->cd(3);
//  MEMultReweighted->Draw("COLZ");
//  c1->cd(4);
//  MEMult->Draw("COLZ");
//
//  const char* can2Name=Form("c2%s",SE->GetName());
//  TCanvas *c2 = new TCanvas(can2Name,can2Name,2000,1000);
//  c2->Divide(2,2);
//  c2->cd(1);
//  SE->Draw();
//  c2->cd(2);
//  PairReweighted->GetSEDist()->Draw("same");
//
//  c2->cd(3);
//  ME->Scale(1./ME->Integral());
//  ME->Draw();
//  c2->cd(4);
//
//  MEReweighted->SetLineColor(2);
//  MEReweighted->Scale(1./MEReweighted->Integral());
//  MEReweighted->Draw("same");
}

