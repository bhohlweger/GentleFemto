/*
 * DreamPair.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamDist.h"

#include <iostream>
DreamDist::DreamDist()
:fSE(nullptr)
,fSEMult(nullptr)
,fME(nullptr)
,fMEMult(nullptr)
{
}
DreamDist::DreamDist(DreamDist* pair, const char* name )
:fSE(nullptr)
,fSEMult(nullptr)
,fME(nullptr)
,fMEMult(nullptr)
{
  this->SetSEDist(pair->GetSEDist(),name);
  this->SetSEMultDist(pair->GetSEMultDist(),name);
  this->SetMEDist(pair->GetMEDist(),name);
  this->SetMEMultDist(pair->GetMEMultDist(),name);
}

DreamDist::~DreamDist()
{
}

