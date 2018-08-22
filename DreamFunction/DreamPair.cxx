/*
 * DreamPair.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamPair.h"
#include <iostream>
DreamPair::DreamPair(const char* name )
:fName(name)
,fSE(nullptr)
,fSEMult(nullptr)
,fME(nullptr)
,fMEMult(nullptr)
{
}
DreamPair::DreamPair(DreamPair* pair, const char* name )
:fName(name)
,fSE(nullptr)
,fSEMult(nullptr)
,fME(nullptr)
,fMEMult(nullptr)
{
  this->SetSEDist(pair->GetSEDist());
  this->SetSEMultDist(pair->GetSEMultDist());
  this->SetMEDist(pair->GetMEDist());
  this->SetMEMultDist(pair->GetMEMultDist());
}

DreamPair::~DreamPair()
{
}

