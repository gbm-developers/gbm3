//  GBM by Greg Ridgeway  Copyright (C) 2003

#include <algorithm>
#include "dataset.h"

CDataset::CDataset()
{
    fHasOffset = false;
    adX = NULL;
    aiXOrder = NULL;
    adY = NULL;
    adOffset = NULL;
    adWeight = NULL;

    cRows = 0;
    cCols = 0;
}


CDataset::~CDataset()
{
}


void CDataset::SetData
(
    double *adX,
    int *aiXOrder,
    double *adY,
    double *adOffset,
    double *adWeight,
    double *adMisc,
    int cRows,
    int cCols,
    int *acVarClasses,
    int *alMonotoneVar
)
{
  if (!(adX && adY)) {
    throw GBM::invalid_argument();
  }
  
  this->cRows = cRows;
  this->cCols = cCols;
  
  this->adX = adX;
  this->aiXOrder = aiXOrder;
  this->adY = adY;
  this->adOffset = adOffset;
  this->adWeight = adWeight;
  this->acVarClasses = acVarClasses;
  this->alMonotoneVar = alMonotoneVar;
  
  if((adOffset != NULL) && !ISNA(*adOffset))
    {
      this->adOffset = adOffset;
      fHasOffset = true;
    }
  else
    {
      this->adOffset = NULL;
      fHasOffset = false;
    }
  if((adMisc != NULL) && !ISNA(*adMisc))
    {
      this->adMisc = adMisc;
    }
  else
    {
      this->adMisc = NULL;
    }
}




