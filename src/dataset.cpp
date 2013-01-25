//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "dataset.h"

CDataset::CDataset()
{
    fHasOffset = false;
    adX = NULL;
    aiXOrder = NULL;
    adXTemp4Order = NULL;
    adY = NULL;
    adOffset = NULL;
    adWeight = NULL;
    apszVarNames = NULL;

    cRows = 0;
    cCols = 0;
}


CDataset::~CDataset()
{
}




GBMRESULT CDataset::ResetWeights()
{
    GBMRESULT hr = GBM_OK;
    int i = 0;

    if(adWeight == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    for(i=0; i<cRows; i++)
    {
        adWeight[i] = 1.0;
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



GBMRESULT CDataset::SetData
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
    GBMRESULT hr = GBM_OK;

    if((adX == NULL) || (adY == NULL))
    {
        hr = GBM_INVALIDARG;
        goto Error;
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

Cleanup:
   return hr;
Error:
    goto Cleanup;
}




