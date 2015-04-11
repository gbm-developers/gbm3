//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       dataset.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Dataset class
//
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef DATASET_H
#define DATASET_H

#include <vector>

typedef std::vector<int> bag;

#include "buildinfo.h"
#include "gbmexcept.h"

class CDataset
{
public:
	CDataset();
	~CDataset();

	void SetData(double *adX,
		     int *aiXOrder,
		     double *adY,
		     double *adOffset,
		     double *adWeight,
		     double *adMisc,
		     int cRows,
		     int cCols,
		     int *acVarClasses,
		     int *alMonotoneVar);

    bool fHasOffset;
    double *adX;
    int *aiXOrder;

    double *adY;
    double *adOffset;
    double *adWeight;
    double *adMisc;
    int *acVarClasses;
    int *alMonotoneVar;

    int cRows;
    int cCols;
private:

};

#endif // DATASET_H


