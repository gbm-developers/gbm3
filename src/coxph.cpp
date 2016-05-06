//-----------------------------------
//
// File: coxph.cpp
//
// Description: Cox partial hazards model.
//
//-----------------------------------

//-----------------------------------
// Definitions
//-----------------------------------
#define frac .00000001
#define recenter 50

//-----------------------------------
// Includes
//-----------------------------------
#include "coxph.h"
#include "origCoxState.h"
#include "censoredCoxState.h"
#include "countingCoxState.h"
#include <Rinternals.h>
#include <math.h>

//----------------------------------------
// Function Members - Private
//----------------------------------------
CCoxPH::CCoxPH(double* stats, int* sortedEnd, int* sortedSt, int* strats, bool isStartStop, int tiedMethod):
startStopCase(isStartStop), sortedEndTimes(sortedEnd), sortedStartTimes(sortedSt), strata(strats)
{
	status = stats;
	tiedTimesMethod = tiedMethod;

	// Set up which methods CoxPh will use
	if(tiedMethod >= 0)
	{
		if(startStopCase)
		{
			coxStateMethods = new CountingCoxState(this);
		}
		else
		{
			coxStateMethods = new CensoredCoxState(this);
		}
	}
	else
	{
		coxStateMethods = new OrigCoxState(this);
	}

}


//----------------------------------------
// Function Members - Public
//----------------------------------------
CDistribution* CCoxPH::Create(DataDistParams& distParams)
{

	// Initialize variables to pass to constructor
	double* stat = 0;
	int* sortedSt = NULL;
	int* sortedEnd = NULL;
	bool isStartStop = false;
	int tiesMethod;

	// Switch on misc to set up ties method
	std::string miscString = Rcpp::as<std::string>(distParams.misc[0]);
	if(miscString == "effron") 
	{
		tiesMethod = 1;
	}
	else if(miscString == "breslow")
	{
	  tiesMethod = 0;
	}
	else
	{
	  tiesMethod = -1; 
	}
	
	// Set up strata
	Rcpp::IntegerVector strats(distParams.strata);

	// Check if start/stop case or not
	Rcpp::IntegerMatrix sortMatrix(distParams.sorted);
	if(distParams.respY.ncol() > 2)
	{
		isStartStop=true;
		stat = distParams.respY(Rcpp::_, 2).begin();
		sortedEnd = sortMatrix(Rcpp::_, 1).begin();
		sortedSt = sortMatrix(Rcpp::_, 0).begin();

		return new CCoxPH(stat, sortedEnd, sortedSt, strats.begin(), isStartStop, tiesMethod);

	}

	// If not start/stop
	stat = distParams.respY(Rcpp::_, 1).begin();
	sortedEnd = sortMatrix(Rcpp::_, 0).begin();


	return new CCoxPH(stat, sortedEnd, sortedSt, strats.begin(), isStartStop, tiesMethod);


}

CCoxPH::~CCoxPH()
{
	delete coxStateMethods;
}

void CCoxPH::ComputeWorkingResponse
(
	const CDataset* pData,
    const double *adF,
    double *adZ
)
{
    coxStateMethods->ComputeWorkingResponse(pData, adF, adZ);
}

double CCoxPH::InitF
(
	const CDataset* pData
)
{
    return 0.0;
}


double CCoxPH::Deviance
(
	const CDataset* pData,
    const double *adF,
    bool isValidationSet
)
{
    // Set size and move to validation set if necessary
    long cLength = pData->get_trainSize();
	if(isValidationSet)
	{
		pData->shift_to_validation();
		status = shift_ptr(status, pData->get_trainSize());
		sortedEndTimes = shift_ptr(sortedEndTimes, pData->get_trainSize());
		sortedStartTimes = shift_ptr(sortedStartTimes, pData->get_trainSize());
		strata = shift_ptr(strata, pData->get_trainSize());
		cLength = pData->GetValidSize();
	}

	double returnValue = 0.0;
	returnValue = coxStateMethods->Deviance(cLength, pData, adF);

	// Shift back for future calculations if required
	if(isValidationSet)
	{
		pData->shift_to_train();
		status = shift_ptr(status, -(pData->get_trainSize()));
		sortedEndTimes = shift_ptr(sortedEndTimes, -(pData->get_trainSize()));
		sortedStartTimes = shift_ptr(sortedStartTimes, -(pData->get_trainSize()));
		strata = shift_ptr(strata, -(pData->get_trainSize()));
	}


	return returnValue;

}

void CCoxPH::FitBestConstant
(
	const CDataset* pData,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps* pTreeComps
)
{
   coxStateMethods->FitBestConstant(pData, adF, cTermNodes, adZ, pTreeComps);
}

double CCoxPH::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const bag& afInBag,
  const double shrinkage,
  const double* adFadj
)
{
    return coxStateMethods->BagImprovement(data, adF, afInBag, shrinkage, adFadj);
}

//-----------------------------------
// Function: GetStatusVec
//
// Returns: ptr to double vector
//
// Description: gets pointer to first element of status vector
//
// Parameters: none
//
//-----------------------------------
double* CCoxPH::StatusVec()
{
	return &status[0];
}
const double* CCoxPH::StatusVec() const
{
	return &status[0];
}

//-----------------------------------
// Function: EndTimeIndices
//
// Returns: ptr to int vector
//
// Description: gets pointer to first element of
//	the vector containing the indices which sort
//  the patient's end times.
//
// Parameters: none
//
//-----------------------------------
int* CCoxPH::EndTimeIndices()
{
	return &sortedEndTimes[0];
}
const int* CCoxPH::EndTimeIndices() const
{
	return &sortedEndTimes[0];
}

//-----------------------------------
// Function: StartTimeIndices
//
// Returns: ptr to int vector
//
// Description: gets pointer to first element of
//	the vector containing the indices which sort
//  the patient's start times.
//
// Parameters: none
//
//-----------------------------------
int* CCoxPH::StartTimeIndices()
{
	return &sortedStartTimes[0];
}
const int* CCoxPH::StartTimeIndices() const
{
	return &sortedStartTimes[0];
}

//-----------------------------------
// Function: StrataVec
//
// Returns: ptr to int vector
//
// Description: gets pointer to first element of
//	the vector containing the number of patients
//  in each strata.
//
// Parameters: none
//
//-----------------------------------
int* CCoxPH::StrataVec()
{
	return &strata[0];
}
const int* CCoxPH::StrataVec() const
{
	return &strata[0];
}


//-----------------------------------
// Function: TieApproxMethod
//
// Returns: ptr to int vector
//
// Description: gets int defining which approx method
//  to use in event of tied times.
//
// Parameters: none
//
//-----------------------------------
int CCoxPH::TieApproxMethod()
{
	return tiedTimesMethod;
}
int CCoxPH::TieApproxMethod() const
{
	return tiedTimesMethod;
}





