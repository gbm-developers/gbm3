//-----------------------------------
//
// File: coxph.cpp
//
// Description: Cox partial hazards model.
//
//-----------------------------------

//-----------------------------------
// Includes
//-----------------------------------
#include "coxph.h"
#include "censored_cox_state.h"
#include "counting_cox_state.h"
#include <Rinternals.h>
#include <math.h>

namespace {
  int GetTiesMethod(const std::string& selection)
  {
    if (selection == "efron")
    {
      return 1;
    }
    else if (selection == "breslow")
    {
      return 0;
    }

    throw GBM::InvalidArgument("unknown tie-handling method");
  }
}

//----------------------------------------
// Function Members - Private
//----------------------------------------
CCoxPH::CCoxPH(double* stats, int* sortedEnd, int* sortedSt, int* strats,
		bool isStartStop, int tiedMethod, double priorCoeff):
kStartStopCase_(isStartStop), sortedendtimes_(sortedEnd), sortedstarttimes_(sortedSt), strata(strats), kPriorCoeffVariation_(priorCoeff)
{
	status_ = stats;
	tiedtimesmethod_ = tiedMethod;

	// Set up which methods CoxPh will use
	if(kStartStopCase_)
	{
		coxstate_methods_ = new CountingCoxState(this);
	}
	else
	{
		coxstate_methods_ = new CensoredCoxState(this);
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
	int tiesMethod = GetTiesMethod(Rcpp::as<string>(distParams.misc[0]));

	// Switch on misc to set up ties method
	std::string miscString = Rcpp::as<std::string>(distParams.misc[0]);

	// Set up strata
	Rcpp::IntegerVector strats(distParams.strata);

	// Check if start/stop case or not
	Rcpp::IntegerMatrix sortMatrix(distParams.sorted);
	if(distParams.response.ncol() > 2)
	{
		isStartStop=true;
		stat = distParams.response(Rcpp::_, 2).begin();
		sortedEnd = sortMatrix(Rcpp::_, 1).begin();
		sortedSt = sortMatrix(Rcpp::_, 0).begin();

		return new CCoxPH(stat, sortedEnd, sortedSt, strats.begin(), isStartStop, tiesMethod, distParams.prior_coefficient_variation);

	}

	// If not start/stop
	stat = distParams.response(Rcpp::_, 1).begin();
	sortedEnd = sortMatrix(Rcpp::_, 0).begin();


	return new CCoxPH(stat, sortedEnd, sortedSt, strats.begin(), isStartStop, tiesMethod, distParams.prior_coefficient_variation);


}

CCoxPH::~CCoxPH()
{
	delete coxstate_methods_;
}

void CCoxPH::ComputeWorkingResponse
(
	const CDataset& data,
    const double *adF,
    double *adZ
)
{
    coxstate_methods_->ComputeWorkingResponse(data, adF, adZ);
}

double CCoxPH::InitF
(
	const CDataset& data
)
{
    return 0.0;
}


double CCoxPH::Deviance
(
	const CDataset& data,
    const double *adF,
    bool isValidationSet
)
{
    // Set size and move to validation set if necessary
    long cLength = data.get_trainsize();
	if(isValidationSet)
	{
		data.shift_to_validation();
		status_ = shift_ptr(status_, data.get_trainsize());
		sortedendtimes_ = shift_ptr(sortedendtimes_, data.get_trainsize());
		sortedstarttimes_ = shift_ptr(sortedstarttimes_, data.get_trainsize());
		strata = shift_ptr(strata, data.get_trainsize());
		cLength = data.get_validsize();
	}

	double returnValue = 0.0;
	returnValue = coxstate_methods_->Deviance(cLength, data, adF);

	// Shift back for future calculations if required
	if(isValidationSet)
	{
		data.shift_to_train();
		status_ = shift_ptr(status_, -(data.get_trainsize()));
		sortedendtimes_ = shift_ptr(sortedendtimes_, -(data.get_trainsize()));
		sortedstarttimes_ = shift_ptr(sortedstarttimes_, -(data.get_trainsize()));
		strata = shift_ptr(strata, -(data.get_trainsize()));
	}


	return returnValue;

}

void CCoxPH::FitBestConstant
(
	const CDataset& data,
    const double *adF,
    unsigned long cTermNodes,
    double* adZ,
    CTreeComps& treeComps
)
{
   coxstate_methods_->FitBestConstant(data, adF, cTermNodes, adZ, treeComps);
}

double CCoxPH::BagImprovement
(
	const CDataset& data,
    const double *adF,
    const double shrinkage,
  const double* adFadj
)
{
    return coxstate_methods_->BagImprovement(data, adF, shrinkage, adFadj);
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
	return &status_[0];
}
const double* CCoxPH::StatusVec() const
{
	return &status_[0];
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
	return &sortedendtimes_[0];
}
const int* CCoxPH::EndTimeIndices() const
{
	return &sortedendtimes_[0];
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
	return &sortedstarttimes_[0];
}
const int* CCoxPH::StartTimeIndices() const
{
	return &sortedstarttimes_[0];
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
// Returns: int
//
// Description: gets int defining which approx method
//  to use in event of tied times.
//
// Parameters: none
//
//-----------------------------------
int CCoxPH::TieApproxMethod() const
{
	return tiedtimesmethod_;
}

//-----------------------------------
// Function: PriorCoeffVar
//
// Returns: double
//
// Description: gets double defining the expected hazard
//   of a node with 0 events and expected events.
//
// Parameters: none
//
//-----------------------------------
double CCoxPH::PriorCoeffVar() const
{
	return kPriorCoeffVariation_;
}


