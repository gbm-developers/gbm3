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
CCoxPH::CCoxPH(double* stats, int* sorted_end, int* sorted_start, int* strats,
		bool is_startstop, int tiesmethod, double priorcoeff):
kStartStopCase_(is_startstop), sortedendtimes_(sorted_end), sortedstarttimes_(sorted_start), strata(strats), kPriorCoeffVariation_(priorcoeff)
{
	status_ = stats;
	tiedtimesmethod_ = tiesmethod;

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
CDistribution* CCoxPH::Create(DataDistParams& distparams)
{

	// Initialize variables to pass to constructor
	double* stat = 0;
	int* sortedst = NULL;
	int* sortedend = NULL;
	bool isstartstop = false;
	int tiesmethod = GetTiesMethod(Rcpp::as<string>(distparams.misc[0]));

	// Switch on misc to set up ties method
	std::string miscstring = Rcpp::as<std::string>(distparams.misc[0]);

	// Set up strata
	Rcpp::IntegerVector strats(distparams.strata);

	// Check if start/stop case or not
	Rcpp::IntegerMatrix sortMatrix(distparams.sorted);
	if(distparams.response.ncol() > 2)
	{
		isstartstop=true;
		stat = distparams.response(Rcpp::_, 2).begin();
		sortedend = sortMatrix(Rcpp::_, 1).begin();
		sortedst = sortMatrix(Rcpp::_, 0).begin();

		return new CCoxPH(stat, sortedend, sortedst, strats.begin(), isstartstop, tiesmethod, distparams.prior_coefficient_variation);

	}

	// If not start/stop
	stat = distparams.response(Rcpp::_, 1).begin();
	sortedend = sortMatrix(Rcpp::_, 0).begin();


	return new CCoxPH(stat, sortedend, sortedst, strats.begin(), isstartstop, tiesmethod, distparams.prior_coefficient_variation);


}

CCoxPH::~CCoxPH()
{
	delete coxstate_methods_;
}

void CCoxPH::ComputeWorkingResponse
(
	const CDataset& kData,
    const double* kFuncEstimate,
    double* residuals
)
{
    coxstate_methods_->ComputeWorkingResponse(kData, kFuncEstimate, residuals);
}

double CCoxPH::InitF
(
	const CDataset& kData
)
{
    return 0.0;
}


double CCoxPH::Deviance
(
	const CDataset& kData,
    const double* kFuncEstimate,
    bool is_validationset
)
{
    // Set size and move to validation set if necessary
    long num_rows_in_set = kData.get_trainsize();
	if(is_validationset)
	{
		kData.shift_to_validation();
		status_ = shift_ptr(status_, kData.get_trainsize());
		sortedendtimes_ = shift_ptr(sortedendtimes_, kData.get_trainsize());
		sortedstarttimes_ = shift_ptr(sortedstarttimes_, kData.get_trainsize());
		strata = shift_ptr(strata, kData.get_trainsize());
		num_rows_in_set = kData.get_validsize();
	}

	double returnvalue = 0.0;
	returnvalue = coxstate_methods_->Deviance(num_rows_in_set, kData, kFuncEstimate);

	// Shift back for future calculations if required
	if(is_validationset)
	{
		kData.shift_to_train();
		status_ = shift_ptr(status_, -(kData.get_trainsize()));
		sortedendtimes_ = shift_ptr(sortedendtimes_, -(kData.get_trainsize()));
		sortedstarttimes_ = shift_ptr(sortedstarttimes_, -(kData.get_trainsize()));
		strata = shift_ptr(strata, -(kData.get_trainsize()));
	}


	return returnvalue;

}

void CCoxPH::FitBestConstant
(
	const CDataset& kData,
    const double* kFuncEstimate,
    unsigned long num_terminalnodes,
    double* residuals,
    CCARTTree& tree
)
{
   coxstate_methods_->FitBestConstant(kData, kFuncEstimate, num_terminalnodes, residuals, tree);
}

double CCoxPH::BagImprovement
(
	const CDataset& kData,
    const double* kFuncEstimate,
    const double kShrinkage,
  const double* kDeltaEstimate
)
{
    return coxstate_methods_->BagImprovement(kData, kFuncEstimate, kShrinkage, kDeltaEstimate);
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


