// Replacement for nodeParameters

#ifndef __nodeSplitParams_h__
#define __nodeSplitParams_h__

#include <Rcpp.h>

class nodeSplitParams
{
public:
	// Constructor
	nodeSplitParams()
	{
		RightNumObs = 0;
		RightTotalWeight = 0.0;
		RightWeightResiduals = 0.0;

		LeftWeightResiduals   = 0.0;
		LeftTotalWeight = 0.0;
		LeftNumObs     = 0;

		MissingWeightResiduals   = 0.0;
		MissingTotalWeight = 0.0;
		MissingNumObs     = 0;

		ImprovedResiduals = 0.0;
		SplitValue = 0.0;
		SplitClass = 0;
		SplitVar = 0;
		aiBestCategory.resize(1024);
	}

	// Destructor
	~nodeSplitParams();

	/*// Public Functions - inline
	 void ResetSplitProperties(double weightedResiduals, double trainingWeight,
								   unsigned long numObs, double splitValue, unsigned long variableClasses, unsigned long splitVar);
	 void UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement);
	 void UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement);
	 void NodeGradResiduals();*/
private:
	// Left Node Definition
	double LeftWeightResiduals;
	double LeftTotalWeight;
	long LeftNumObs;

	// Right Node Definition
	double RightWeightResiduals;
	double RightTotalWeight;
	long RightNumObs;

	// Missing Node Definition
	double MissingWeightResiduals;
	double MissingTotalWeight;
	long MissingNumObs;

	// Splitting values
	double SplitValue; // Continuous Split Value
	unsigned long SplitVar; // Which feature to split on
	unsigned long SplitClass; // Categorical Split Value
	std::vector<int> aiBestCategory; // Vector of levels ordering
	double ImprovedResiduals;
};

/*inline void nodeSplit::ResetSplitProperties(double weightedResiduals, double trainingWeight,
									   unsigned long numObs, double splitValue, unsigned long variableClasses, unsigned long splitVar)
{

		RightWeightResiduals   = weightedResiduals;
		RightTotalWeight = trainingWeight;
		RightNumObs     = numObs;


		LeftWeightResiduals   = 0.0;
		LeftTotalWeight = 0.0;
		LeftNumObs     = 0;


		MissingWeightResiduals   = 0.0;
		MissingTotalWeight = 0.0;
		MissingNumObs     = 0;


		SplitVar = splitVar;
		SplitValue = splitValue;
		ImprovedResiduals = 0.0;
		SplitClass = variableClasses;
}


inline void nodeSplit::UpdateMissingNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to missing
	MissingWeightResiduals += predIncrement;
	MissingTotalWeight += trainWIncrement;
	MissingNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;
}

inline void nodeSplit::UpdateLeftNode(double predIncrement, double trainWIncrement, long numIncrement)
{
	// Move data point from right node to left node
	LeftWeightResiduals += predIncrement;
	LeftTotalWeight += trainWIncrement;
	LeftNumObs += numIncrement;

	RightWeightResiduals -= predIncrement;
	RightTotalWeight -= trainWIncrement;
	RightNumObs -= numIncrement;

}

inline void nodeSplit::NodeGradResiduals()
{
	// Returns weighted

	double dTemp = 0.0;
	double dResult = 0.0;

	// Only need to look at left and right
	if(MissingNumObs == 0.0)
	{
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult = LeftTotalWeight*RightTotalWeight*dTemp*dTemp/(LeftTotalWeight+RightTotalWeight);
	}
	else
	{
		// Grad - left/right
		dTemp = LeftWeightResiduals/LeftTotalWeight - RightWeightResiduals/RightTotalWeight;
		dResult += LeftTotalWeight*RightTotalWeight*dTemp*dTemp;

		// Grad - left/missing
		dTemp = LeftWeightResiduals/LeftTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += LeftTotalWeight*MissingTotalWeight*dTemp*dTemp;

		// Grad - right/missing
		dTemp = RightWeightResiduals/RightTotalWeight - MissingWeightResiduals/MissingTotalWeight;
		dResult += RightTotalWeight*MissingTotalWeight*dTemp*dTemp;
		dResult /= (LeftTotalWeight + RightTotalWeight + MissingTotalWeight);
	}

	// Update current residuals
	ImprovedResiduals = dResult;
}*/
#endif // __nodeSplitParams_h__
