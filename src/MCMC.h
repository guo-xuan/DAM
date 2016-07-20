/*
 * MCMC.h
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#ifndef MCMC_H_
#define MCMC_H_

#include "Config.h"
#include "DependentModule.h"
#include "HyperGroup.h"
#include "GwasData.h"

class MCMC {
public:
	MCMC();
	~MCMC();

	//only this array has index larger than Hyper Group Types
	char * acAssociationTypesForAllVariants;
	map<char,DependentModule*> mcpDependentAssociation;
	vector<int> viIndependentVariant;
	int** iFrequencyForAllVairants;

	static int iNumBurnins;
	static int iNumIterations;
	static int iNumVariants;
	static int iNumGroups;
	//{{1,2},{3}}
	static vector<vector<vector<int>>> vvviHyperGroup;
	static int iNumHyperGroupTypes;
	static GwasData * gwasData;
	//prior, the index 0 is the one showing no group difference
	static double * dPriorIndependent;
	static double * dPriorDependent;

	double dLogPosterior;
	double dTemperature;
	double dProbabilitySwitch;
	double dProbabilitySwitchIn;
	double dProbabilitySwitchOut;
	double dProbabilityUpdateDependentIndependent;
	double dProbabilityUpdateDependentIntra;
	double dProbabilityUpdateDependentNoise;
	double dProbabilityUpdateIndependentIntra;
	double dProbabilityUpdateIndependentNoise;

	bool initilizeMCMC();
	bool isConveraged();
	bool runMCMC();
	static bool setupMCMC(double _iNumPriorIndependentVariants, double _iNumPriorDependentVariants);
	bool updateDependentIndependent();
	bool updateDependentIntra();
	bool updateDependentNoise();
	bool updateIndependentIntra();
	bool updateIndependentNoise();

};

#endif /* MCMC_H_ */
