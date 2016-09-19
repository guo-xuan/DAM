/*
 * MCMC.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "MCMC.h"

int MCMC::iNumVariants = 0;
int MCMC::iNumBurnins = 0;
int MCMC::iNumIterations = 0;
int MCMC::iNumHyperGroupTypes = 0;
int MCMC::iNumGroups = 0;
vector<vector<vector<int>>> MCMC::vvviHyperGroup;
GwasData * MCMC::gwasData = NULL;
double * MCMC::dPriorDependent = NULL;
double * MCMC::dPriorIndependent = NULL;
IndependentModule * MCMC::oIndependentModule = NULL;
UINT32 MCMC::iNumIndependentVariants = 0;
UINT32 MCMC::iNumDependentVariants = 0;
UINT32 MCMC::iMaxVariants = 0;

MCMC::MCMC() {
	dLogPosterior = 0;
	dTemperature = 1;
	dProbabilitySwitch = 0.33;
	dProbabilitySwitchIn = 0.34;
	dProbabilitySwitchOut = 0.33;
	dProbabilityUpdateDependentIndependent = 0.3;
	dProbabilityUpdateDependentIntra = 0.1;
	dProbabilityUpdateDependentNoise = 0.3;
	dProbabilityUpdateIndependentIntra = 0.2;
	dProbabilityUpdateIndependentNoise = 0.1;
	iFrequencyForAllVairants = new int*[iNumVariants]();
	for(int i = 0;i < iNumVariants;i++) {
		//dependent and independent, so times 2
		iFrequencyForAllVairants[i] = new int[iNumHyperGroupTypes * 2]();
	}
	acAssociationTypesForAllVariants = NULL;
}

MCMC::~MCMC() {
	for(int i = 0;i < iNumVariants;i++) {
		delete[] iFrequencyForAllVairants[i];
	}
	delete[] iFrequencyForAllVairants;
	if(acAssociationTypesForAllVariants != NULL) {
		delete[] acAssociationTypesForAllVariants;
	}
}

/**
 *
 */
bool MCMC::initilizeMCMC() {
	acAssociationTypesForAllVariants = new char[iNumVariants]();
	dLogPosterior = 0;
	// randomly assign the association types to all variants
	UINT32 id = 0;
	vviIndependentVariants.push_back(vector<UINT32>());
	for(int i = 1;i < iNumHyperGroupTypes;i++) {
		vviIndependentVariants.push_back(vector<UINT32>());
		for(UINT32 j = 0;j < iNumIndependentVariants;) {
			id = (UINT32) ((((double) rand()) / ((double) RAND_MAX)) * ((double) iNumVariants));
			if(acAssociationTypesForAllVariants[id] == 0) {
				acAssociationTypesForAllVariants[id] = i;
				j++;
				vviIndependentVariants.at(i).push_back(id);
				dLogPosterior += dPriorIndependent[i];
				dLogPosterior += oIndependentModule->getIndependentPosterior(id, i);
			}
		}
	}
	vector<UINT32> viVariants;
	for(int i = 1 + iNumHyperGroupTypes;i < iNumHyperGroupTypes * 2;i++) {
		mcpDependentAssociation[i] = new DependentModule(iMaxVariants, gwasData,
				vvviHyperGroup.at(i - iNumHyperGroupTypes));
		viVariants.clear();
		for(UINT32 j = 0;j < iNumDependentVariants;) {
			id = (UINT32) ((((double) rand()) / ((double) RAND_MAX)) * ((double) iNumVariants));
			if(acAssociationTypesForAllVariants[id] == 0) {
				acAssociationTypesForAllVariants[id] = i;
				j++;
				dLogPosterior += dPriorDependent[i];
				viVariants.push_back(id);
			}
		}
		mcpDependentAssociation[i]->addVariants(viVariants);
		dLogPosterior = mcpDependentAssociation[i]->getPosterior();
	}
	return true;
}

bool MCMC::runMCMC() {

	return true;
}

/**
 * setup all the static variables in MCMC
 */
bool MCMC::setupMCMC(double _iNumPriorIndependentVariants, double _iNumPriorDependentVariants) {
	iNumIndependentVariants = _iNumPriorIndependentVariants;
	iNumDependentVariants = _iNumPriorDependentVariants;
	vvviHyperGroup = HyperGroup::getHyperGroup(iNumGroups);
	iNumHyperGroupTypes = vvviHyperGroup.size();
	dPriorIndependent = new double[iNumHyperGroupTypes];
	dPriorDependent = new double[iNumHyperGroupTypes * 2];
	for(int i = 1;i < iNumHyperGroupTypes;i++) {
		dPriorIndependent[i] = _iNumPriorIndependentVariants / ((double) iNumVariants);
		dPriorIndependent[0] += dPriorIndependent[i];
		dPriorDependent[i + iNumHyperGroupTypes] = _iNumPriorDependentVariants / ((double) iNumVariants);
		dPriorDependent[0 + iNumHyperGroupTypes] += dPriorDependent[i];
	}
	dPriorIndependent[0] = 1 - dPriorIndependent[0];
	dPriorDependent[0 + iNumHyperGroupTypes] = 1 - dPriorDependent[0 + iNumHyperGroupTypes];
	for(int i = 0;i < iNumHyperGroupTypes;i++) {
		dPriorIndependent[i] = log(dPriorIndependent[i]);
		dPriorDependent[i + iNumHyperGroupTypes] = log(dPriorDependent[i + iNumHyperGroupTypes]);
	}
	// Initialize the posterior probability for each variant
	oIndependentModule = new IndependentModule();
	oIndependentModule->initial(gwasData, vvviHyperGroup);
	iMaxVariants = Config::iMaxVariants;
	return true;
}

bool MCMC::destroyMCMC() {

	delete[] dPriorIndependent;
	delete[] dPriorDependent;
	delete oIndependentModule;

	return true;
}
