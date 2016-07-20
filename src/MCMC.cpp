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
GwasData * MCMC::gwasData=NULL;

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
	for(int i=0;i<iNumVariants;i++){
		//dependent and independent, so times 2
		iFrequencyForAllVairants[i] = new int[iNumHyperGroupTypes*2]();
	}
}

MCMC::~MCMC() {
	for(int i=0;i<iNumVariants;i++){
		delete[] iFrequencyForAllVairants[i];
	}
	delete[] iFrequencyForAllVairants;
}

bool MCMC::initilizeMCMC(){
	acAssociationTypesForAllVariants = new char[iNumVariants]();

	return true;
}

bool MCMC::runMCMC() {



	return true;
}

bool MCMC::setupMCMC(double _iNumPriorIndependentVariants, double _iNumPriorDependentVariants){
	vvviHyperGroup = HyperGroup::getHyperGroup(iNumGroups);
	iNumHyperGroupTypes = vvviHyperGroup.size();
	dPriorIndependent = new double[iNumHyperGroupTypes];
	dPriorDependent = new double[iNumHyperGroupTypes];
	for(int i = 1;i<iNumHyperGroupTypes;i++){
		dPriorIndependent[i] = _iNumPriorIndependentVariants/((double)iNumVariants);
		dPriorIndependent[0] += dPriorIndependent[i];
		dPriorDependent[i] = _iNumPriorDependentVariants/((double)iNumVariants);
		dPriorDependent[0] += dPriorDependent[i];
	}
	return true;
}
