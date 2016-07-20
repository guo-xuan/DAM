//============================================================================
// Name        : DAM.cpp
// Author      : Xuan Guo
// Version     :
// Copyright   : Your copyright notice
// Description : Genome-wide Association Mapping for Multiple Groups/Populations
//============================================================================

#include <iostream>

//DAM header
#include "Config.h"
#include "GwasData.h"
#include "MCMC.h"

using namespace std;

bool setConfiguration(GwasData * gwasData) {
	Config::iNumVariants = gwasData->iNumVariants;
	Config::iNumGroups = gwasData->iNumGroups;
	MCMC::iNumVariants = Config::iNumVariants;
	MCMC::iNumGroups = Config::iNumGroups;
	MCMC::iNumBurnins = MCMC::iNumVariants * 400;
	MCMC::iNumIterations = MCMC::iNumVariants * (MCMC::iNumVariants > 100000 ? 100000 : MCMC::iNumVariants);
	MCMC::setupMCMC(Config::iNumIndependentVariants, Config::iNumDependentVariants);
	return true;
}

int main(int argc, char **argv) {

	// parse command line options:
	if (!Config::setConfig(argc, argv)) {
		cout << "Please follow the above help information." << endl;
		return false;
	}

	GwasData * gwasData = new GwasData();
	gwasData->readInput(Config::vsInputVariantFiles);
	setConfiguration(gwasData);
	//string sFileOut = "test.txt";
	//gwasData->writeOutput(sFileOut);
	//create chains
	vector<MCMC *> vMcmc;
	for (size_t i = 0; i < Config::iNumChains; i++) {
		MCMC * pMcmc = new MCMC();
		vMcmc.push_back(pMcmc);
	}
	//run chains

	//free memory
	delete gwasData;
	for (size_t i = 0; i < Config::iNumChains; i++) {
		delete vMcmc.at(i);
	}
	return 0;
}
