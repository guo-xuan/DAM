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

bool setupConfig(GwasData * gwasData) {
	Config::iNumVariants = gwasData->iNumVariants;
	Config::iNumGroups = gwasData->iNumGroups;
	return true;
}

int main(int argc, char **argv) {

	// parse command line options:
	if(!Config::setConfig(argc, argv)) {
		cout << "Please follow the above help information." << endl;
		return false;
	}

	GwasData * gwasData = new GwasData();
	gwasData->readInput(Config::vsInputVariantFiles);
	setupConfig(gwasData);

	//string sFileOut = "test.txt";
	//gwasData->writeOutput(sFileOut);

	// set the gwas data in the config and also MCMC
	MCMC::setupMCMC(gwasData);
	// create chains
	vector<MCMC *> vMcmc;
	for(size_t i = 0;i < Config::iNumChains;i++) {
		MCMC * pMcmc = new MCMC();
		// pMcmc->initilizeMCMC();
		pMcmc->test();
		vMcmc.push_back(pMcmc);
	}
	//run chains

	//free memory
	delete gwasData;
	for(size_t i = 0;i < Config::iNumChains;i++) {
		delete vMcmc.at(i);
	}
	MCMC::destroyMCMC();
	return 0;
}
