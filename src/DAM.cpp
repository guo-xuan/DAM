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
#include "AssociationEvaluation.h"
#include "GammaXG.h"

using namespace std;

bool setupConfig(GwasData * gwasData) {
	Config::iNumVariants = gwasData->iNumVariants;
	Config::iNumGroups = gwasData->iNumGroups;
	return true;
}

int main(int argc, char **argv) {

/*	vector<UINT32> vCombination;
	vCombination.push_back(0);
	vCombination.push_back(1);
	vCombination.push_back(1);
	vector<Interaction> vInteraction;
	AssociationEvaluation obj;
	obj.combinationGenerator(1, 20, vInteraction, 1000, vCombination);
	return 0;*/
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
		pMcmc->initilizeMCMC();
		// pMcmc->test();
		vMcmc.push_back(pMcmc);
	}
	//run chains
	for(size_t i = 0;i < Config::iNumChains;i++) {
		vMcmc.at(i)->runMCMC();
	}
	// sum up all the frequencies from all chains, choose the top 200 variants for next phase
	// vector<UINT32> viVariants = {1, 2, 3, 4, 5, 6};

	//free memory
	for(size_t i = 0;i < Config::iNumChains;i++) {
		delete vMcmc.at(i);
	}

	MCMC::destroyMCMC();

	vector<UINT32> viVariants = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	// step-wise evaluation
	AssociationEvaluation * pAssociationEvaluation = new AssociationEvaluation();
	pAssociationEvaluation->initialize(viVariants, gwasData, Config::iMaxVariantsInInteraction);
	pAssociationEvaluation->Evaluation();
	string file_result = Config::sOutputFolder + "result.txt";
	pAssociationEvaluation->WriteResults(file_result);
	// delete pAssociationEvaluation;
	// delete gwasData;
	cout << "Done." << endl;
	return 0;
}
