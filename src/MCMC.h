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
#include "IndependentModule.h"
#include <stdlib.h>
#include <random>
#include <iostream>

class MCMC {
public:
	MCMC();
	~MCMC();

	//only this array has index larger than Hyper Group Types
	char * acAssociationTypesForAllVariants;
	bool * abSelectedStatusForAllVariants;
	map<char, DependentModule*> mcpDependentAssociation;
	/**
	 * Level 1: association type
	 * Level 2: ids of variants with that association type
	 */
	// vector<vector<UINT32>> vviIndependentVariants;
	vector<UINT32> viIndependentVariants;
	UINT32* iFrequencyForAllVairants;

	static UINT64 iNumBurnins;
	static UINT64 iNumIterations;
	static UINT32 iNumVariants;
	static UINT32 iNumGroups;
	static UINT32 iNumTotalSamples;
	static UINT32 iVariantTypes;
	static UINT64 STEP;
	// {{1,2},{3}}
	static vector<vector<vector<int>>> vvviHyperGroup;
	// number of hyper group types: 3 groups, 5 hyper group types
	static UINT32 iNumHyperGroupTypes;
	static GwasData * gwasData;
	/**
	 * log(prior), the index 0 is the one showing no group difference
	 * len(dPriorIndependent) = len(dPriorDependent) = number of hyper groups
	 */
	static double * dPriorIndependent;
	static double * dPriorDependent;
	// the posterior probability for each variant
	static IndependentModule * oIndependentModule;
	// user-defined number of independent variants
	static UINT32 iPriorNumIndependentVariants;
	// user-defined number of dependent Variants
	static UINT32 iPriorNumDependentVariants;
	// number maximum dependent variants
	static UINT32 iMaxVariants;

	double dLogPosterior;

	static bool destroyMCMC();
	void getTopDependentVariant(vector<UINT32> & _viVariant, UINT32 _iNum);
	void getTopIndependentVariant(vector<UINT32> & _viVariant, UINT32 _iNum);
	bool initilizeMCMC();
	bool isConveraged(double _dLogPosterior);
	void recordFrequency();
	bool runMCMC();
	static bool setupMCMC(GwasData * _gwasData);
	void shiftPosition();
	bool test();
	bool updateDependentIndependent(double _p0, double _p1, double _p2, double _T, UINT32 _iMinDependentVariants);
	bool updateDependentIntra(double _T);
	bool updateDependentNoise(double _p0, double _p1, double _p2, double _T, UINT32 _iMinDependentVariants);
	bool updateIndependentIntra(double _T);
	bool updateIndependentNoise(double _p0, double _p1, double _p2, double _T);

private:

	double dTemperature;
	double dPosterior;
	double dProbabilitySwitch;
	double dProbabilitySwitchIn;
	double dProbabilitySwitchOut;
	double dProbabilityUpdateDependentIndependent;
	double dProbabilityUpdateDependentIntra;
	double dProbabilityUpdateDependentNoise;
	double dProbabilityUpdateIndependentIntra;
	double dProbabilityUpdateIndependentNoise;

	UINT32 iIndexTypeDependentAssoication;
	UINT32 iIndexTypeDependentAssoicationSwitch;
	UINT32 iNumSelectVariants;
	UINT32 iNumSelectVariantsSwitch;
	UINT32 iNumRestVariants;
	UINT32 iNumNewTotalVariants;
	UINT32 iSelectedVariantInner;
	UINT32 iSelectedVariantInnerSwitch;
	UINT32 iSelectedVariantOuter;
	UINT32 iSelectedVariantOuterSwitch;
	double dProbabilityMove;
	double dPosteriorDifference;
	double dProbabilityRandom;
	vector<UINT32> viInnerSelectedVariantIds;
	vector<UINT32> viInnerSelectedVariantIdsSwitch;
	vector<UINT32> viOuterSelectedVariantIds;
	vector<UINT32> viOuterSelectedVariantIdsSwitch;
	vector<UINT32> viOuterUnselectedVariantIds;
	vector<UINT32> viShuffleSelect;
	vector<char> vcAssociationTypes;

	UINT32 iNumLimitDependentVariants;
	UINT32 iNumTrueLimit;

	double logratio;

	UINT32 iNumDependentVariants;
	UINT32 iNumIndependentVariants;

	uniform_real_distribution<double> distribution;
	default_random_engine generator;
	double dRandomNumber;

	DependentModule* pDependentModule;
	DependentModule* pDependentModuleSwitch;

	int iMiniDistance;

	vector<double> vdLogPosterior;
	UINT32 iLogPosteriorSize;
	UINT32 iLogPosteriorIndex;

	vector<UINT32> viNearbyVariants;

	bool isTooClose(UINT32 _iOuterId, vector<UINT32> & _viSelectedVariants);
	void get_close_variant(UINT32 _iOuterId, vector<UINT32> & _viExistVariants, vector<UINT32> & _viCandidateVariants);
	double get_random_number();
	UINT32 get_random_number(UINT32 _iMax);
	void shuffleSelect(vector<UINT32> & _viShuffleSelect, UINT32 _iNumSelect, vector<UINT32> & _viInnerSelectedVariantIds);
	void shuffleSelect(UINT32 _iNum, UINT32 _iNumSelect, vector<UINT32> & _viInnerSelectedVariantIds);

};

#endif /* MCMC_H_ */
