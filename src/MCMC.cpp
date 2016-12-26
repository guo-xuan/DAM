/*
 * MCMC.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "MCMC.h"

UINT32 MCMC::iNumVariants = 0;
UINT64 MCMC::iNumBurnins = 0;
UINT64 MCMC::iNumIterations = 0;
UINT64 MCMC::STEP = 1000;
UINT32 MCMC::iNumHyperGroupTypes = 0;
UINT32 MCMC::iNumGroups = 0;
UINT32 MCMC::iNumTotalSamples = 0;
UINT32 MCMC::iVariantTypes = 0;
vector<vector<vector<int>>> MCMC::vvviHyperGroup;
GwasData * MCMC::gwasData = NULL;
double * MCMC::dPriorDependent = NULL;
double * MCMC::dPriorIndependent = NULL;
IndependentModule * MCMC::oIndependentModule = NULL;
UINT32 MCMC::iPriorNumIndependentVariants = 0;
UINT32 MCMC::iPriorNumDependentVariants = 0;
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
	iFrequencyForAllVairants = new UINT32[iNumVariants * iNumHyperGroupTypes * 2]();
	acAssociationTypesForAllVariants = NULL;
	dPosterior = 0;
	iNumSelectVariants = 0;
	iNumRestVariants = 0;
	dPosteriorDifference = 0;
	distribution = uniform_real_distribution<double>(0.0, 1.0);
	dRandomNumber = 0;

	iLogPosteriorSize = 2000;
	iLogPosteriorIndex = 0;

	logratio = 0;
}

MCMC::~MCMC() {
	delete[] iFrequencyForAllVairants;
	if(acAssociationTypesForAllVariants != NULL) {
		delete[] acAssociationTypesForAllVariants;
	}
}

bool MCMC::destroyMCMC() {

	delete[] dPriorIndependent;
	delete[] dPriorDependent;
	delete oIndependentModule;

	return true;
}

void MCMC::getTopDependentVariant(vector<UINT32> & _viVariant, UINT32 _iNum){

}

/**
 *
 */
bool MCMC::initilizeMCMC() {
	acAssociationTypesForAllVariants = new char[iNumVariants]();
	abSelectedStatusForAllVariants = new bool[iNumVariants]();
	dLogPosterior = 0;
	// randomly assign the association types to all variants
	UINT32 id = 0;
	// vviIndependentVariants.push_back(vector<UINT32>());
	for(UINT32 i = 1;i < iNumHyperGroupTypes;i++) {
		// vviIndependentVariants.push_back(vector<UINT32>());
		for(UINT32 j = 0;j < iPriorNumIndependentVariants;) {
			id = (UINT32) ((((double) rand()) / ((double) RAND_MAX)) * ((double) iNumVariants));
			if(acAssociationTypesForAllVariants[id] == 0) {
				acAssociationTypesForAllVariants[id] = i;
				j++;
				viIndependentVariants.push_back(id);
				// vviIndependentVariants.at(i).push_back(id);
				dLogPosterior += dPriorIndependent[i];
				dLogPosterior += oIndependentModule->getIndependentPosterior(id, i);
			}
		}
	}
	iNumIndependentVariants = viIndependentVariants.size();
	vector<UINT32> viVariants;
	for(UINT32 i = 1 + iNumHyperGroupTypes;i < iNumHyperGroupTypes * 2;i++) {
		mcpDependentAssociation[i] = new DependentModule(iMaxVariants, gwasData,
				vvviHyperGroup.at(i - iNumHyperGroupTypes), i);
		viVariants.clear();
		for(UINT32 j = 0;j < iPriorNumDependentVariants;) {
			id = (UINT32) ((((double) rand()) / ((double) RAND_MAX)) * ((double) iNumVariants));
			if(acAssociationTypesForAllVariants[id] == 0) {
				acAssociationTypesForAllVariants[id] = i;
				j++;
				dLogPosterior += dPriorDependent[i];
				viVariants.push_back(id);
			}
		}
		mcpDependentAssociation[i]->addVariants(viVariants);
		dLogPosterior = mcpDependentAssociation[i]->calPosterior();
	}
	return true;
}

bool MCMC::isConveraged(double _dLogPosterior) {
	UINT32 i, len = iLogPosteriorSize / 2;
	if(vdLogPosterior.size() < iLogPosteriorSize) {
		vdLogPosterior.push_back(_dLogPosterior);
		iLogPosteriorIndex += 1;
	} else {
		iLogPosteriorIndex = iLogPosteriorIndex % iLogPosteriorSize;
		vdLogPosterior.at(iLogPosteriorIndex) = _dLogPosterior;
		double sum1 = 0, sum2 = 0, s1 = 0, s2 = 0;
		for(i = 1;i <= len;++i) {
			sum1 += vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize);
		}
		for(;i <= iLogPosteriorSize;++i) {
			sum2 += vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize);
		}
		sum1 /= ((double) len);
		sum2 /= ((double) len);
		for(i = 1;i <= len;++i) {
			s1 += (vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize) - sum1)
					* (vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize) - sum1);
		}
		for(;i <= iLogPosteriorSize;++i) {
			s2 += (vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize) - sum2)
					* (vdLogPosterior.at((i + iLogPosteriorIndex) % iLogPosteriorSize) - sum2);
		}
		if(s2 / s1 > 0.9 && s2 / s1 < 1.1) {
			return true;
		}
	}
	return false;
}

void MCMC::recordFrequency() {
	UINT32 iUnit = iNumHyperGroupTypes * 2;
	for(UINT32 i = 0;i < iNumIndependentVariants;i++) {
		iFrequencyForAllVairants[viIndependentVariants.at(i) * iUnit
				+ acAssociationTypesForAllVariants[viIndependentVariants.at(i)]] += 1;
	}
	map<char, DependentModule*>::iterator it;
	vector<UINT32> * pvi;
	for(it = mcpDependentAssociation.begin();it != mcpDependentAssociation.end();it++) {
		pDependentModule = it->second;
		pvi = pDependentModule->getVariants();
		for(UINT32 i = 0;i < pvi->size();i++) {
			iFrequencyForAllVairants[pvi->at(i) * iUnit + it->first] += 1;
		}
	}
}

bool MCMC::runMCMC() {
	iNumLimitDependentVariants = (UINT32) (log(iNumTotalSamples) / log(iVariantTypes) - 1.0);
	iNumTrueLimit = (UINT32) (log(((double) iNumTotalSamples) / 10.0) / log(iVariantTypes) - 1);
	if(iNumLimitDependentVariants <= 0) {
		iNumLimitDependentVariants = 1;
	}
	if(iNumTrueLimit <= 1) {
		iNumTrueLimit = 2;
	}
	UINT32 iMinDependentVariants = 0;
	double T = 1;
	double sp0 = 0.3;	//sp: frequency for different moves
	double sp1 = 0.3;
	double sp2 = 0.1;
	double p0 = 1. / 10.; //p: frequency for updating interaction group between other groups
	double p1 = 8. / 10.;
	double p2 = 1. - p0 - p1;
	double dRandomNumber = 0;

	UINT64 iTotalIteration = iNumBurnins + iNumIterations;
	for(UINT64 i = 0;i < iTotalIteration;i++) {
		if(i % 1000ULL == 0) {
			cout << ".";
		}
		iMinDependentVariants = 0;
		sp0 = 0.7;	//sp: frequency for different moves
		sp1 = 0.15;
		sp2 = 0.1;
		p0 = 1. / 10.; //p: frequency for updating interaction group between other groups
		p1 = 8. / 10.;
		p2 = 1. - p0 - p1;
		if(i < iNumBurnins) {
			p0 = 1. / 3. * (double) (i + 1) / (double) iNumBurnins;
			p1 = 1. - p0 * 2.;
			p2 = 1. - p0 - p1;
			iMinDependentVariants = (int) ((double) (iNumLimitDependentVariants - 1)
					* (1. - exp(-(double) (iNumBurnins - i) / iNumBurnins)) / (1. - exp(-1.))) + 1;
		}
		T = 1;
		if(i < iNumBurnins / 1.5) {
			T += log(iNumTotalSamples * 3) * ((double) iNumBurnins / 1.5 - i) * 1.5 / (double) iNumBurnins;
		}

		dRandomNumber = ((double) rand()) / ((double) RAND_MAX);

		if(dRandomNumber < sp0) {
			updateDependentNoise(p0, p1, p2, T, iMinDependentVariants);
		} else if(dRandomNumber < sp0 + sp1) {
			updateIndependentNoise(p0, p1, p2, T);
		} else if(dRandomNumber < sp0 + sp1 + sp2) {
			updateDependentIndependent(p0, p1, p2, T, iMinDependentVariants);
		} else {
			dRandomNumber = ((double) rand()) / ((double) RAND_MAX);
			if(dRandomNumber < 0.8) {
				updateDependentIntra(T);
			} else if(dRandomNumber < 1.0) {
				updateIndependentIntra(T);
			}
		}

		if(i >= iNumBurnins) {
			if((i - iNumBurnins) % STEP == 0) {
				shiftPosition();
				recordFrequency();
				if(isConveraged(dLogPosterior)) {
					return true;
				}
			}
		}
	}
	return true;
}

/**
 * setup all the static variables in MCMC
 */
bool MCMC::setupMCMC(GwasData * _gwasData) {
	iNumVariants = Config::iNumVariants;
	iNumGroups = Config::iNumGroups;
	iNumBurnins = iNumVariants * 400;
	iNumIterations = ((UINT64) iNumVariants) * (iNumVariants > 100000 ? 100000 : ((UINT64) iNumVariants));
	iNumTotalSamples = _gwasData->iNumTotalSamplesAcrossGroups;
	iVariantTypes = _gwasData->getVariantNumTypes();
	gwasData = _gwasData;
	iPriorNumIndependentVariants = Config::iNumIndependentVariants;
	iPriorNumDependentVariants = Config::iNumDependentVariants;

	vvviHyperGroup = HyperGroup::getHyperGroup(iNumGroups);
	iNumHyperGroupTypes = vvviHyperGroup.size();

	dPriorIndependent = new double[iNumHyperGroupTypes];
	dPriorDependent = new double[iNumHyperGroupTypes * 2];

	for(UINT32 i = 1;i < iNumHyperGroupTypes;i++) {
		dPriorIndependent[i] = iPriorNumIndependentVariants / ((double) iNumVariants);
		dPriorIndependent[0] += dPriorIndependent[i];
		dPriorDependent[i + iNumHyperGroupTypes] = iPriorNumDependentVariants / ((double) iNumVariants);
		dPriorDependent[0 + iNumHyperGroupTypes] += dPriorDependent[i];
	}
	dPriorIndependent[0] = 1 - dPriorIndependent[0];

	dPriorDependent[0 + iNumHyperGroupTypes] = 1 - dPriorDependent[0 + iNumHyperGroupTypes];
	for(UINT32 i = 0;i < iNumHyperGroupTypes;i++) {
		dPriorIndependent[i] = log(dPriorIndependent[i]);
		dPriorDependent[i + iNumHyperGroupTypes] = log(dPriorDependent[i + iNumHyperGroupTypes]);
	}

	// Initialize the posterior probability for each variant
	oIndependentModule = new IndependentModule();
	oIndependentModule->initial(gwasData, vvviHyperGroup);
	iMaxVariants = Config::iMaxVariants;

	return true;
}

/**
 * find the position with the largest probability
 */
void MCMC::shiftPosition() {
	map<char, DependentModule*>::iterator it;
	vector<UINT32> * viVariants;
	UINT32 iChromsomeIndex1;
	UINT32 i, j, k;
	int iPos1, iPos2;
	bool bTooClose = false;
	vector<UINT32> viTemp, viTemp2;
	vector<double> vdPosterior;
	double maxPosterior;
	UINT32 iVariantId, iOldVariantId;
	for(it = mcpDependentAssociation.begin();it != mcpDependentAssociation.end();++it) {
		if(it->second->iNumVariants == 0) {
			continue;
		}
		viVariants = it->second->getVariants();
		for(i = 0;i < viVariants->size();++i) {
			iOldVariantId = viVariants->at(i);
			gwasData->getVariantNearby(iOldVariantId, iMiniDistance, viNearbyVariants);
			// remove variants that are too close to the existing one
			viTemp.clear();
			for(j = 0;j < viNearbyVariants.size();++j) {
				bTooClose = false;
				iPos1 = gwasData->getVariantPosition(viNearbyVariants.at(j));
				iChromsomeIndex1 = gwasData->getVariantChromsome(viNearbyVariants.at(j));
				for(k = 0;k < viVariants->size();++k) {
					if(k != i && iChromsomeIndex1 == gwasData->getVariantChromsome(viVariants->at(k))) {
						iPos2 = gwasData->getVariantPosition(viVariants->at(k));
						if(abs(iPos1 - iPos2) < iMiniDistance) {
							bTooClose = true;
							break;
						}
					}
				}
				if(!bTooClose) {
					viTemp.push_back(viNearbyVariants.at(j));
				}
			}
			if(viTemp.empty()) {
				continue;
			}
			viNearbyVariants = viTemp;
			// try all the nearby variants
			vdPosterior.clear();
			viTemp.clear();
			viTemp.push_back(i);
			for(j = 0;j < viNearbyVariants.size();j++) {
				viTemp2.clear();
				viTemp2.push_back(viNearbyVariants.at(j));
				it->second->replaceVariants(viTemp, viTemp2);
				vdPosterior.push_back(it->second->calPosterior());
			}
			iVariantId = 0;
			maxPosterior = vdPosterior.at(0);
			for(j = 1;j < vdPosterior.size();++j) {
				if(vdPosterior.at(j) > maxPosterior) {
					maxPosterior = vdPosterior.at(j);
					iVariantId = j;
				}
			}
			// apply if the probability is smaller
			if(maxPosterior > it->second->getPosterior()) {
				dLogPosterior += (maxPosterior - it->second->getPosterior());
				viTemp.clear();
				viTemp.push_back(i);
				viTemp2.clear();
				viTemp2.push_back(viNearbyVariants.at(iVariantId));
				it->second->replaceVariants(viTemp, viTemp2);
				it->second->calPosterior();
				it->second->apply();
				acAssociationTypesForAllVariants[iOldVariantId] = acAssociationTypesForAllVariants[viNearbyVariants.at(
						iVariantId)];
				acAssociationTypesForAllVariants[viNearbyVariants.at(iVariantId)] = it->first;
			}
		}
	}
}

bool MCMC::test() {
	DependentModule * p = new DependentModule(20, gwasData, vvviHyperGroup.at(1), 1);
	vector<UINT32> viTemp;
	viTemp.push_back(101);
	viTemp.push_back(110);
	viTemp.push_back(210);
	viTemp.push_back(310);
	viTemp.push_back(410);
	viTemp.push_back(510);
	viTemp.push_back(710);
	viTemp.push_back(810);
	viTemp.push_back(910);
	double begin = omp_get_wtime();
	for(UINT32 i = 0;i < 50;i++) {
		p->addVariants(viTemp);
		p->cleanHashtable();
	}
	viTemp.clear();
	viTemp.push_back(10);
	viTemp.push_back(110);
	viTemp.push_back(783);
	p->addVariants(viTemp);
	p->pHashTable->print();
	p->calPosterior();
	p->cleanHashtable();
	double end = omp_get_wtime();
	cout << "finished in " << double(end - begin) << " Seconds." << endl << endl;
	delete p;
	return true;
}

/**
 * C(n, m)
 */
double combination(double _n, double _m) {
	double dResult = 1;
	for(double i = 0;i < _m;i++) {
		dResult *= (_n - i) / (i + 1);
	}
	return dResult;
}

bool MCMC::updateDependentIndependent(double _p0, double _p1, double _p2, double _T, UINT32 _iMinDependentVariants) {
	iIndexTypeDependentAssoication = get_random_number((iNumHyperGroupTypes - 1)) + (iNumHyperGroupTypes + 1);
	UINT32 iNum = 0;
	pDependentModule = mcpDependentAssociation[iIndexTypeDependentAssoication];
	iNum = pDependentModule->iNumVariants;
	if(iNum == 0) { // Independent -> Dependent
		dProbabilityMove = _p0 + _p1 + (((double) rand()) / ((double) RAND_MAX)) * _p2;
	} else if(iNum < _iMinDependentVariants) { // Independent -> Dependent
		dProbabilityMove = _p0 + _p1 + (((double) rand()) / ((double) RAND_MAX)) * _p2;
	} else if(iNum > iNumLimitDependentVariants) { // Dependent -> Independent
		dProbabilityMove = (((double) rand()) / ((double) RAND_MAX)) * _p0;
	} else if(_iMinDependentVariants == iNumLimitDependentVariants) { // Dependent <-> Independent
		dProbabilityMove = _p0 + (((double) rand()) / ((double) RAND_MAX)) * _p1;
	} else if(iNum == _iMinDependentVariants) { // Dependent <-> Independent, Independent -> Dependent
		dProbabilityMove = _p0 + (((double) rand()) / ((double) RAND_MAX)) * (_p1 + _p2);
	} else if(iNum == iNumLimitDependentVariants) { // Dependent <-> Independent, Independent -> Dependent
		dProbabilityMove = ((double) rand()) / ((double) RAND_MAX) * (_p0 + _p1);
	} else {
		dProbabilityMove = ((double) rand()) / ((double) RAND_MAX);
	}
	// Dependent -> Independent
	if(dProbabilityMove < _p0) {
		do {
			iNumSelectVariants = (get_random_number(iNum) + 1);
			iNumRestVariants = iNum - iNumSelectVariants;
		} while(iNumRestVariants < iNumLimitDependentVariants);

		vcAssociationTypes.clear();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			vcAssociationTypes.push_back(get_random_number(iNumHyperGroupTypes - 1) + 1);
		}

		if(iNum < iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log((_p2 * _p0) / (_p1 + _p2)); // (_p2 / (_p1 + _p2)) / _p0
			} else {
				logratio = log(_p2 / _p0);
			}
		} else if(iNum == iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log((_p2 * (_p0 + _p1)) / (_p0 * (_p1 + _p2))); // (_p2/(_p1+_p2)) / (_p0/(_p0+_p1))
			} else {
				logratio = log(_p2 / _p0);
			}
		} else if(iNum > iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = 0;
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = 0;
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log(_p2 / (_p1 + _p2)); // (_p2/(_p1+_p2)) / 1.0)
			} else {
				logratio = log(_p2); // _p2 / 1.0
			}
		}

		logratio += log(
				(combination(iNum, iNumSelectVariants))
						/ combination(iNumIndependentVariants + iNumSelectVariants, iNumSelectVariants));

		viInnerSelectedVariantIds.clear();
		viOuterSelectedVariantIds.clear();
		shuffleSelect(iNum, iNumSelectVariants, viInnerSelectedVariantIds);
		pDependentModule->getVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIds);
		pDependentModule->delVariants(viInnerSelectedVariantIds);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i),
					vcAssociationTypes.at(i));
			dPosteriorDifference += (dPriorIndependent[(UINT32) vcAssociationTypes.at(i)]
					- dPriorDependent[iIndexTypeDependentAssoication]);
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = vcAssociationTypes.at(i);
				viIndependentVariants.push_back(viOuterSelectedVariantIds.at(i));
			}
			iNumDependentVariants -= iNumSelectVariants;
			iNumIndependentVariants += iNumSelectVariants;
		} else {
			pDependentModule->rollBack();
		}
		// Dependent <-> Independent
	} else if(dProbabilityMove < _p0 + _p1) {
		if(iNumIndependentVariants == 0) {
			return true;
		}
		iNumSelectVariants = (get_random_number(iNum) + 1);
		viInnerSelectedVariantIds.clear();
		viOuterSelectedVariantIds.clear();
		shuffleSelect(iNum, iNumSelectVariants, viInnerSelectedVariantIds);
		pDependentModule->getVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIds);
		pDependentModule->getRestVariants(viInnerSelectedVariantIds, viOuterUnselectedVariantIds);
		viOuterSelectedVariantIdsSwitch.clear();
		shuffleSelect(viIndependentVariants, iNumSelectVariants, viOuterSelectedVariantIdsSwitch);
		dPosteriorDifference = 0;
		pDependentModule->replaceVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIdsSwitch);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += (oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i),
					acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)])
					- oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIdsSwitch.at(i),
							acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)]));
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T >= 0 || dPosteriorDifference / _T >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			viIndependentVariants.erase(viIndependentVariants.begin() + (iNumIndependentVariants - iNumSelectVariants),
					viIndependentVariants.end());
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = viOuterSelectedVariantIdsSwitch.at(
						i);
				acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] =
						iIndexTypeDependentAssoication;
				viIndependentVariants.push_back(viOuterSelectedVariantIds.at(i));
			}
		} else {
			pDependentModule->rollBack();
		}
		// Independent -> Dependent
	} else {
		if(iNumIndependentVariants == 0) {
			return true;
		}
		if(iNumLimitDependentVariants + 1 - iNum > iNumIndependentVariants) {
			iNumSelectVariants = (get_random_number(iNumIndependentVariants) + 1);
		} else {
			iNumSelectVariants = (get_random_number(iNumLimitDependentVariants + 1 - iNum) + 1);
		}

		iNumNewTotalVariants = iNum + iNumSelectVariants;
		if(iNum == 0) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0);
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 + _p1) / _p0);
			} else {
				logratio = 0;
			}
		} else if(iNum < _iMinDependentVariants) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0);
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 + _p1) / _p0);
			} else {
				logratio = 0;
			}
		} else if(iNum == _iMinDependentVariants) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log((_p0 * (_p1 + _p2)) / _p2); // _p0 / (_p2/(_p1+_p2))
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 * (_p1 + _p2)) / (_p2 * (_p0 + _p1))); // (_p0/(_p0+_p1)) / (_p2/(_p1+_p2))
			} else {
				logratio = log((_p1 + _p2) / _p2);
			}
		} else {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0 / _p2); // _p0 / (_p2/(_p1+_p2))
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log(_p0 / ((_p0 + _p1) * _p2)); // (_p0/(_p0+_p1)) / (_p2/(_p1+_p2))
			} else {
				logratio = log(1.0 / _p2);
			}
		}
		logratio += log(
				(combination(iNumIndependentVariants, iNumSelectVariants))
						/ (combination(iNumNewTotalVariants, iNumSelectVariants)));

		shuffleSelect(viIndependentVariants, iNumSelectVariants, viOuterSelectedVariantIds);

		pDependentModule->addVariants(viOuterSelectedVariantIds);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += (-oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i),
					acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)]));
			dPosteriorDifference += (dPriorDependent[iIndexTypeDependentAssoication]
					- dPriorIndependent[(UINT32) acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)]]);
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = iIndexTypeDependentAssoication;
			}
			viIndependentVariants.erase(viIndependentVariants.begin() + (iNumIndependentVariants - iNumSelectVariants),
					viIndependentVariants.end());
			iNumDependentVariants += iNumSelectVariants;
			iNumIndependentVariants -= iNumSelectVariants;
		} else {
			pDependentModule->rollBack();
		}
	}
	return true;
}

bool MCMC::updateDependentIntra(double _T) {
	iIndexTypeDependentAssoication = get_random_number((iNumHyperGroupTypes - 1)) + (iNumHyperGroupTypes + 1);
	do {
		iIndexTypeDependentAssoicationSwitch = get_random_number((iNumHyperGroupTypes - 1)) + (iNumHyperGroupTypes + 1);
	} while(iIndexTypeDependentAssoication == iIndexTypeDependentAssoicationSwitch);
	pDependentModule = mcpDependentAssociation[iIndexTypeDependentAssoication];
	pDependentModuleSwitch = mcpDependentAssociation[iIndexTypeDependentAssoicationSwitch];
	UINT32 iNum = 0, iNumSwitch = 0;
	iNum = pDependentModule->iNumVariants;
	iNumSwitch = pDependentModuleSwitch->iNumVariants;
	if(iNum == 0 && iNumSwitch == 0) {
		return true;
	}
	if(iNum != 0) {
		iNumSelectVariants = (get_random_number(iNum) + 1);
		shuffleSelect(iNum, iNumSelectVariants, viInnerSelectedVariantIds);
		pDependentModule->getVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIds);
	}
	if(iNumSwitch != 0) {
		vector<UINT32> * pviVariants = pDependentModuleSwitch->getVariants();
		viShuffleSelect.clear();
		for(UINT32 i = 0;i < pviVariants->size();i++) {
			if(!isTooClose(pviVariants->at(i), viOuterSelectedVariantIds)) {
				viShuffleSelect.push_back(i);
			}
		}
		iNumSelectVariantsSwitch = (get_random_number(viShuffleSelect.size()) + 1);
		shuffleSelect(viShuffleSelect, iNumSelectVariantsSwitch, viInnerSelectedVariantIdsSwitch);
		pDependentModuleSwitch->getVariants(viInnerSelectedVariantIdsSwitch, viOuterSelectedVariantIdsSwitch);
	}
	if(iNumSelectVariants == iNumSelectVariantsSwitch) {
		pDependentModule->replaceVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIdsSwitch);
		pDependentModuleSwitch->replaceVariants(viInnerSelectedVariantIdsSwitch, viOuterSelectedVariantIds);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		dPosteriorDifference += pDependentModuleSwitch->calPosterior() - pDependentModuleSwitch->getPosterior();
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T >= 0 || dPosteriorDifference / _T >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			pDependentModuleSwitch->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] =
						iIndexTypeDependentAssoicationSwitch;
				acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] =
						iIndexTypeDependentAssoication;
			}
		} else {
			pDependentModule->rollBack();
			pDependentModuleSwitch->rollBack();
		}
	} else {
		pDependentModule->delVariants(viInnerSelectedVariantIds);
		pDependentModule->addVariants(viOuterSelectedVariantIdsSwitch);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		pDependentModuleSwitch->delVariants(viInnerSelectedVariantIdsSwitch);
		pDependentModuleSwitch->addVariants(viOuterSelectedVariantIds);
		dPosteriorDifference += pDependentModuleSwitch->calPosterior() - pDependentModuleSwitch->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += (dPriorDependent[iIndexTypeDependentAssoicationSwitch]
					- dPriorDependent[iIndexTypeDependentAssoication]);
		}
		for(UINT32 i = 0;i < iNumSelectVariantsSwitch;i++) {
			dPosteriorDifference += (dPriorDependent[iIndexTypeDependentAssoication]
					- dPriorDependent[iIndexTypeDependentAssoicationSwitch]);
		}
		if(dPosteriorDifference / _T >= 0 || dPosteriorDifference / _T >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			pDependentModuleSwitch->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] =
						iIndexTypeDependentAssoicationSwitch;
				acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] =
						iIndexTypeDependentAssoication;
			}
			for(UINT32 i = 0;i < iNumSelectVariantsSwitch;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] =
						iIndexTypeDependentAssoication;
			}
		} else {
			pDependentModule->rollBack();
			pDependentModuleSwitch->rollBack();
		}
	}
	return true;
}

bool MCMC::updateDependentNoise(double _p0, double _p1, double _p2, double _T, UINT32 _iMinDependentVariants) {
	iIndexTypeDependentAssoication = get_random_number((iNumHyperGroupTypes - 1)) + (iNumHyperGroupTypes + 1);
	UINT32 iNum = 0;
	pDependentModule = mcpDependentAssociation[iIndexTypeDependentAssoication];
	iNum = pDependentModule->iNumVariants;
	if(iNum == 0) { // noise -> dependent
		dProbabilityMove = _p0 + _p1 + (((double) rand()) / ((double) RAND_MAX)) * _p2;
	} else if(iNum < _iMinDependentVariants) { // noise -> dependent
		dProbabilityMove = _p0 + _p1 + (((double) rand()) / ((double) RAND_MAX)) * _p2;
	} else if(iNum > iNumLimitDependentVariants) { // dependent -> noise
		dProbabilityMove = (((double) rand()) / ((double) RAND_MAX)) * _p0;
	} else if(_iMinDependentVariants == iNumLimitDependentVariants) { // dependent <-> noise
		dProbabilityMove = _p0 + (((double) rand()) / ((double) RAND_MAX)) * _p1;
	} else if(iNum == _iMinDependentVariants) { // dependent <-> noise, noise -> dependent
		dProbabilityMove = _p0 + (((double) rand()) / ((double) RAND_MAX)) * (_p1 + _p2);
	} else if(iNum == iNumLimitDependentVariants) { // dependent <-> noise, noise -> dependent
		dProbabilityMove = ((double) rand()) / ((double) RAND_MAX) * (_p0 + _p1);
	} else {
		dProbabilityMove = ((double) rand()) / ((double) RAND_MAX);
	}
	// Dependent -> noise
	if(dProbabilityMove < _p0) {
		do {
			iNumSelectVariants = (get_random_number(iNum) + 1);
			iNumRestVariants = iNum - iNumSelectVariants;
		} while(iNumRestVariants < iNumLimitDependentVariants);

		if(iNum < iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log((_p2 * _p0) / (_p1 + _p2)); // (_p2 / (_p1 + _p2)) / _p0
			} else {
				logratio = log(_p2 / _p0);
			}
		} else if(iNum == iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = log(1. / _p0);
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log((_p2 * (_p0 + _p1)) / (_p0 * (_p1 + _p2))); // (_p2/(_p1+_p2)) / (_p0/(_p0+_p1))
			} else {
				logratio = log(_p2 / _p0);
			}
		} else if(iNum > iNumLimitDependentVariants) {
			if(iNumRestVariants == 0) {
				logratio = 0;
			} else if(iNumRestVariants < _iMinDependentVariants) {
				logratio = 0;
			} else if(iNumRestVariants == _iMinDependentVariants) {
				logratio = log(_p2 / (_p1 + _p2)); // (_p2/(_p1+_p2)) / 1.0)
			} else {
				logratio = log(_p2); // _p2 / 1.0
			}
		}

		logratio += log(
				(combination(iNum, iNumSelectVariants))
						/ combination(
								iNumVariants - iNumDependentVariants - iNumIndependentVariants + iNumSelectVariants,
								iNumSelectVariants));

		viInnerSelectedVariantIds.clear();
		viOuterSelectedVariantIds.clear();
		shuffleSelect(iNum, iNumSelectVariants, viInnerSelectedVariantIds);
		pDependentModule->getVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIds);
		pDependentModule->delVariants(viInnerSelectedVariantIds);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i), 0);
			dPosteriorDifference += (dPriorIndependent[0] - dPriorDependent[iIndexTypeDependentAssoication]);
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = 0;
			}
			iNumDependentVariants -= iNumSelectVariants;
		} else {
			pDependentModule->rollBack();
		}
		// Dependent <-> noise
	} else if(dProbabilityMove < _p0 + _p1) {
		iNumSelectVariants = (get_random_number(iNum) + 1);
		viInnerSelectedVariantIds.clear();
		viOuterSelectedVariantIds.clear();
		shuffleSelect(iNum, iNumSelectVariants, viInnerSelectedVariantIds);
		pDependentModule->getVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIds);
		pDependentModule->getRestVariants(viInnerSelectedVariantIds, viOuterUnselectedVariantIds);
		UINT32 iSelectInd = 0;
		viOuterSelectedVariantIdsSwitch.clear();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			do {
				iSelectInd = get_random_number(iNumVariants);
			} while(acAssociationTypesForAllVariants[iSelectInd] != 0
					|| abSelectedStatusForAllVariants[iSelectInd] == true
					|| isTooClose(iSelectInd, viOuterUnselectedVariantIds)
					|| isTooClose(iSelectInd, viOuterSelectedVariantIdsSwitch));
			viOuterSelectedVariantIdsSwitch.push_back(iSelectInd);
			abSelectedStatusForAllVariants[iSelectInd] = true;
		}
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			abSelectedStatusForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] = false;
		}
		dPosteriorDifference = 0;
		pDependentModule->replaceVariants(viInnerSelectedVariantIds, viOuterSelectedVariantIdsSwitch);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += (oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i), 0)
					- oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIdsSwitch.at(i), 0));
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = 0;
				acAssociationTypesForAllVariants[viOuterSelectedVariantIdsSwitch.at(i)] =
						iIndexTypeDependentAssoication;
			}
		} else {
			pDependentModule->rollBack();
		}
		// Dependent <- noise
	} else {
		iNumSelectVariants = (get_random_number(iNumLimitDependentVariants + 1 - iNum) + 1);
		iNumNewTotalVariants = iNum + iNumSelectVariants;
		if(iNum == 0) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0);
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 + _p1) / _p0);
			} else {
				logratio = 0;
			}
		} else if(iNum < _iMinDependentVariants) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0);
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 + _p1) / _p0);
			} else {
				logratio = 0;
			}
		} else if(iNum == _iMinDependentVariants) {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log((_p0 * (_p1 + _p2)) / _p2); // _p0 / (_p2/(_p1+_p2))
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log((_p0 * (_p1 + _p2)) / (_p2 * (_p0 + _p1))); // (_p0/(_p0+_p1)) / (_p2/(_p1+_p2))
			} else {
				logratio = log((_p1 + _p2) / _p2);
			}
		} else {
			if(iNumNewTotalVariants < iNumLimitDependentVariants) {
				logratio = log(_p0 / _p2); // _p0 / (_p2/(_p1+_p2))
			} else if(iNumNewTotalVariants == iNumLimitDependentVariants) {
				logratio = log(_p0 / ((_p0 + _p1) * _p2)); // (_p0/(_p0+_p1)) / (_p2/(_p1+_p2))
			} else {
				logratio = log(1.0 / _p2);
			}
		}
		logratio += log(
				(combination(iNumVariants - iNumDependentVariants - iNumIndependentVariants, iNumSelectVariants))
						/ (combination(iNumNewTotalVariants, iNumSelectVariants)));
		UINT32 iSelectInd = 0;
		viOuterSelectedVariantIds.clear();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			do {
				iSelectInd = get_random_number(iNumVariants);
			} while(acAssociationTypesForAllVariants[iSelectInd] != 0
					|| abSelectedStatusForAllVariants[iSelectInd] == true
					|| isTooClose(iSelectInd, viOuterSelectedVariantIds)
					|| isTooClose(iSelectInd, (*(pDependentModule->getVariants()))));
			viOuterSelectedVariantIds.push_back(iSelectInd);
			abSelectedStatusForAllVariants[iSelectInd] = true;
		}
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			abSelectedStatusForAllVariants[viOuterSelectedVariantIds.at(iSelectInd)] = false;
		}
		pDependentModule->addVariants(viOuterSelectedVariantIds);
		dPosteriorDifference = pDependentModule->calPosterior() - pDependentModule->getPosterior();
		for(UINT32 i = 0;i < iNumSelectVariants;i++) {
			dPosteriorDifference += (-oIndependentModule->getIndependentPosterior(viOuterSelectedVariantIds.at(i), 0));
			dPosteriorDifference += (dPriorDependent[iIndexTypeDependentAssoication] - dPriorIndependent[0]);
		}
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			pDependentModule->apply();
			dLogPosterior += dPosteriorDifference;
			for(UINT32 i = 0;i < iNumSelectVariants;i++) {
				acAssociationTypesForAllVariants[viOuterSelectedVariantIds.at(i)] = iIndexTypeDependentAssoication;
			}
			iNumDependentVariants += iNumSelectVariants;
		} else {
			pDependentModule->rollBack();
		}
	}
	return true;
}

bool MCMC::updateIndependentIntra(double _T) {
	if(iNumIndependentVariants == 0) {
		return true;
	}
	iSelectedVariantInner = get_random_number(iNumIndependentVariants);
	iSelectedVariantOuter = viIndependentVariants.at(iSelectedVariantInner);
	do {
		iIndexTypeDependentAssoication = get_random_number(iNumHyperGroupTypes) + 1;
	} while(iIndexTypeDependentAssoication == (UINT32) acAssociationTypesForAllVariants[iSelectedVariantOuter]);
	dPosteriorDifference = oIndependentModule->getIndependentPosterior(iSelectedVariantOuter,
			iIndexTypeDependentAssoication)
			- oIndependentModule->getIndependentPosterior(iSelectedVariantOuter,
					acAssociationTypesForAllVariants[iSelectedVariantOuter]);

	dProbabilityRandom = get_random_number();
	if(dPosteriorDifference / _T >= 0 || dPosteriorDifference / _T >= log(dProbabilityRandom)) {
		dLogPosterior += dPosteriorDifference;
		acAssociationTypesForAllVariants[iSelectedVariantOuter] = iIndexTypeDependentAssoication;
	}
	return true;
}

bool MCMC::updateIndependentNoise(double _p0, double _p1, double _p2, double _T) {
	if(iNumIndependentVariants == 0) {
		dProbabilityMove = _p0 + _p1 + get_random_number() * _p2;
	} else {
		dProbabilityMove = get_random_number();
	}
	// Independent -> Noise
	if(dProbabilityMove < _p0) {
		if(iNumIndependentVariants == 1) {
			logratio = log(1.0 / _p0);
		} else {
			logratio = log(_p2 / _p0);
		}

		logratio += log(
				((double) iNumIndependentVariants)
						/ (((double) (iNumHyperGroupTypes - 1)
								* (iNumVariants + 1 - iNumDependentVariants - iNumIndependentVariants))));
		iSelectedVariantInner = get_random_number(iNumIndependentVariants);
		iSelectedVariantOuter = viIndependentVariants.at(iSelectedVariantInner);
		dPosteriorDifference = oIndependentModule->getIndependentPosterior(iSelectedVariantOuter, 0)
				- oIndependentModule->getIndependentPosterior(iSelectedVariantOuter,
						acAssociationTypesForAllVariants[iSelectedVariantOuter]);
		dPosteriorDifference += (dPriorIndependent[0]
				- dPriorIndependent[(UINT32) acAssociationTypesForAllVariants[iSelectedVariantOuter]]);
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			dLogPosterior += dPosteriorDifference;
			viIndependentVariants.erase(viIndependentVariants.begin() + iSelectedVariantInner);
			acAssociationTypesForAllVariants[iSelectedVariantOuter] = 0;
			iNumIndependentVariants -= 1;
		}
		// Independent <-> Noise
	} else if(dProbabilityMove < _p0 + _p1) {
		iSelectedVariantInner = get_random_number(iNumIndependentVariants);
		iSelectedVariantOuter = viIndependentVariants.at(iSelectedVariantInner);
		do {
			iSelectedVariantOuterSwitch = get_random_number(iNumVariants);
		} while(acAssociationTypesForAllVariants[iSelectedVariantOuterSwitch] != 0);
		dPosteriorDifference = oIndependentModule->getIndependentPosterior(iSelectedVariantOuter, 0)
				- oIndependentModule->getIndependentPosterior(iSelectedVariantOuter,
						acAssociationTypesForAllVariants[iSelectedVariantOuter]);
		dPosteriorDifference += (oIndependentModule->getIndependentPosterior(iSelectedVariantOuterSwitch,
				acAssociationTypesForAllVariants[iSelectedVariantOuter])
				- oIndependentModule->getIndependentPosterior(iSelectedVariantOuterSwitch, 0));
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T >= 0 || dPosteriorDifference / _T >= log(dProbabilityRandom)) {
			dLogPosterior += dPosteriorDifference;
			viIndependentVariants.at(iSelectedVariantInner) = iSelectedVariantOuterSwitch;
			acAssociationTypesForAllVariants[iSelectedVariantOuterSwitch] =
					acAssociationTypesForAllVariants[iSelectedVariantOuter];
			acAssociationTypesForAllVariants[iSelectedVariantOuter] = 0;
		}
		// Independent <- Noise
	} else {
		if(iNumIndependentVariants == 0) {
			logratio = log(_p0);
		} else {
			logratio = log(_p0 / _p2);
		}
		logratio +=
				log(
						(((double) (iNumHyperGroupTypes - 1)
								* (iNumVariants - iNumDependentVariants - iNumIndependentVariants)))
								/ (((double) iNumIndependentVariants + 1)));
		do {
			iSelectedVariantOuter = get_random_number(iNumVariants);
		} while(acAssociationTypesForAllVariants[iSelectedVariantOuter] != 0);
		iIndexTypeDependentAssoication = get_random_number(iNumHyperGroupTypes) + 1;
		dPosteriorDifference = oIndependentModule->getIndependentPosterior(iSelectedVariantOuter,
				iIndexTypeDependentAssoication) - oIndependentModule->getIndependentPosterior(iSelectedVariantOuter, 0);
		dPosteriorDifference += (dPriorIndependent[iIndexTypeDependentAssoication] - dPriorIndependent[0]);
		dProbabilityRandom = get_random_number();
		if(dPosteriorDifference / _T + logratio >= 0
				|| dPosteriorDifference / _T + logratio >= log(dProbabilityRandom)) {
			dLogPosterior += dPosteriorDifference;
			viIndependentVariants.push_back(iSelectedVariantOuter);
			acAssociationTypesForAllVariants[iSelectedVariantOuter] = iIndexTypeDependentAssoication;
			iNumIndependentVariants += 1;
		}
	}
	return true;
}

// private method

/**
 * check if Outer id is close to the given list of variants
 */
bool MCMC::isTooClose(UINT32 _iOuterId, vector<UINT32> & _viSelectedVariants) {
	UINT32 iChromsomeIndex = gwasData->getVariantChromsome(_iOuterId);
	int iPos = gwasData->getVariantPosition(_iOuterId);
	for(UINT32 i = 0;i < _viSelectedVariants.size();i++) {
		if(iChromsomeIndex == gwasData->getVariantChromsome(_viSelectedVariants.at(i))) {
			if(abs(iPos - gwasData->getVariantPosition(_viSelectedVariants.at(i))) < iMiniDistance) {
				return true;
			}
		}
	}
	return false;
}

void MCMC::get_close_variant(UINT32 _iOuterId, vector<UINT32> & _viExistVariants,
		vector<UINT32> & _viCandidateVariants) {

}

double MCMC::get_random_number() {
	return distribution(generator);
}

UINT32 MCMC::get_random_number(UINT32 _iMax) {
	dRandomNumber = distribution(generator);
	return (UINT32) dRandomNumber * _iMax;
}

void MCMC::shuffleSelect(vector<UINT32> & _viShuffleSelect, UINT32 _iNumSelect,
		vector<UINT32> & _viInnerSelectedVariantIds) {
	UINT32 iSelectInd = 0;
	_viInnerSelectedVariantIds.clear();
	for(UINT32 i = 0;i < _iNumSelect;i++) {
		iSelectInd = get_random_number(_viShuffleSelect.size() - i);
		_viInnerSelectedVariantIds.push_back(_viShuffleSelect.at(iSelectInd));
		iter_swap(_viShuffleSelect.begin() + iSelectInd, _viShuffleSelect.begin() + (_viShuffleSelect.size() - i - 1));
	}
}

void MCMC::shuffleSelect(UINT32 _iNum, UINT32 _iNumSelect, vector<UINT32> & _viInnerSelectedVariantIds) {
	viShuffleSelect.clear();
	for(UINT32 i = 0;i < _iNum;i++) {
		viShuffleSelect.push_back(i);
	}

	UINT32 iSelectInd = 0;
	_viInnerSelectedVariantIds.clear();
	for(UINT32 i = 0;i < _iNumSelect;i++) {
		iSelectInd = get_random_number(viShuffleSelect.size() - i);
		_viInnerSelectedVariantIds.push_back(viShuffleSelect.at(iSelectInd));
		iter_swap(viShuffleSelect.begin() + iSelectInd, viShuffleSelect.begin() + (viShuffleSelect.size() - i - 1));
	}

}
