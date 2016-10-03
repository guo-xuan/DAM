/*
 * DependentModule.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "DependentModule.h"
#include "GwasData.h"

DependentModule::DependentModule(UINT32 _nMaxVariantAllowed, GwasData * _pGwasData, vector<vector<int>> _vvGroupInfo,
		UINT32 _iAssociationType) {
	nMaxVariantAllowed = _nMaxVariantAllowed;
	pGwasData = _pGwasData;
	iNumGroups = pGwasData->iNumGroups;
	nTotalSamples = pGwasData->iNumTotalSamplesAcrossGroups;
	viGroupSizeList = &(pGwasData->viNumSamplesPerGroup);
	vvGroupInfo = _vvGroupInfo;
	iAssociationType = _iAssociationType;
	pHashTable = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups, *viGroupSizeList);
	pHashTableBackup = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups, *viGroupSizeList);
	// pMapTable = new MapTable(nMaxVariantAllowed, nTotalSamples, iNumGroups, *viGroupSizeList);
	dPosteriorProbability = 0;
	dPosteriorProbabilityBackup = 0;
	iNumVariants = 0;
	iNumVariantsBackup = 0;
}

DependentModule::~DependentModule() {
	// delete pMapTable;
	delete pHashTable;
	delete pHashTableBackup;
}

bool DependentModule::addVariants(vector<UINT32> & _viOuterIndex) {
	// string implementation
	for(UINT32 i = 0;i < _viOuterIndex.size();i++) {
		for(UINT32 j = 0;j < iNumGroups;j++) {
			pGwasData->getValue(_viOuterIndex.at(i), j, viTemp);
			// pMapTable->add(i, j, viTemp);
			pHashTable->add(i, j, viTemp);
		}
		viVariants.push_back(_viOuterIndex.at(i));
		vpvviHyperGroupFrequency.push_back(
				(&(pGwasData->getVariantHyperGroupFrequency(_viOuterIndex.at(i))->at(iAssociationType))));
	}
	iNumVariants += _viOuterIndex.size();
	// pMapTable->hashing();
	pHashTable->setNumVariants(iNumVariants);
	pHashTable->clean();
	pHashTable->hashing();
	// pHashTable->print();
	return true;
}

bool DependentModule::apply() {
	pHashTableBackup->copyFrom((*pHashTable));
	dPosteriorProbabilityBackup = dPosteriorProbability;
	iNumVariantsBackup = iNumVariants;
	viVariantsBackup = viVariants;
	vpvviHyperGroupFrequencyBackup = vpvviHyperGroupFrequency;
	return true;
}

double DependentModule::calPosterior() {
	dPosteriorProbability = pHashTable->calPosteriorProbability(vvGroupInfo, vpvviHyperGroupFrequency);
	return dPosteriorProbability;
}

void DependentModule::cleanHashtable() {
	// pMapTable->clean();
	pHashTable->clean();
	iNumVariants = 0;
}

bool DependentModule::delVariants(vector<UINT32> & _viInnerIndex) {
	for(UINT32 i = 0;i < _viInnerIndex.size();i++) {
		pHashTable->del(_viInnerIndex.at(i));
		viVariants.erase(viVariants.begin() + _viInnerIndex.at(i));
		vpvviHyperGroupFrequency.erase(vpvviHyperGroupFrequency.begin() + _viInnerIndex.at(i));
	}
	iNumVariants -= _viInnerIndex.size();
	pHashTable->setNumVariants(iNumVariants);
	pHashTable->clean();
	pHashTable->hashing();
	return true;
}

vector<UINT32> * DependentModule::getVariants(){
	return (& viVariants);
}

void DependentModule::getVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex) {
	_viOuterIndex.clear();
	for(UINT32 i = 0;i < _viInnerIndex.size();i++) {
		_viOuterIndex.push_back(viVariants.at(_viInnerIndex.at(i)));
	}
}

double DependentModule::getPosterior() {
	// dPosteriorProbability = pHashTable->calPosteriorProbability(vvGroupInfo, vpvviHyperGroupFrequency);
	return dPosteriorProbabilityBackup;
}

void DependentModule::getRestVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex) {
	_viOuterIndex.clear();
	for(UINT32 i = 0;i < viVariants.size();i++) {
		if(find(begin(_viInnerIndex), end(_viInnerIndex), viVariants.at(i)) == end(_viInnerIndex)) {
			_viOuterIndex.push_back(viVariants.at(i));
		}
	}
}

bool DependentModule::replaceVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex) {
	if(_viInnerIndex.size() != _viOuterIndex.size()) {
		cout << "Error replace variants." << endl;
	}
	for(UINT32 i = 0;i < _viInnerIndex.size();i++) {
		for(UINT32 j = 0;j < iNumGroups;j++) {
			pGwasData->getValue(_viOuterIndex.at(i), j, viTemp);
			pHashTable->replace(_viInnerIndex.at(i), j, viTemp);
		}
		viVariants.at(_viInnerIndex.at(i)) = _viOuterIndex.at(i);
		vpvviHyperGroupFrequency.at(_viInnerIndex.at(i)) = (&(pGwasData->getVariantHyperGroupFrequency(
				_viOuterIndex.at(i))->at(iAssociationType)));
	}
	pHashTable->clean();
	pHashTable->hashing();
	return true;
}

bool DependentModule::rollBack() {
	pHashTable->copyFrom((*pHashTableBackup));
	dPosteriorProbability = dPosteriorProbabilityBackup;
	iNumVariants = iNumVariantsBackup;
	viVariants = viVariantsBackup;
	vpvviHyperGroupFrequency = vpvviHyperGroupFrequencyBackup;
	return true;
}
