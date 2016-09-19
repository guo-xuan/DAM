/*
 * DependentModule.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "DependentModule.h"
#include "GwasData.h"

DependentModule::DependentModule(UINT32 _nMaxVariantAllowed, GwasData * _pGwasData, vector<vector<int>> _vvGroupInfo) {
	nMaxVariantAllowed = _nMaxVariantAllowed;
	pGwasData = _pGwasData;
	iNumGroups = pGwasData->iNumGroups;
	nTotalSamples = pGwasData->iNumTotalSamplesAcrossGroups;
	viGroupSizeList = &(pGwasData->viNumSamplesPerGroup);
	for(UINT32 i = 0;i < iNumGroups;i++) {
		char * p = new char[_nMaxVariantAllowed * viGroupSizeList->at(i)];
		vcGwasData.push_back(p);
		p = new char[_nMaxVariantAllowed * viGroupSizeList->at(i)];
		vcGwasDataBackup.push_back(p);
	}
	vvGroupInfo = _vvGroupInfo;
	iNumKeys = 0;
	iNumKeysBackup = 0;
	pHashTable = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups, *viGroupSizeList);
	pHashTableBackup = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups, *viGroupSizeList);
	dPosteriorProbability = 0;
	// string implementation
	iNextAvailableCounts = 0;
	viNumCounts.resize(iNumGroups * nTotalSamples, 0);
	for(UINT32 i = 0;i < iNumGroups;i++) {
		vvsVariants.push_back(vector<string>());
		for(UINT32 j = 0;j < viGroupSizeList->at(i);j++) {
			vvsVariants.at(i).push_back(string());
		}
	}
}

DependentModule::~DependentModule() {

}

bool DependentModule::addVariants(vector<UINT32> & _viOuterIndex) {
	// string implementation
	for(UINT32 i = 0;i < _viOuterIndex.size();i++) {
		for(UINT32 j = 0;j < iNumGroups;j++) {
			pGwasData->getValue(_viOuterIndex.at(i), j, viTemp);
			for(UINT32 k = 0;k < viTemp.size();k++) {
				vvsVariants.at(j).at(k).push_back(viTemp.at(k));
			}
		}
	}
	msuHashTable.clear();
	for(UINT32 i = 0;i < iNumGroups;i++) {
		for(UINT32 j = 0;j < viGroupSizeList->at(i);j++) {
			it = msuHashTable.find(vvsVariants.at(i).at(j));
			if(it != msuHashTable.end()) {
				viNumCounts.at(it->second*iNumGroups+i)++;
			} else {
				msuHashTable[vvsVariants.at(i).at(j)] = iNextAvailableCounts;
				viNumCounts.at(iNextAvailableCounts*iNumGroups+i)++;
				iNextAvailableCounts++;
			}
		}
	}
	// qiuming's hashtable implementation
	return true;
}

void DependentModule::cleanHashtable(){
	fill(viNumCounts.begin(), viNumCounts.begin()+(iNextAvailableCounts*iNumGroups), 0);
	iNextAvailableCounts = 0;
}

void DependentModule::cleanVariants(){
	for(UINT32 i = 0;i < iNumGroups;i++) {
		for(UINT32 j = 0;j < viGroupSizeList->at(i);j++) {
			vvsVariants.at(i).at(j).clear();
		}
	}
}

double DependentModule::getPosterior() {
	return dPosteriorProbability;
}
