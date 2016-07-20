/*
 * DependentModule.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "DependentModule.h"
#include "GwasData.h"

DependentModule::DependentModule(UINT32 _nMaxVariantAllowed, GwasData * _pGwasData,
		vector<vector<UINT32>> _vvGroupInfo) {
	nMaxVariantAllowed = _nMaxVariantAllowed;
	pGwasData = _pGwasData;
	iNumGroups = pGwasData->iNumGroups;
	nTotalSamples = pGwasData->iNumTotalSamplesAcrossGroups;
	viGroupSizeList = &(pGwasData->viGroupSizeList);
	for (UINT32 i = 0; i < iNumGroups; i++) {
		char * p = new char[_nMaxVariantAllowed * viGroupSizeList->at(i)];
		vcGwasData.push_back(p);
		p = new char[_nMaxVariantAllowed * viGroupSizeList->at(i)];
		vcGwasDataBackup.push_back(p);
	}
	vvGroupInfo = _vvGroupInfo;
	iNumKeys = 0;
	iNumKeysBackup = 0;
	pHashTable = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups);
	pHashTableBackup = new HashTable(nMaxVariantAllowed, nTotalSamples, iNumGroups);
}

DependentModule::~DependentModule() {
	// TODO Auto-generated destructor stub
}

bool DependentModule::addVariants(vector<UINT32> & _viOuterIndex) {
	for (UINT32 i = 0; i < iNumGroups; i++) {

	}
	return true;
}
