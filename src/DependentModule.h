/*
 * DependentModule.h
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#ifndef DEPENDENTMODULE_H_
#define DEPENDENTMODULE_H_

#include "Config.h"
#include "HashTable.h"

class DependentModule {
public:

	HashTable * pHashTable;
	HashTable * pHashTableBackup;
	vector<vector<UINT32>> vvGroupInfo;

	//1D: group, 2D: sample, 3D: variant
	vector<char*> vcGwasData;
	UINT32 iNumKeys;
	vector<char*> vcGwasDataBackup;
	UINT32 iNumKeysBackup;
	UINT32 nMaxVariantAllowed;
	UINT32 nTotalSamples;
	UINT32 iNumGroups;
	GwasData * pGwasData;
	vector<UINT32> * viGroupSizeList;

	DependentModule(UINT32 _nMaxVariantAllowed, GwasData * GwasData, vector<vector<UINT32>> _vvGroupInfo);
	~DependentModule();

	bool addVariants(vector<UINT32> & _viOuterIndex);
	bool delVariants(vector<UINT32> & _viInnerIndex);
	bool replaceVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex);
	bool apply();
	bool rollBack();
	double getPosterior();
};

#endif /* DEPENDENTMODULE_H_ */
