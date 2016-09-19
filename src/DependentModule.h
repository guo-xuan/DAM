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
#include "GwasData.h"

class DependentModule {
public:

	double dPosteriorProbability;

	HashTable * pHashTable;
	HashTable * pHashTableBackup;
	vector<vector<int>> vvGroupInfo;

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

	// string implementation
	vector<UINT32> viNumCounts;
	UINT32 iNextAvailableCounts;
	vector<vector<string>> vvsVariants;
	map<string, UINT32> msuHashTable;
	vector<UINT32> viTemp;
	map<string, UINT32>::iterator it;

	DependentModule(UINT32 _nMaxVariantAllowed, GwasData * GwasData, vector<vector<int>> _vvGroupInfo);
	~DependentModule();

	bool addVariants(vector<UINT32> & _viOuterIndex);
	bool apply();
	void cleanHashtable();
	void cleanVariants();
	bool delVariants(vector<UINT32> & _viInnerIndex);
	double getPosterior();
	bool replaceVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex);
	bool rollBack();
};

#endif /* DEPENDENTMODULE_H_ */
