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
#include "MapTable.h"

class DependentModule {
public:

	double dPosteriorProbability;
	double dPosteriorProbabilityBackup;

	HashTable * pHashTable;
	HashTable * pHashTableBackup;
	vector<vector<int>> vvGroupInfo;
	UINT32 iAssociationType;

	UINT32 nMaxVariantAllowed;
	UINT32 nTotalSamples;
	UINT32 iNumGroups;
	GwasData * pGwasData;
	vector<UINT32> * viGroupSizeList;
	vector<UINT32> viTemp;

	// MapTable * pMapTable;

	UINT32 iNumVariants;
	vector<UINT32> viVariants;
	vector<vector<vector<double>> *> vpvviHyperGroupFrequency;

	UINT32 iNumVariantsBackup;
	vector<UINT32> viVariantsBackup;
	vector<vector<vector<double>> *> vpvviHyperGroupFrequencyBackup;

	DependentModule(UINT32 _nMaxVariantAllowed, GwasData * GwasData, vector<vector<int>> _vvGroupInfo, UINT32 _iAssociationType);
	~DependentModule();

	bool addVariants(vector<UINT32> & _viOuterIndex);
	bool apply();
	double calPosterior();
	void cleanHashtable();
	bool delVariants(vector<UINT32> & _viInnerIndex);
	vector<UINT32> * getVariants();
	void getVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex);
	double getPosterior();
	void getRestVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex);
	bool replaceVariants(vector<UINT32> & _viInnerIndex, vector<UINT32> & _viOuterIndex);
	bool rollBack();
};

#endif /* DEPENDENTMODULE_H_ */
