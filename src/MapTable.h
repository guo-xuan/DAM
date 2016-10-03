/*
 * MapTable.h
 *
 *  Created on: Sep 20, 2016
 *      Author: naux
 */

#ifndef MAPTABLE_H_
#define MAPTABLE_H_

#include "Config.h"
#include "GwasData.h"

class MapTable {
public:

	UINT32 iNumGroups;
	UINT32 iTotalSamples;
	UINT32 iNumVariants;
	UINT32 i;
	UINT32 j;
	UINT32 k;
	vector<UINT32> viGroupsSize;

	vector<UINT32> viNumCounts;
	UINT32 iNextAvailableCounts;
	vector<vector<string>> vvsVariants;
	map<string, UINT32> msuHashTable;
	map<string, UINT32>::iterator it;

	MapTable(UINT32 _iKeySize, UINT32 _iTotalSamples, UINT32 _iNumGroups, vector<UINT32> & _viGroupsSize);
	~MapTable();

	void add(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData);
	void clean();
	void del(UINT32 _iInnerIndex);
	void extract(vector<vector<UINT32>> & _vviHyperGroup, vector<UINT32> & _viCount);
	void hashing();
	void replace(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData);
	void setNumVariants(UINT32 _iNumVariants);

};

#endif /* MAPTABLE_H_ */
