/*
 * HashTable.h
 *
 *  Created on: May 15, 2016
 *      Author: naux
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"

class HashTable {
public:

	static char EMPTY;
	static UINT32 CHARSIZE;
	static UINT32 UINT32SIZE;

	UINT32 iKeySize; //number of integers in a key
	vector<UINT32> viGroupsSize;
	char ** ppcData;
	char * pcHashTable;
	UINT32 * piCount;
	UINT32 iNumGroups;
	UINT32 iTotalSamples;
	UINT32 iHashTableSize;
	UINT32 iHashTableSizeMinusOne;
	UINT32 iNumVariants;
	UINT64 iHashCode;
	UINT32 i;
	UINT32 j;
	UINT32 k;
	UINT32 l;
	UINT64 sum1;
	UINT64 sum2;
	bool bEqual;
	vector<UINT32> viSet;

	double ni;
	double alphai;
	double n_sum;
	double alpha_sum;
	double posterior;

	// vector<vector<char>> vvcData;
	// vector<char> vcHashTable;
	// vector<UINT32> viCount;

	HashTable(UINT32 _iKeySize, UINT32 _iTotalSamples, UINT32 _iNumGroups, vector<UINT32> & _viGroupsSize);
	~HashTable();
	void add(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData);
	double calPosteriorProbability(vector<vector<int>> & _vviHyperGroup, vector<vector<vector<double>> *> & _vpvviHyperGroupFrequency);
	void clean();
	void copyFrom(HashTable & _hashtable);
	void del(UINT32 _iInnerIndex);
	void hashing();
	void print();
	void replace(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData);
	void setNumVariants(UINT32 _iNumVariants);

};

#endif /* HASHTABLE_H_ */
