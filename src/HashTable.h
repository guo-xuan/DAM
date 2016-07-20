/*
 * HashTable.h
 *
 *  Created on: May 15, 2016
 *      Author: naux
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"

struct KeyValue {
	static UINT32 iNumGroups;
	UINT64 * piKey;
	UINT32 iPiSize;
	UINT32 * iFrequency;

	KeyValue(UINT32 _nKey) {
		piKey = new UINT64[_nKey];
		iPiSize = 0;
		iFrequency = new UINT32[iNumGroups];
	}

	~KeyValue() {
		delete[] iFrequency;
		delete[] piKey;
	}

	KeyValue* operator=(const KeyValue& _a) {
		iPiSize = _a.iPiSize;
		for (size_t i = 0; i < iNumGroups; i++) {
			iFrequency[i] = _a.iFrequency[i];
		}
		for (size_t i = 0; i < iPiSize; i++) {
			piKey[i] = _a.piKey[i];
		}
		return this;
	}

	bool compare(const UINT64* _piKey, UINT32 _iPiSize) {
		if (iPiSize != _iPiSize) {
			return false;
		}
		for (size_t i = 0; i < iPiSize; i++) {
			if (piKey[i] != _piKey[i]) {
				return false;
			}
		}
		return true;
	}

	void clean() {
		iPiSize = 0;
		for (size_t i = 0; i < iNumGroups; i++) {
			iFrequency[i] = 0;
		}
	}

	void set(UINT64 * _piKey, UINT32 _iPiSize) {
		for (size_t i = 0; i < _iPiSize; i++) {
			piKey[i] = _piKey[i];
		}
		iPiSize = _iPiSize;
	}
};

class HashTable {
public:
	UINT32 iKeySize; //number of integers in a key
	vector<KeyValue *> *vkTable;
	vector<UINT32> *viSet;
	UINT32 iSetSize;
	UINT32 iHashCode;
	size_t i;
	size_t j;
	size_t k;
	KeyValue * pKeyValue;
	UINT32 iHashTableSize;
	UINT32 iHashTableSizeMinusOne;

	HashTable(UINT32 _iKeySize, UINT32 _iSampleSize, UINT32 _iNumGroups);
	~HashTable();

	void clean();
	UINT32 extract(vector<vector<UINT32>> & _vviN, vector<vector<double>> & _vvdPrior,
			vector<vector<UINT32>> & _vGroupSet, vector<vector<UINT32*>> _vviVariantFrequence);
	UINT32 hashCode(UINT64 * _piKey, UINT32 _iPiSize);
	void operator=(const HashTable & _a);
	UINT32 push(UINT64 * _piKey, UINT32 _iPiSize);
};

#endif /* HASHTABLE_H_ */
