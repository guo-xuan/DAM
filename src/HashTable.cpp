/*
 * HashTable.cpp
 *
 *  Created on: May 15, 2016
 *      Author: naux
 */

#include "HashTable.h"

UINT32 KeyValue::iNumGroups = 0;

HashTable::HashTable(UINT32 _iKeySize, UINT32 _iSampleSize, UINT32 _iNumGroups) {
	iKeySize = _iKeySize;
	iHashTableSize = _iSampleSize;
	iHashTableSizeMinusOne = iHashTableSize - 1;
	viSet = new vector<UINT32>(iHashTableSize);
	iSetSize = 0;
	vkTable = new vector<KeyValue *>(_iSampleSize);
	for (size_t i = 0; i < vkTable->size(); i++) {
		vkTable->at(i) = new KeyValue(iKeySize);
	}
	iHashCode = 0;
	i = 0;
	j = 0;
	k = 0;
	pKeyValue = NULL;
	if (KeyValue::iNumGroups == 0) {
		KeyValue::iNumGroups = _iNumGroups;
	} else if (KeyValue::iNumGroups != _iNumGroups) {
		cerr << "Group number is inconsistent" << endl;
		exit(1);
	}
}

HashTable::~HashTable() {
	for (i = 0; i < vkTable->size(); i++) {
		delete vkTable->at(i);
	}
	delete vkTable;
	delete viSet;
}

void HashTable::clean() {
	for (i = 0; i < iSetSize; i++) {
		vkTable->at(viSet->at(i))->clean();
	}
	iSetSize = 0;
}

UINT32 HashTable::extract(vector<vector<UINT32>> & _vviN, vector<vector<double>> & _vvdPrior,
		vector<vector<UINT32>> & _vGroupSet, vector<vector<UINT32*>> _vviVariantFrequence) {
	UINT32 temp = 0;
	for (i = 0; i < _vGroupSet.size(); i++) { //super group
		_vviN.at(i).clear();
		_vvdPrior.at(i).clear();
	}
	for (i = 0; i < iSetSize; i++) {
		pKeyValue = vkTable->at(viSet->at(i));
		for (j = 0; j < _vGroupSet.size(); j++) { //super group
			temp = 0;
			for (k = 0; k < _vGroupSet.at(j).size(); k++) {
				temp += pKeyValue->iFrequency[_vGroupSet.at(j).at(k)];
			}
			if (temp > 0) {
				_vviN.at(j).push_back(temp);

			}
		}
	}
	return iSetSize;
}

UINT32 HashTable::hashCode(UINT64 * _piKey, UINT32 _iPiSize) {
	iHashCode = 0;
	iHashCode = _piKey[0];
	for (i = 1; i < _iPiSize; i++) {
		iHashCode += _piKey[i];
	}
	iHashCode %= iHashTableSize;
	return iHashCode;
}

void HashTable::operator=(const HashTable & _a) {
	iSetSize = _a.iSetSize;
	for (i = 0; i < iSetSize; i++) {
		j = _a.viSet->at(i);
		viSet->at(i) = j;
		pKeyValue = vkTable->at(j);
		(*pKeyValue) = (*(_a.vkTable->at(j)));
	}
}

UINT32 HashTable::push(UINT64 * _piKey, UINT32 _iPiSize) {
	iHashCode = hashCode(_piKey, _iPiSize);
	while (1) {
		pKeyValue = vkTable->at(iHashCode);
		if (pKeyValue->iPiSize == 0) {
			pKeyValue->set(_piKey, _iPiSize);
			viSet->at(iSetSize) = iHashCode;
			iSetSize++;
			break;
		} else if (pKeyValue->compare(_piKey, _iPiSize)) {
			break;
		} else {
			iHashCode = (iHashCode == iHashTableSizeMinusOne) ? 0 : iHashCode + 1;
		}
	}
	return iHashCode;
}
