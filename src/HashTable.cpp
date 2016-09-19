/*
 * HashTable.cpp
 *
 *  Created on: May 15, 2016
 *      Author: naux
 */

#include "HashTable.h"

char HashTable::EMPTY = 127;

HashTable::HashTable(UINT32 _iKeySize, UINT32 _iTotalSamples, UINT32 _iNumGroups, vector<UINT32> & _viGroupsSize) {
	iKeySize = _iKeySize;
	iNumGroups = _iNumGroups;
	viGroupsSize = _viGroupsSize;
	iTotalSamples = _iTotalSamples;
	ppcData = new char*[iNumGroups];
	for(i = 0;i < iNumGroups;i++) {
		// vvcData.push_back(vector<char>());
		// vvcData.at(i).resize(iKeySize * viGroupsSize.at(i), 0);
		ppcData[i] = new char[iKeySize * viGroupsSize.at(i)];
		memset(ppcData[i], 0, iKeySize * viGroupsSize.at(i) * sizeof(char));
	}
	iHashTableSize = iTotalSamples * 2;
	iHashTableSizeMinusOne = iHashTableSize - 1;
	// vcHashTable.resize(iHashTableSize * iKeySize, EMPTY);
	pcHashTable = new char[iHashTableSize * iKeySize];
	memset(pcHashTable, 0, iHashTableSize * iKeySize * sizeof(char));
	// viCount.resize(iHashTableSize * iNumGroups, 0);
	piCount = new UINT32[iHashTableSize * iNumGroups];
	memset(piCount, 0, iHashTableSize * iNumGroups * sizeof(UINT32));
	iNumVariants = 0;
	iHashCode = 0;
	i = 0;
	j = 0;
	k = 0;
	bEqual = false;
}

HashTable::~HashTable() {
	for(i = 0;i < iNumGroups;i++) {
		delete[] ppcData[i];
	}
	delete[] ppcData;
	delete[] pcHashTable;
	delete[] piCount;
}

void HashTable::add(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData) {
	if(_viData.size() > viGroupsSize.at(_iGroupIndex) || _iGroupIndex > iNumGroups) {
		cout << "Error: Hash Table add." << endl;
	}
	for(i = 0;i < _viData.size();i++) {
		// vvcData.at(_iGroupIndex).at(i * iKeySize + _iInnerIndex) = _viData.at(i);
		ppcData[_iGroupIndex][i * iKeySize + _iInnerIndex] = _viData.at(i);
	}
}

void HashTable::clean() {
	for(i = 0;i < viSet.size();i++) {
		memset(pcHashTable + (viSet.at(i) * iKeySize), 0, iNumVariants * sizeof(char));
		memset(piCount + (viSet.at(i) * iNumGroups), 0, iNumGroups * sizeof(UINT32));
	}
	viSet.clear();
}

void HashTable::hashing() {
	for(i = 0;i < iNumGroups;i++) {
		for(j = 0;j < viGroupsSize.at(i);j++) {
			// iHashCode = vvcData.at(i).at(j * iKeySize + 0);
			iHashCode = ppcData[i][j * iKeySize + 0];
			for(k = 1;k < iNumVariants;k++) {
				// iHashCode += vvcData.at(i).at(j * iKeySize + k);
				iHashCode += ppcData[i][j * iKeySize + k];
			}
			iHashCode %= iHashTableSize;
			// while(vcHashTable.at(iHashCode * iKeySize) != EMPTY) {
			while(pcHashTable[iHashCode * iKeySize] != EMPTY) {
				// check if the codes are the same
				bEqual = true;
				if(memcmp(pcHashTable + (iHashCode * iKeySize + k), ppcData[i] + (j * iKeySize + k),
						iNumVariants * sizeof(char)) == 0) {
					bEqual = true;
				} else {
					bEqual = false;
				}
				/*for(k = 0;k < iNumVariants;k++) {
				 // if(vcHashTable.at(iHashCode * iKeySize + k) != vvcData.at(i).at(j * iKeySize + k)) {
				 if(pcHashTable[iHashCode * iKeySize + k] != ppcData[i][j * iKeySize + k]) {
				 bEqual = false;
				 break;
				 }
				 }*/
				if(bEqual) {
					break;
				}
				iHashCode = (iHashCode == iHashTableSizeMinusOne) ? 0 : iHashCode + 1;
			}
			// if(vcHashTable.at(iHashCode * iKeySize) == EMPTY) {
			if(pcHashTable[iHashCode * iKeySize] == EMPTY) {
				memcpy(pcHashTable + (iHashCode * iKeySize + k), ppcData[i] + (j * iKeySize + k), iNumVariants * sizeof(char) );
				/*
				for(k = 0;k < iNumVariants;k++) {
					// vcHashTable.at(iHashCode * iKeySize + k) = vvcData.at(i).at(j * iKeySize + k);
					pcHashTable[iHashCode * iKeySize + k] = ppcData[i][j * iKeySize + k];
				}*/
				viSet.push_back(iHashCode);
			}
			// viCount.at(iHashCode * iNumGroups + i) += 1;
			piCount[iHashCode * iNumGroups + i] += 1;
		}
	}
}

void HashTable::setNumVariants(UINT32 _iNumVariants) {
	iNumVariants = _iNumVariants;
}
