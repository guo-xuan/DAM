/*
 * HashTable.cpp
 *
 *  Created on: May 15, 2016
 *      Author: naux
 */

#include "HashTable.h"

char HashTable::EMPTY = 127;
UINT32 HashTable::CHARSIZE = sizeof(char);
UINT32 HashTable::UINT32SIZE = sizeof(UINT32);

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
		memset(ppcData[i], 0, iKeySize * viGroupsSize.at(i) * CHARSIZE);
	}
	iHashTableSize = iTotalSamples * 2;
	iHashTableSizeMinusOne = iHashTableSize - 1;
	// vcHashTable.resize(iHashTableSize * iKeySize, EMPTY);
	pcHashTable = new char[iHashTableSize * iKeySize];
	memset(pcHashTable, EMPTY, iHashTableSize * iKeySize * CHARSIZE);
	// viCount.resize(iHashTableSize * iNumGroups, 0);
	piCount = new UINT32[iHashTableSize * iNumGroups];
	memset(piCount, 0, iHashTableSize * iNumGroups * UINT32SIZE);
	iNumVariants = 0;
	iHashCode = 0;
	i = 0;
	j = 0;
	k = 0;
	l = 0;
	bEqual = false;
	sum1 = 1;
	sum2 = 1;
	ni = 0;
	alphai = 0;
	n_sum = 0;
	alpha_sum = 0;
	posterior = 0;
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

double HashTable::calPosteriorProbability(vector<vector<int>> & _vviHyperGroup, vector<vector<vector<double>> *> & _vpvviHyperGroupFrequency) {
	posterior = 0;
	for(i = 0;i < _vviHyperGroup.size();i++) {
		n_sum = 0;
		alpha_sum = 0;
		for(j = 0;j < viSet.size();j++) {
			ni = 0;
			for(k = 0;k < _vviHyperGroup.at(i).size();k++) {
				ni += piCount[viSet.at(j) * iNumGroups + _vviHyperGroup.at(i).at(k)];
			}
			alphai = 1;
			for(l = 0;l < iNumVariants;l++) {
				alphai *= _vpvviHyperGroupFrequency.at(l)->at(i).at(pcHashTable[viSet.at(j) * iKeySize + l]);
			}
			posterior += lgamma(ni + alphai);
			posterior -= lgamma(alphai);
			n_sum += ni;
			alpha_sum += alphai;
		}
		posterior += lgamma(alpha_sum);
		posterior -= lgamma(n_sum + alpha_sum);
	}
	return posterior;
}

void HashTable::clean() {
	for(i = 0;i < viSet.size();i++) {
		memset(pcHashTable + (viSet.at(i) * iKeySize), EMPTY, iNumVariants * CHARSIZE);
		memset(piCount + (viSet.at(i) * iNumGroups), 0, iNumGroups * UINT32SIZE);
	}
	viSet.clear();
}

void HashTable::copyFrom(HashTable & _hashtable) {
	for(i = 0;i < iNumGroups;i++) {
		memcpy(ppcData[i], _hashtable.ppcData[i], iKeySize * viGroupsSize.at(i) * CHARSIZE);
	}
	memcpy(piCount, _hashtable.piCount, iHashTableSize * iNumGroups * UINT32SIZE);
	viSet = _hashtable.viSet;
	for(i = 0;i < viSet.size();i++) {
		memcpy(pcHashTable + (viSet.at(i) * iKeySize), _hashtable.pcHashTable + (viSet.at(i) * iKeySize), iKeySize * CHARSIZE);
	}
	posterior = _hashtable.posterior;
}

void HashTable::del(UINT32 _iInnerIndex) {
	if(_iInnerIndex == (iNumVariants - 1)) {
		iNumVariants--;
		return;
	}
	k = iNumVariants - _iInnerIndex - 1;
	for(i = 0;i < iNumGroups;i++) {
		for(j = 0;j < viGroupsSize.at(i);j++) {
			memmove(ppcData[i] + (j * iKeySize + _iInnerIndex), ppcData[i] + (j * iKeySize + _iInnerIndex + 1),
					k * CHARSIZE);
		}
	}
}

void HashTable::hashing() {
	for(i = 0;i < iNumGroups;i++) {
		for(j = 0;j < viGroupsSize.at(i);j++) {
			// iHashCode = vvcData.at(i).at(j * iKeySize + 0);
			// iHashCode = ppcData[i][j * iKeySize + 0];
			sum1 = 1;
			sum2 = 1;
			for(k = 0;k < iNumVariants;k++) {
				if(k < 32) {
					sum1 = (sum1 << 2) | (((int) (ppcData[i][j * iKeySize + k])) & 0X03);
				} else {
					sum2 = (sum2 << 2) | (((int) (ppcData[i][j * iKeySize + k])) & 0X03);
				}
				// iHashCode += vvcData.at(i).at(j * iKeySize + k);
				// iHashCode += ppcData[i][j * iKeySize + k];
			}
			iHashCode = ((sum1 % iHashTableSize) * (sum2 % iHashTableSize)) % iHashTableSize;
			// iHashCode %= iHashTableSize;
			// while(vcHashTable.at(iHashCode * iKeySize) != EMPTY) {
			while(pcHashTable[iHashCode * iKeySize] != EMPTY) {
				// check if the codes are the same
				bEqual = true;
				if(memcmp(pcHashTable + (iHashCode * iKeySize), ppcData[i] + (j * iKeySize),
						iNumVariants * CHARSIZE) == 0) {
					bEqual = true;
				} else {
					bEqual = false;
				}
				if(bEqual) {
					break;
				}
				iHashCode = (iHashCode == iHashTableSizeMinusOne) ? 0 : iHashCode + 1;
			}
			// if(vcHashTable.at(iHashCode * iKeySize) == EMPTY) {
			if(pcHashTable[iHashCode * iKeySize] == EMPTY) {
				memcpy(pcHashTable + (iHashCode * iKeySize), ppcData[i] + (j * iKeySize), iNumVariants * CHARSIZE);
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

void HashTable::print() {
	/*	for(i = 0;i < iNumVariants;i++) {
	 for(j = 0;j < iNumGroups;j++) {
	 for(k = 0;k < viGroupsSize.at(j);k++) {
	 cout << ((int)ppcData[j][k * iKeySize + i]) << endl;
	 }
	 }
	 }*/

	for(i = 0;i < viSet.size();i++) {
		for(j = 0;j < iNumVariants;j++) {
			cout << ((int) pcHashTable[viSet.at(i) * iKeySize + j]) << " ";
		}
		cout << "\t";
		for(j = 0;j < iNumGroups;j++) {
			cout << piCount[viSet.at(i) * iNumGroups + j] << "\t";
		}
		cout << endl;
	}
}

void HashTable::replace(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData) {
	if(_viData.size() > viGroupsSize.at(_iGroupIndex) || _iGroupIndex > iNumGroups) {
		cout << "Error: Hash Table replace." << endl;
	}
	for(i = 0;i < _viData.size();i++) {
		// vvcData.at(_iGroupIndex).at(i * iKeySize + _iInnerIndex) = _viData.at(i);
		ppcData[_iGroupIndex][i * iKeySize + _iInnerIndex] = _viData.at(i);
	}
}

void HashTable::setNumVariants(UINT32 _iNumVariants) {
	iNumVariants = _iNumVariants;
}
