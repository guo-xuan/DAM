/*
 * MapTable.cpp
 *
 *  Created on: Sep 20, 2016
 *      Author: naux
 */

#include "MapTable.h"

MapTable::MapTable(UINT32 _iKeySize, UINT32 _iTotalSamples, UINT32 _iNumGroups, vector<UINT32> & _viGroupsSize) {

	iTotalSamples = _iTotalSamples;
	iNumGroups = _iNumGroups;
	viGroupsSize = _viGroupsSize;
	iNumVariants = 0;
	iNextAvailableCounts = 0;
	i = 0;
	j = 0;
	k = 0;
	viNumCounts.resize(iNumGroups * iTotalSamples, 0);
	for(i = 0;i < iNumGroups;i++) {
		vvsVariants.push_back(vector<string>());
		for(j = 0;j < viGroupsSize.at(i);j++) {
			vvsVariants.at(i).push_back(string());
		}
	}

}

MapTable::~MapTable() {
	// TODO Auto-generated destructor stub
}

void MapTable::add(UINT32 _iInnerIndex, UINT32 _iGroupIndex, vector<UINT32> & _viData) {
	for(UINT32 k = 0;k < _viData.size();k++) {
		vvsVariants.at(_iGroupIndex).at(k).push_back(_viData.at(k));
	}
}

void MapTable::clean() {
	fill(viNumCounts.begin(), viNumCounts.begin() + (iNextAvailableCounts * iNumGroups), 0);
	iNextAvailableCounts = 0;
	for(i = 0;i < iNumGroups;i++) {
		for(j = 0;j < viGroupsSize.at(i);j++) {
			vvsVariants.at(i).at(j).clear();
		}
	}
}

void MapTable::hashing() {
	msuHashTable.clear();
	for(i = 0;i < iNumGroups;i++) {
		for(j = 0;j < viGroupsSize.at(i);j++) {
			it = msuHashTable.find(vvsVariants.at(i).at(j));
			if(it != msuHashTable.end()) {
				viNumCounts.at(it->second * iNumGroups + i)++;}
			else {
				msuHashTable[vvsVariants.at(i).at(j)] = iNextAvailableCounts;
				viNumCounts.at(iNextAvailableCounts * iNumGroups + i)++;iNextAvailableCounts
				++;
			}
		}
	}
}
