/*
 * ConditionalChisquare.cpp
 *
 *  Created on: Oct 16, 2016
 *      Author: naux
 */

#include "ConditionalChisquare.h"

ConditionalChisquare::ConditionalChisquare() {

}

ConditionalChisquare::~ConditionalChisquare() {
	delete[] pdFrequency;
	delete[] piKeyPairs;
}

void ConditionalChisquare::collectPartialTablePerLevel() {
	UINT64 iKeyGiven, iKeyLeft, iKey;
	UINT32 i, j, k;
	memset(pdFrequency, 0, sizeof(double) * iSizeFrequency);
	umiiTable.clear();
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	UINT32 idFrequency = 0;
	for(i = 0;i < iNumGroups;++i) {
		for(j = 0;j < vGroupSize.at(i);++j) {
			iKeyGiven = 0;
			iKeyLeft = 0;
			iKey = 0;
			for(k = 0;k < viExistVariants.size();++k) {
				if(viExistVariants.at(k)) {
					iKeyGiven = iKeyGiven << iMaxNumVariantTypes;
					iKeyGiven += pData[i][viVariants.at(k) * vGroupSize.at(j) + j];
				} else {
					iKeyLeft = iKeyLeft << iMaxNumVariantTypes;
					iKeyLeft += pData[i][viVariants.at(k) * vGroupSize.at(j) + j];
				}
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][viVariants.at(k) * vGroupSize.at(j) + j];
			}
			// do the contingency table
			itFrequency = umiiTable.find(iKey);
			if(itFrequency == umiiTable.end()) {
				umiiTable[iKey] = idFrequency;
				piKeyPairs[idFrequency*2] = iKeyLeft;
				piKeyPairs[idFrequency*2+1] = iKeyGiven;
				++pdFrequency[idFrequency * iNumGroups + i];
				++idFrequency;
			} else {
				++pdFrequency[itFrequency->second * iNumGroups + i];
			}
		}
	}
}

void ConditionalChisquare::initilize(int _iMaxVariant, GwasData * _gwasData, char ** _pData,
		UINT32 _iMaxNumVariantTypes) {
	iNumVariantSize = _iMaxVariant;
	iNumGroups = _gwasData->iNumGroups;
	vGroupSize = _gwasData->viNumSamplesPerGroup;
	pData = _pData;
	iMaxNumVariantTypes = _iMaxNumVariantTypes;
	iSizeFrequency = 0;
	for(size_t i = 0;i < vGroupSize.size();++i) {
		iSizeFrequency += vGroupSize.at(i);
	}
	pdFrequency = new double[iSizeFrequency * iNumGroups];
	piKeyPairs = new UINT64[iSizeFrequency * 2];
}
