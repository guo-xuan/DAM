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
	delete[] pdKeyFrequency;
	delete[] pdLeftKeyFrequency;
	delete[] pdRightKeyFrequency;
	delete[] piKeyPairs;
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
	pdKeyFrequency = new double[iSizeFrequency * iNumGroups];
	pdLeftKeyFrequency = new double[iSizeFrequency * iNumGroups];
	pdRightKeyFrequency = new double[iSizeFrequency * iNumGroups];
	piKeyPairs = new UINT64[iSizeFrequency * 2];
}

bool ConditionalChisquare::isSignificant() {
	// try all association types first
	// if significant, then check conditional chi-square
	// if also significant, return true after set the p-value and the association type

	return false;
}

void ConditionalChisquare::setVariants(vector<UINT32> & _viVariants, vector<bool> & _viExistVariants) {

}

void ConditionalChisquare::setVariants(vector<UINT32> & _viVariants,
		vector<vector<Interaction> > & _viExistInteractions) {

}

/**
 * Private Methods
 */

double pValue(double _degree, double _critical) {
	double p = 0;
	if(_critical / 2 > (_degree / 2 + 1)) {
		p = GammaXG::regularizedGammaQ(_degree / 2, _critical / 2, 10e-50, INT_MAX);
	} else {
		p = 1 - GammaXG::regularizedGammaP(_degree / 2, _critical / 2, 10e-50, INT_MAX);
	}
	return p;
}

double ConditionalChisquare::calculateChiSquare(){
	double dChiSquare = 0;

	return dChiSquare;
}

double ConditionalChisquare::calculateConditionalChiSquare(vector<UINT32> & _viVariantsLeft,
		vector<UINT32> & _viVariantsRight, vector<vector<int>> & _vviAssociationTypes) {
	double dChiSquare = 0;
	this->collectPartialTablePerLevel(_viVariantsLeft, _viVariantsRight);
	UINT32 iNum = umiiKeyTable.size(), i, j, k, g, s = _viVariantsLeft.size();
	UINT64 iKeyGiven, iKeyLeft, iKey;
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	unordered_map<UINT64, UINT32>::const_iterator itLeftKeyFrequency;
	unordered_map<UINT64, UINT32>::const_iterator itRightKeyFrequency;
	double dExpect, dExpectLeft, dExpectRight, dObserved;
	double dNumRow = 0, dNumCol = 0, df = 0;
	for(i = 0;i < _vviAssociationTypes.size();++i) {
		for(j = 0;j < iNum;j++) {
			dExpectLeft = 0;
			dExpectRight = 0;
			dExpect = 0;
			iKeyLeft = piKeyPairs[i * 2];
			iKeyGiven = piKeyPairs[i * 2 + 1];
			iKey = (iKeyLeft << (iMaxNumVariantTypes * s)) + iKeyGiven;
			itLeftKeyFrequency = umiiLeftKeyTable.find(iKeyLeft);
			itRightKeyFrequency = umiiRightKeyTable.find(iKeyGiven);
			itFrequency = umiiKeyTable.find(iKey);
			for(k = 0;k < _vviAssociationTypes.at(i).size();++k) {
				g = _vviAssociationTypes.at(i).at(k);
				dExpectLeft += pdLeftKeyFrequency[itLeftKeyFrequency->second * iNumGroups + g];
				dExpectRight += pdRightKeyFrequency[itRightKeyFrequency->second * iNumGroups + g];
				dExpect += vGroupSize.at(g);
				dObserved += pdLeftKeyFrequency[itFrequency->second * iNumGroups + g];
			}
			if(dExpectLeft != 0) {
				dNumRow++;
			}
			if(dExpectRight != 0) {
				dNumCol++;
			}
			dExpect = (dExpectLeft * dExpectRight) / dExpect;
			if(dExpect == 0){
				continue;
			}
			dChiSquare += (dObserved - dExpect) * (dObserved - dExpect) / dExpect;
		}
	}
	df = (dNumRow - 1) * (dNumCol - 1);
	return pValue(df, dChiSquare);
}

void ConditionalChisquare::collectPartialTablePerLevel(vector<UINT32> & _viVariantsLeft,
		vector<UINT32> & _viVariantsRight) {
	UINT64 iKeyGiven, iKeyLeft, iKey;
	UINT32 i, j, k;
	memset(pdKeyFrequency, 0, sizeof(double) * iSizeFrequency);
	memset(pdLeftKeyFrequency, 0, sizeof(double) * iSizeFrequency);
	memset(pdRightKeyFrequency, 0, sizeof(double) * iSizeFrequency);
	umiiKeyTable.clear();
	umiiLeftKeyTable.clear();
	umiiRightKeyTable.clear();
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	UINT32 idKeyFrequency = 0, idLeftKeyFrequency = 0, idRightKeyFrequency = 0;
	for(i = 0;i < iNumGroups;++i) {
		for(j = 0;j < vGroupSize.at(i);++j) {
			iKeyGiven = 0;
			iKeyLeft = 0;
			iKey = 0;
			for(k = 0;k < _viVariantsLeft.size();++k) {
				iKeyLeft = iKeyLeft << iMaxNumVariantTypes;
				iKeyLeft += pData[i][_viVariantsLeft.at(k) * vGroupSize.at(j) + j];
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][_viVariantsLeft.at(k) * vGroupSize.at(j) + j];
			}
			for(k = 0;k < _viVariantsRight.size();++k) {
				iKeyGiven = iKeyGiven << iMaxNumVariantTypes;
				iKeyGiven += pData[i][_viVariantsRight.at(k) * vGroupSize.at(j) + j];
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][_viVariantsRight.at(k) * vGroupSize.at(j) + j];
			}
			// do the contingency table
			itFrequency = umiiKeyTable.find(iKey);
			if(itFrequency == umiiKeyTable.end()) {
				umiiKeyTable[iKey] = idKeyFrequency;
				piKeyPairs[idKeyFrequency * 2] = iKeyLeft;
				piKeyPairs[idKeyFrequency * 2 + 1] = iKeyGiven;
				++pdKeyFrequency[idKeyFrequency * iNumGroups + i];
				++idKeyFrequency;
			} else {
				++pdKeyFrequency[itFrequency->second * iNumGroups + i];
			}
			// for the left key
			itFrequency = umiiLeftKeyTable.find(iKeyLeft);
			if(itFrequency == umiiLeftKeyTable.end()) {
				umiiLeftKeyTable[iKeyLeft] = idLeftKeyFrequency;
				++pdLeftKeyFrequency[idLeftKeyFrequency * iNumGroups + i];
				++idLeftKeyFrequency;
			} else {
				++pdLeftKeyFrequency[itFrequency->second * iNumGroups + i];
			}
			// for the right key
			itFrequency = umiiRightKeyTable.find(iKeyLeft);
			if(itFrequency == umiiRightKeyTable.end()) {
				umiiRightKeyTable[iKeyGiven] = idRightKeyFrequency;
				++pdRightKeyFrequency[idRightKeyFrequency * iNumGroups + i];
				++idRightKeyFrequency;
			} else {
				++pdRightKeyFrequency[itFrequency->second * iNumGroups + i];
			}
		}
	}
}
