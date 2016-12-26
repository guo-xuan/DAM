/*
 * ConditionalChisquare.cpp
 *
 *  Created on: Oct 16, 2016
 *      Author: naux
 */

#include "ConditionalChisquare.h"

ConditionalChisquare::ConditionalChisquare() {
	pdKeyFrequency = NULL;
	pdLeftKeyFrequency = NULL;
	pdRightKeyFrequency = NULL;
	piKeyPairs = NULL;
	dCriticalPvalue = 0.05;
	iTotalVariants = 0;
	pData = NULL;
	iCode = 0;
	iCodeGroup = 0;
	iCodeUpperBound = 0;
}

ConditionalChisquare::~ConditionalChisquare() {
	if(pdKeyFrequency != NULL) {
		delete[] pdKeyFrequency;
	}
	if(pdLeftKeyFrequency != NULL) {
		delete[] pdLeftKeyFrequency;
	}
	if(pdRightKeyFrequency != NULL) {
		delete[] pdRightKeyFrequency;
	}
	if(piKeyPairs != NULL) {
		delete[] piKeyPairs;
	}
}

void ConditionalChisquare::clean() {
	vInteractions.clear();
}

vector<Interaction> & ConditionalChisquare::getInteractions() {
	return vInteractions;
}

void ConditionalChisquare::initilize(int _iMaxVariant, GwasData * _gwasData, char ** _pData,
		UINT32 _iMaxNumVariantTypes) {
	iNumVariantSize = _iMaxVariant;
	iNumGroups = _gwasData->iNumGroups;
	vGroupSize = _gwasData->viNumSamplesPerGroup;
	iTotalVariants = _gwasData->iNumVariants;
	pData = _pData;
	iMaxNumVariantTypes = _iMaxNumVariantTypes;
	iSizeFrequency = 0;
	for(size_t i = 0;i < vGroupSize.size();++i) {
		if(iSizeFrequency < vGroupSize.at(i)){
			iSizeFrequency = vGroupSize.at(i);
		}
	}
	iSizeFrequency = iSizeFrequency * iNumGroups;
	pdKeyFrequency = new double[iSizeFrequency];
	pdLeftKeyFrequency = new double[iSizeFrequency];
	pdRightKeyFrequency = new double[iSizeFrequency];
	piKeyPairs = new UINT64[iSizeFrequency * 2];
	vvviHyperGroup = HyperGroup::getHyperGroup(iNumGroups);
}

vector<Interaction> & ConditionalChisquare::calculateSignificance() {
	// try all association types first
	// if significant, then check conditional chi-square (all combination of split ways)
	// if also significant, return true after set the p-value and the association type
	int i = 0, iAssociationType = -1;
	double dpValue = 0, dpValueConditional, dMinValueConditional = MINPVALUE, dMaxValueConditional = MAXPVALUE;
	double dpMinValue = 100;
	bool isConditionalSignificant = true;
	for(i = 1;i < (int) vvviHyperGroup.size();++i) {
		dpValue = calculateChiSquare(viVariants, vvviHyperGroup.at(i));
		if(dpValue < dpMinValue) {
			dpMinValue = dpValue;
			iAssociationType = i;
		}
	}
	dMinValueConditional = MINPVALUE;
	dMaxValueConditional = MAXPVALUE;
	if(dpMinValue * dNumComparisons <= dCriticalPvalue) { // significant
		this->setVariantSplit(viVariants);
		isConditionalSignificant = true;
		while(getNextVariantSplit(viVariantsLeft, viVariantsRight)) {
			dpValueConditional = calculateConditionalChiSquare(viVariantsLeft, viVariantsRight,
					vvviHyperGroup.at(iAssociationType));
			// check if this conditional chi-square test is significant
			if(dpValueConditional * dNumComparisons > dCriticalPvalue) { // non-significant
				isConditionalSignificant = false;
				break;
			}
			if(dpValueConditional < dMinValueConditional) {
				dMinValueConditional = dpValueConditional;
			}
			if(dpValueConditional > dMaxValueConditional) {
				dMaxValueConditional = dpValueConditional;
			}
		}
		if(isConditionalSignificant) { // save this interaction
			vInteractions.push_back(Interaction());
			vInteractions.back().viInnerVariantIds = viVariants;
			vInteractions.back().iAssociationType = iAssociationType;
			vInteractions.back().dpValue = dpMinValue;
			vInteractions.back().dMaxpValueConditional = dMaxValueConditional;
			vInteractions.back().dMinpValueConditional = dMinValueConditional;
		}
	}

	return vInteractions;
}

void ConditionalChisquare::setVariants(vector<UINT32> & _viVariants) {
	viVariants = _viVariants;
	dNumComparisons = 1;
	for(int i = 0;i < (int) _viVariants.size();++i) {
		dNumComparisons *= (iTotalVariants - i) / (i + 1);
	}
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

double ConditionalChisquare::calculateChiSquare(vector<UINT32> & _viVariants,
		vector<vector<int>> & _vviAssociationTypes) {
	double dChiSquare = 0;
	this->collectTablePerLevel(_viVariants);
	UINT32 iNum = umiiKeyTable.size(), i, j, k, g, l;
	double dExpect, dExpectVariant, dObserved, dSum = 0;
	UINT64 iKey;
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	double dNumRow = iNum, dNumCol = _vviAssociationTypes.size(), df = 0;
	for(l = 0;l < iNumGroups;++l) {
		dSum += vGroupSize.at(l);
	}
	for(i = 0;i < _vviAssociationTypes.size();++i) {
		for(j = 0;j < iNum;++j) {
			iKey = piKeyPairs[j];
			itFrequency = umiiKeyTable.find(iKey);
			dExpectVariant = 0;
			dExpect = 0;
			dObserved = 0;
			for(l = 0;l < iNumGroups;++l) {
				dExpectVariant += pdKeyFrequency[itFrequency->second * iNumGroups + l];
			}
			for(k = 0;k < _vviAssociationTypes.at(i).size();++k) {
				g = _vviAssociationTypes.at(i).at(k);
				dObserved += pdKeyFrequency[itFrequency->second * iNumGroups + g];
				dExpect += vGroupSize.at(g);
			}
			dExpect = (dExpect * dExpectVariant) / dSum;
			if(dExpect == 0) {
				continue;
			}
			dChiSquare += (dObserved - dExpect) * (dObserved - dExpect) / dExpect;
		}
	}
	df = (dNumRow - 1) * (dNumCol - 1);
	return pValue(df, dChiSquare);
}

double ConditionalChisquare::calculateConditionalChiSquare(vector<UINT32> & _viVariantsLeft,
		vector<UINT32> & _viVariantsRight, vector<vector<int>> & _vviAssociationTypes) {
	double dChiSquare = 0;
	this->collectPartialTablePerLevel(_viVariantsLeft, _viVariantsRight);
	UINT32 iNum = umiiKeyTable.size(), i, j, k, g, s = _viVariantsRight.size();
	UINT64 iKeyGiven, iKeyLeft, iKey;
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	unordered_map<UINT64, UINT32>::const_iterator itLeftKeyFrequency;
	unordered_map<UINT64, UINT32>::const_iterator itRightKeyFrequency;
	double dExpect, dExpectLeft, dExpectRight, dObserved;
	double dNumRow = 0, dNumCol = 0, df = 0, dMinpValue = 100, dpValue;
	for(i = 0;i < _vviAssociationTypes.size();++i) {
		siLeftKeySet.clear();
		siRightKeySet.clear();
		for(j = 0;j < iNum;j++) {
			dExpectLeft = 0;
			dExpectRight = 0;
			dExpect = 0;
			dObserved = 0;
			iKeyLeft = piKeyPairs[j * 2];
			iKeyGiven = piKeyPairs[j * 2 + 1];
			iKey = (iKeyLeft << (iMaxNumVariantTypes * s)) + iKeyGiven;
			itLeftKeyFrequency = umiiLeftKeyTable.find(iKeyLeft);
			itRightKeyFrequency = umiiRightKeyTable.find(iKeyGiven);
			itFrequency = umiiKeyTable.find(iKey);
			for(k = 0;k < _vviAssociationTypes.at(i).size();++k) {
				g = _vviAssociationTypes.at(i).at(k);
				dExpectLeft += pdLeftKeyFrequency[itLeftKeyFrequency->second * iNumGroups + g];
				dExpectRight += pdRightKeyFrequency[itRightKeyFrequency->second * iNumGroups + g];
				dExpect += vGroupSize.at(g);
				dObserved += pdKeyFrequency[itFrequency->second * iNumGroups + g];
			}
			if(dExpectLeft != 0) {
				siLeftKeySet.insert(iKeyLeft);
			}
			if(dExpectRight != 0) {
				siRightKeySet.insert(iKeyGiven);
			}
			dExpect = (dExpectLeft * dExpectRight) / dExpect;
			if(dExpect == 0) {
				continue;
			}
			dChiSquare += (dObserved - dExpect) * (dObserved - dExpect) / dExpect;
		}
		dNumRow = siLeftKeySet.size();
		dNumCol = siRightKeySet.size();
		df = (dNumRow - 1) * (dNumCol - 1);
		dpValue = pValue(df, dChiSquare);
		if(dpValue < dMinpValue) {
			dMinpValue = dpValue;
		}
	}

	return dMinpValue;
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
				iKeyLeft += pData[i][_viVariantsLeft.at(k) * vGroupSize.at(i) + j];
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][_viVariantsLeft.at(k) * vGroupSize.at(i) + j];
			}
			for(k = 0;k < _viVariantsRight.size();++k) {
				iKeyGiven = iKeyGiven << iMaxNumVariantTypes;
				iKeyGiven += pData[i][_viVariantsRight.at(k) * vGroupSize.at(i) + j];
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][_viVariantsRight.at(k) * vGroupSize.at(i) + j];
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
			itFrequency = umiiRightKeyTable.find(iKeyGiven);
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

void ConditionalChisquare::collectTablePerLevel(vector<UINT32> & _viVariants) {
	UINT32 idKeyFrequency = 0;
	UINT64 iKey;
	UINT32 i, j, k;
	unordered_map<UINT64, UINT32>::const_iterator itFrequency;
	memset(pdKeyFrequency, 0, sizeof(double) * iSizeFrequency);
	umiiKeyTable.clear();
	for(i = 0;i < iNumGroups;++i) {
		for(j = 0;j < vGroupSize.at(i);++j) {
			iKey = 0;
			for(k = 0;k < _viVariants.size();++k) {
				iKey = iKey << iMaxNumVariantTypes;
				iKey += pData[i][_viVariants.at(k) * vGroupSize.at(i) + j];
			}
			itFrequency = umiiKeyTable.find(iKey);
			if(itFrequency == umiiKeyTable.end()) {
				umiiKeyTable[iKey] = idKeyFrequency;
				piKeyPairs[idKeyFrequency] = iKey;
				++pdKeyFrequency[idKeyFrequency * iNumGroups + i];
				++idKeyFrequency;
			} else {
				++pdKeyFrequency[itFrequency->second * iNumGroups + i];
			}
		}
	}
}

bool ConditionalChisquare::getNextVariantSplit(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight) {
	bool isGoodSplit = false;
	do {
		_viVariantsLeft.clear();
		_viVariantsRight.clear();
		iCode++;
		if(iCode >= iCodeUpperBound) {
			return false;
		}
		if(siCodeSet.find(iCode) != siCodeSet.end()) {
			isGoodSplit = false;
			iCode++;
			continue;
		}
		int temp = iCode;
		for(int i = 0;i < iCodeGroup;++i) {
			if((temp & 1) == 0) { // to left
				_viVariantsLeft.push_back(viVariants.at(i));
				// _viVariantsLeft.insert(_viVariantsLeft.end(), viVariants.begin() + i, viVariants.begin() + i + 1);
			} else { // to right
				_viVariantsRight.push_back(viVariants.at(i));
				// _viVariantsRight.insert(_viVariantsRight.end(), viVariants.begin() + i, viVariants.begin() + i + 1);
			}
			temp = temp >> 1;
		}
		siCodeSet.insert(iCode);
		temp = (~iCode) & iMask;
		siCodeSet.insert(temp);
		isGoodSplit = true;
	} while(!isGoodSplit);

	return true;
}

void ConditionalChisquare::setVariantSplit(vector<UINT32> & _viVariants) {
	iCode = 0;
	siCodeSet.clear();
	iCodeGroup = _viVariants.size();
	iMask = pow(2, iCodeGroup) - 1;
	iCodeUpperBound = iMask;
}
