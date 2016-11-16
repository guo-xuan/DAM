/*
 * AssociationEvaluation.cpp
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#include "AssociationEvaluation.h"

AssociationEvaluation::AssociationEvaluation() {
	iNumGroups = 0;
	pData = NULL;

}

AssociationEvaluation::~AssociationEvaluation() {
	for(UINT32 j = 0;j < iNumGroups;++j) {
		delete[] pData[j];
	}
	delete[] pData;
}

/**
 *
 */
void AssociationEvaluation::initialize(vector<UINT32> & _viVariants, GwasData * _gwasData,
		UINT32 _iMaxInteractionSize) {
	viVariants = _viVariants;
	iNumVariantCandidates = viVariants.size();
	iNumGroups = _gwasData->iNumGroups;
	viNumSamplesPerGroup = _gwasData->viNumSamplesPerGroup;
	pData = new char*[iNumGroups];
	for(UINT32 i = 0;i < iNumGroups;++i) {
		pData[i] = new char[viNumSamplesPerGroup.at(i) * _viVariants.size()];
	}
	vector<UINT32> vi;
	iMaxNumVariantTypes = 0;
	for(UINT32 i = 0;i < _viVariants.size();++i) {
		if(iMaxNumVariantTypes < _gwasData->getVariantNumTypes(_viVariants.at(i))) {
			iMaxNumVariantTypes = _gwasData->getVariantNumTypes(_viVariants.at(i));
		}
		for(UINT32 j = 0;j < iNumGroups;++j) {
			_gwasData->getValue(_viVariants.at(i), j, vi);
			for(UINT32 k = 0;k < vi.size();++k) {
				pData[j][i * viNumSamplesPerGroup.at(i) + k] = vi.at(k);
			}
		}
	}
	vvviHyperGroup = HyperGroup::getHyperGroup(iNumGroups);
	iBatchSize = 2000;
	// for multi-threading
	iMaxNumThreads = omp_get_max_threads();
	iMaxInteractionSize = _iMaxInteractionSize;
	for(UINT32 i = 0;i < iMaxNumThreads;++i) {
		vChiSquareKits.push_back(ConditionalChisquare());
		vChiSquareKits.back().initilize(iMaxInteractionSize, _gwasData, pData, iMaxNumVariantTypes);
	}
}

void AssociationEvaluation::Evaluation() {
	UINT32 i = 0, j, iNum;

	// start the interaction with 1 variant
	vector<UINT32> vCombination;
	for(i = 1;i <= iMaxInteractionSize;++i) {
		vCombination.clear();
		for(j = 0;j < i;++j) {
			vCombination.push_back(j);
		}
		--vCombination.back();
		while(true) {
			combinationGenerator(i, iNumVariantCandidates, vInteractionCandidates, iBatchSize, vCombination);
			if(vInteractionCandidates.empty()) {
				break;
			}
			// openMP parallel evaluation interaction
			iNum = vInteractionCandidates.size();
#pragma omp for schedule(guided)
			for(UINT32 k = 0;k < iNum;++k) {
				int thread_id = omp_get_thread_num();
				vChiSquareKits.at(thread_id).setVariants(vInteractionCandidates.at(k).viInnerVariantIds, viAlreadySelectedVariants);
				if(vChiSquareKits.at(thread_id).isSignificant()){
					// set the association type, the chi-square value or p-value

				}
			}
		}
	}
}

void AssociationEvaluation::WriteResults(string & _sFilename){

}

/**
 * private methods
 */

bool AssociationEvaluation::combinationGenerator(int _iSize, int _iNumCandidate, vector<Interaction> & _vInteraction,
		int _iBatchSize, vector<UINT32> & _vCombination) {
	int i = 0;
	while(_iBatchSize > 0) {
		++_vCombination.at(_iSize - 1);
		if(((int) _vCombination.at(_iSize - 1)) >= _iNumCandidate) {
			// need to increase higher digit
			int iPos = _iSize - 2;
			for(;iPos >= 0;--iPos) {
				++_vCombination.at(iPos);
				if(((int) _vCombination.at(iPos)) < _iNumCandidate - (_iSize - iPos)) {
					break;
				}
			}
			if(iPos < 0 || ((int) _vCombination.at(iPos)) >= _iNumCandidate - (_iSize - iPos)) {
				// reach the end
				return false;
			}
			for(i = iPos + 1;i < _iSize;++i) {
				_vCombination.at(i) = (_vCombination.at(i - 1) + 1);
			}
		}
		_vInteraction.push_back(Interaction());
		for(i = 0;i < _iSize;++i) {
			_vInteraction.back().viInnerVariantIds.push_back(_vCombination.at(i));
			cout << _vCombination.at(i) << "\t";
		}
		cout << endl;
		--_iBatchSize;
	}
	return true;
}
