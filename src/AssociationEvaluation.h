/*
 * AssociationEvaluation.h
 *
 *  Created on: May 5, 2016
 *      Author: naux
 */

#ifndef ASSOCIATIONEVALUATION_H_
#define ASSOCIATIONEVALUATION_H_

#include "GwasData.h"
#include "Config.h"
#include "HyperGroup.h"
#include "ConditionalChisquare.h"

#include <unordered_set>

struct Interaction {
	vector<UINT32> viInnerVariantIds;
	UINT32 iAssociationType;

};

class AssociationEvaluation {
public:
	AssociationEvaluation();
	~AssociationEvaluation();

	void initialize(vector<UINT32> & _viVariants, GwasData * _gwasData, UINT32 _iMaxInteractionSize);
	void Evaluation();

private:
	/**
	 * level 1: group
	 * level 2: variant
	 * level 3: sample
	 */
	char ** pData;
	vector<UINT32> viVariants;
	int iNumVariantCandidates;
	UINT32 iNumGroups;
	vector<UINT32> viNumSamplesPerGroup;
	vector<vector<vector<int>>> vvviHyperGroup;

	vector<Interaction> vInteractions;
	vector<Interaction> vInteractionCandidates;
	unordered_set<UINT32> sUsedVariantInnerIds;

	int iBatchSize;
	UINT32 iMaxNumThreads;
	UINT32 iMaxInteractionSize;

	vector<ConditionalChisquare> vChiSquareKits;

	UINT32 iMaxNumVariantTypes;


	bool combinationGenerator(int _iSize, int _iNumCandidate, vector<Interaction> & _vInteraction, int _iBatchSize, vector<UINT32> & _vCombination);
};

#endif /* ASSOCIATIONEVALUATION_H_ */
