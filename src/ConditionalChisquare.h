/*
 * ConditionalChisquare.h
 *
 *  Created on: Oct 16, 2016
 *      Author: naux
 */

#ifndef CONDITIONALCHISQUARE_H_
#define CONDITIONALCHISQUARE_H_

#include "GwasData.h"
#include "Config.h"
#include <unordered_map>
#include <utility>
#include "AssociationEvaluation.h"
#include "GammaXG.h"
#include <unordered_set>

struct Interaction;

class ConditionalChisquare {
public:
	ConditionalChisquare();
	~ConditionalChisquare();

	void clean();
	vector<Interaction> & getInteractions();
	void initilize(int _iMaxVariant, GwasData * _gwasData, char ** _pData, UINT32 _iMaxNumVariantTypes);
	vector<Interaction> & calculateSignificance();
	void setVariants(vector<UINT32> & _viVariants);


private:
	vector<vector<vector<int>>> vvviHyperGroup;
	UINT32 iNumGroups;
	vector<UINT32> vGroupSize;
	UINT32 iNumVariantSize;
	double iTotalVariants;
	double dNumComparisons;
	double dCriticalPvalue;

	vector<Interaction> vInteractions;

	/**
	 * level 1: group
	 * level 2: variant
	 * level 3: sample
	 */
	char ** pData;

	vector<UINT32> viVariants;
	vector<bool> viExistVariants;

	UINT32 iMaxNumVariantTypes;

	// inner contingency table
	unordered_map<UINT64, UINT32> umiiKeyTable;
	unordered_map<UINT64, UINT32> umiiLeftKeyTable;
	unordered_map<UINT64, UINT32> umiiRightKeyTable;
	// all frequencies, level 1: group, level 2: left key + right key
	double * pdKeyFrequency;
	// size of pdFrequency
	UINT32 iSizeFrequency;
	double * pdLeftKeyFrequency;
	double * pdRightKeyFrequency;
	// left key & right key
	UINT64 * piKeyPairs;

	// for splitting the variants two left and right
	unordered_set<int> siCodeSet;
	int iCodeUpperBound;
	int iCode;
	int iCodeGroup;
	int iMask;
	vector<UINT32> viVariantsLeft;
	vector<UINT32> viVariantsRight;

	double calculateChiSquare(vector<UINT32> & _viVariants, vector<vector<int>> & _vviAssociationTypes);
	double calculateConditionalChiSquare(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight, vector<vector<int>> & _vviAssociationTypes);
	// collect data and put into a contingency table with given variants
	void collectPartialTablePerLevel(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight);
	// collect data and put into a contingency table
	void collectTablePerLevel(vector<UINT32> & _viVariants);
	// split the variants using the code
	bool getNextVariantSplit(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight);

	void setVariantSplit(vector<UINT32> & _viVariants);
};

#endif /* CONDITIONALCHISQUARE_H_ */
