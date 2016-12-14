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

struct Interaction;

class ConditionalChisquare {
public:
	ConditionalChisquare();
	~ConditionalChisquare();

	void initilize(int _iMaxVariant, GwasData * _gwasData, char ** _pData, UINT32 _iMaxNumVariantTypes);
	bool isSignificant();
	void setVariants(vector<UINT32> & _viVariants, vector<bool> & _viExistVariants);
	void setVariants(vector<UINT32> & _viVariants, vector<vector<Interaction> > & _viExistInteractions);


private:
	UINT32 iNumGroups;
	vector<UINT32> vGroupSize;
	UINT32 iNumVariantSize;
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

	double calculateChiSquare();
	double calculateConditionalChiSquare(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight, vector<vector<int>> & _vviAssociationTypes);
	void collectPartialTablePerLevel(vector<UINT32> & _viVariantsLeft, vector<UINT32> & _viVariantsRight);
};

#endif /* CONDITIONALCHISQUARE_H_ */
