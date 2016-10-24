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

class ConditionalChisquare {
public:
	ConditionalChisquare();
	~ConditionalChisquare();

	bool isSignificant();
	void initilize(int _iMaxVariant, GwasData * _gwasData, char ** _pData, UINT32 _iMaxNumVariantTypes);
	void setVariants(vector<UINT32> & _viVariants, vector<bool> & _viExistVariants);

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
	unordered_map<UINT64, UINT32> umiiTable;
	double * pdFrequency;
	UINT32 iSizeFrequency;
	UINT64 * piKeyPairs;

	double calculateChiSquare();
	double calculateConditionalChiSquare();
	void collectPartialTablePerLevel();
};

#endif /* CONDITIONALCHISQUARE_H_ */
