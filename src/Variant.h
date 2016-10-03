/*
 * Variant.h
 *
 *  Created on: Sep 4, 2016
 *      Author: naux
 */

#ifndef VARIANT_H_
#define VARIANT_H_

#include "Config.h"
#include <stdlib.h>

class Variant {
public:

	static UINT32 iLogMaxInteger;
	static UINT32 iMissingData;
	static UINT32 ** miPower; //base exponent
	static vector<UINT32> viNumSamplePerGroup;

	// string sChromsomeName;
	UINT32 iChromsomeNameIndex;
	string sVariantName;
	UINT32 iVariantPosition;

	// number of types (a bi-allele snp has 3 genotypes)
	UINT32 iNumTypes;
	UINT32 iNumSamplesPerInteger;

	Variant();
	~Variant();

	void getDataInGroup(UINT32 iGroupIndex, vector<UINT32> & _vi);
	vector<vector<vector<double>>> * getHyperGroupFrequency();
	vector<vector<vector<UINT32>>> * getHyperGroupOccurrence();
	void setData(const string & _sLine, char _delim, const vector<UINT32> & _viNumSamplePerGroup,
			const vector<UINT32> & _viSampleGroupInfo, const map<string, UINT32> & msiChromsomeIndex);
	void setData(const string & _sLine, char _delim, const vector<UINT32> & _viNumSamplePerGroup,
			const vector<UINT32> & _viSampleGroupInfo, vector<UINT32>& _vector, map<UINT32, UINT32> & _mii,
			vector<UINT32> & _viGroupCount, const map<string, UINT32> & msiChromsomeIndex);
	void setFrequency(const vector<vector<vector<int>>> & _vvviHyperGroupInfo);

private:
	/**
	 * data in each group of samples
	 * Level 1: group index
	 * Level 2: sample/pack_integer index
	 * Level 3: data
	 */
	vector<vector<UINT32>> vviData;
	/**
	 * Level 1: group index
	 * Level 2: variant type index
	 * Level 3: occurrence of a variant type
	 */
	vector<vector<UINT32>> vviOccurrence;

	/**
	 * Level 1: association type index
	 * Level 2: hyper group index
	 * Level 3: frequency of a variant type
	 */
	vector<vector<vector<double>>> vvviHyperGroupFrequency;

	/**
	 * Level 1: association type index
	 * Level 2: hyper group index
	 * Level 3: occurrence of a variant type
	 */
	vector<vector<vector<UINT32>>> vvviHyperGroupOccurrence;

	UINT32 getNumSamplesPerInteger(double base);
};

#endif /* VARIANT_H_ */
