/*
 * GwasData.h
 *
 *  Created on: May 3, 2016
 *      Author: naux
 */

#ifndef GWASDATA_H_
#define GWASDATA_H_

#include "Config.h"
#include "Variant.h"

class GwasData {
public:

	UINT32 iNumVariants;
	UINT32 iNumTotalSamplesAcrossGroups;
	UINT32 iNumGroups;
	UINT32 iNumVariantTypes;
	// number of samples in each group
	vector<UINT32> viNumSamplesPerGroup;
	// group index for each sample
	vector<UINT32> viSampleGroupInfo;
	// data for each variant
	vector<Variant> vvData;

	vector<string> viGroupNamesSet;

	map<string, UINT32> msiChromsomeIndex;

	// the genotype counts for each variant
	vector<UINT8> viNumVariantTypes;

	static char cDelim;
	static UINT32 iChunkSize;

	GwasData();
	~GwasData();

	void getValue(UINT32 _iVariantIndex, UINT32 _iGroupIndex, vector<UINT32> & _vi);
	UINT32 getVariantChromsome(UINT32 _iVariantIndex);
	vector<vector<vector<double>>> * getVariantHyperGroupFrequency(UINT32 _iVariantIndex);
	vector<vector<vector<UINT32>>> * getVariantHyperGroupOccurrence(UINT32 _iVariantIndex);
	string getVariantName(UINT32 _iVariantIndex);
	void getVariantNearby(UINT32 _iVariantIndex, UINT32 _distance, vector<UINT32> & _vi);
	UINT32 getVariantNumTypes();
	UINT32 getVariantNumTypes(UINT32 _iVariantIndex);
	int getVariantPosition(UINT32 _iVariantIndex);
	bool loadBasicInfo(vector<string> & _vsInputFileList);
	bool loadDataParallel(vector<string> & _vsInputFileList);
	void readInput(vector<string> & _vsInputFileList);
	void setVariantHyperFrequency(UINT32 _iVariantIndex, const vector<vector<vector<int>>>& _vvviHyperGroupInfo);
	void writeOutput(string &_sFilename, const map<string, UINT32> & msiChromsomeIndex);
	void writeOutput(string &_sFilename, vector<UINT32> _viVariantList);
};
#endif /* GWASDATA_H_ */
