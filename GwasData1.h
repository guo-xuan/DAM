/*
 * GwasData.h
 *
 *  Created on: May 3, 2016
 *      Author: naux
 */

#ifndef GWASDATA_H_
#define GWASDATA_H_

#include "Config.h"

class GwasData {
public:
	vector<string> vsChromsomeNamesSet;
	vector<UINT8> viVariantChromsomeNameList;
	vector<UINT32> viVariantPositionList;
	vector<string> vsVariantNameList;

	UINT32 iNumVariants;
	UINT32 iNumTotalSamplesAcrossGroups;
	UINT32 iNumGroups;
	vector<string> viGroupNamesSet;
	vector<UINT32> viGroupSizeList;
	// the genotype counts for each variant
	vector<UINT8> viNumVariantTypes;
	// group -> variants -> variant type -> frequency
	vector<vector<UINT32*> *> vviVariantFrequence;
	// For each variant, the # samples using 32 bits
	vector<UINT8> viNumSamples32Bits;
	// For each variant, the # 32 bits for each group
	vector<UINT16*> viSampleVectorSizePerGroup;
	// UINT32** is data for each group
	vector<UINT32**> vmDataMatrix; //row variant; column samples
	UINT32 ** miPower; //base exponent

	static char cDelim;
	static UINT32 iMissingData;
	static UINT32 iLogMaxInteger;
	static UINT32 iChunkSize;

	GwasData();
	~GwasData();

	UINT32 getValue(UINT32 _iVariantIndex, UINT32 _iSampleIndex);
	void getValue(UINT32 _iVariantIndex, UINT32 _iGroupIndex, vector<UINT32> & _vi);
	void getVariantFrequency(UINT32 _iVariantIndex, vector<int> & _viGroups, vector<UINT32> & _viFrequency);
	bool loadBasicInfo(vector<string> & _vsInputFileList);
	bool loadDataParallel(vector<string> & _vsInputFileList);
	void readInput(vector<string> & _vsInputFileList);
	void writeOutput(string &_sFilename);
	void writeOutput(string &_sFilename, vector<UINT32> _viVariantList);
};
#endif /* GWASDATA_H_ */
