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
	vector<string*> vsChromsomeNames;
	vector<UINT16> viVariantChromsomeList;
	vector<UINT32> viVariantPositionList;
	vector<string> vsVariantNameList;

	UINT32 iNumVariants;
	UINT32 iNumSamples;
	UINT32 iNumGroups;
	vector<string> viSampleNamesList;
	vector<UINT32> viGroupSizeList;
	vector<UINT8> viNumVariantTypesList;
	vector<UINT32**> vmDataMatrix;

	static char cDelim;

	GwasData();
	~GwasData();

	bool loadBasicInfo(vector<string> & _vsInputFileList);
	void readInput(vector<string> & _vsInputFileList);
	bool variantControl(vector<string> & _vsList);
	UINT32 getValue(UINT32 _iVariantIndex, UINT32 _iSampleIndex);
};

#endif /* GWASDATA_H_ */
