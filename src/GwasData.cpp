/*
 * GwasData.cpp
 *
 *  Created on: May 3, 2016
 *      Author: naux
 */

#include "GwasData.h"

char GwasData::cDelim = '\t';
UINT32 GwasData::iChunkSize = 5000;

GwasData::GwasData() {
	iNumVariants = 0;
	iNumGroups = 0;
	iNumTotalSamplesAcrossGroups = 0;
}

GwasData::~GwasData() {

	for (size_t i = 0; i < vviVariantFrequence.size(); i++) {
		for (size_t j = 0; j < vviVariantFrequence.at(i)->size(); j++) {
			delete[] vviVariantFrequence.at(i)->at(j);
		}
		vviVariantFrequence.at(i)->clear();
		delete vviVariantFrequence.at(i);
	}
}

void split(const string& sLine, char delim, vector<string>& _vector) {
	size_t i = 0;
	size_t pos = sLine.find(delim);
	_vector.clear();
	while (pos != string::npos) {
		_vector.push_back(sLine.substr(i, pos - i));
		i = ++pos;
		pos = sLine.find(delim, pos);

		if (pos == string::npos)
			_vector.push_back(sLine.substr(i, sLine.length()));
	}
}

bool split(const string& sLine, char delim, int _iSeperatePos, vector<string> & _vsVariantInfoList,
		vector<UINT32>& _viValueList) {
	size_t i = 0;
	size_t pos = sLine.find(delim);
	_viValueList.clear();
	_vsVariantInfoList.clear();
	int iCount = 0;
	for (; iCount < _iSeperatePos; iCount++) {
		_vsVariantInfoList.push_back(sLine.substr(i, pos - i));
		i = ++pos;
		pos = sLine.find(delim, pos);
		if (pos == string::npos) {
			return false;
		}
	}
	while (pos != string::npos) {
		_viValueList.push_back(stoi(sLine.substr(i, pos - i)));
		i = ++pos;
		pos = sLine.find(delim, pos);
	}
	if (pos == string::npos) {
		_viValueList.push_back(stoi(sLine.substr(i, sLine.length())));
	}
	return true;
}

void GwasData::getValue(UINT32 _iVariantIndex, UINT32 _iGroupIndex, vector<UINT32> & _vi) {
	vvData.at(_iVariantIndex).getDataInGroup(_iGroupIndex, _vi);
}

void GwasData::getVariantFrequency(UINT32 _iVariantIndex, vector<int> & _viGroups, vector<UINT32> & _viFrequency){
	UINT32 iFrequency = 0;
	_viFrequency.clear();
	for(UINT32 i = 0; i<viNumVariantTypes.at(_iVariantIndex);i++){
		iFrequency = 0;
		for(UINT32 j=0;j<_viGroups.size();j++){
			iFrequency += ((vviVariantFrequence.at(_viGroups.at(j)))->at(_iVariantIndex))[i];
		}
		_viFrequency.push_back(iFrequency);
	}
}

/**
 * Level 1: association type index
 * Level 2: hyper group index
 * Level 3: frequency of a variant type
 */
vector<vector<vector<double>>> * GwasData::getVariantHyperGroupFrequency(UINT32 _iVariantIndex){
	return vvData.at(_iVariantIndex).getHyperGroupFrequency();
}

/**
 * Level 1: association type index
 * Level 2: hyper group index
 * Level 3: occurrence of a variant type
 */
vector<vector<vector<UINT32>>> * GwasData::getVariantHyperGroupOccurrence(UINT32 _iVariantIndex){
	return vvData.at(_iVariantIndex).getHyperGroupOccurrence();
}

UINT32 GwasData::getVairantNumTypes(UINT32 _iVariantIndex){
	return vvData.at(_iVariantIndex).iNumTypes;
}

bool GwasData::loadBasicInfo(vector<string> & _vsInputFileList) {
	cout << "Initializing Space" << endl;
	ifstream myFile;
	string text;
	vector<string> vs;
	vector<UINT32> vi;
	//open the first file
	myFile.open(_vsInputFileList.at(0).c_str());
	if (!myFile.is_open()) {
		cout << "Unable to open file: " << _vsInputFileList.at(0) << endl;
		return false;
	}
	//read the first line
	string sGroupInfoHeader;
	getline(myFile, text);
	sGroupInfoHeader = text;
	split(text, cDelim, vs);
	map<string, int> msi;
	for (size_t i = 3; i < vs.size(); i++) {
		if(msi.find(vs.at(i))!=msi.end()){
			msi[vs.at(i)] = msi[vs.at(i)] + 1;
		}else{
			msi[vs.at(i)] = 1;
		}
		iNumTotalSamplesAcrossGroups++;
	}
	iNumGroups = msi.size();
	myFile.close();
	viNumSamplesPerGroup.resize(iNumGroups, 0);
	viSampleGroupInfo.resize(iNumTotalSamplesAcrossGroups, 0);
	UINT32 iValue = 0;
	for (size_t i = 3; i < vs.size(); i++) {
		iValue = atoi(vs.at(i).c_str());
		if(iValue>=viNumSamplesPerGroup.size()){
			cout << "Group Index Error." << endl;
		}
		viSampleGroupInfo.at(i-3) = iValue;
		viNumSamplesPerGroup.at(iValue)++;
	}
	//check the group info and the number of variant types per variant for all input files
	for (size_t i = 0; i < _vsInputFileList.size(); i++) {
		myFile.open(_vsInputFileList.at(i).c_str());
		if (!myFile.is_open()) {
			cout << "Unable to open file: " << _vsInputFileList.at(i) << endl;
			return false;
		}
		//check the group information
		getline(myFile, text);
		if (text != sGroupInfoHeader) {
			cout << "Group info is not consistent in file: " << _vsInputFileList.at(i) << endl;
			myFile.close();
			return false;
		}
		//check the variant information
		while (!myFile.eof()) {
			getline(myFile, text);
			if (text.empty()) {
				continue;
			}
			split(text, cDelim, 3, vs, vi);
			if (vi.size() != iNumTotalSamplesAcrossGroups) {
				cout << "Number of variants is not consistent across files: " << _vsInputFileList.at(i) << endl;
				cout << "Error line: " << text << endl;
				myFile.close();
				return false;
			} else {
				iNumVariants++;
			}

		}
		myFile.close();
	}
	vvData.resize(iNumVariants);
	Variant::viNumSamplePerGroup = viNumSamplesPerGroup;
	return true;
}

bool GwasData::loadDataParallel(vector<string> & _vsInputFileList) {
	ifstream myFile;
	string text;
	UINT32 chunkCounter;
	vector<string>* lineData = new vector<string>();
	UINT32 iVariantIndex = 0;
	vector<string> vsGroupInfo;
	for (size_t i = 0; i < _vsInputFileList.size(); i++) {
		myFile.open(_vsInputFileList.at(i).c_str());
		if (!myFile.is_open()) {
			cout << "Unable to open file: " << _vsInputFileList.at(i) << endl;
			return false;
		}
		// skip the header line
		getline(myFile, text);
		// read each variant
		while (!myFile.eof()) {
			//buffer the data
			lineData->clear();
			for (chunkCounter = 0; !myFile.eof() && chunkCounter < iChunkSize;) {
				getline(myFile, text);
				text.erase(std::remove(text.begin(), text.end(), '\n'), text.end());
				if (text.empty()) {
					break;
				}
				lineData->push_back(text);
				chunkCounter++;
			}
			//parallel parse data
//#pragma omp parallel
			{
				string sLine;
				vector<string>::const_iterator vs_iter;
				vector<UINT32> viGroupCount(iNumGroups, 0);
				vector<UINT32> vi;
				map<UINT32, UINT32> mii;
//#pragma omp for schedule(guided)
				for (UINT32 j = 0; j < chunkCounter; j++) {
					Variant variant;
					vvData.at(iVariantIndex + j) = variant;
					vvData.at(iVariantIndex + j).setData(lineData->at(j), cDelim, viNumSamplesPerGroup, viSampleGroupInfo, vi, mii, viGroupCount);
				}
			}
			iVariantIndex += chunkCounter;
			cout << "Reading variants %" << 100 * ((double) iVariantIndex) / ((double) iNumVariants) << "\r";
		}
		myFile.close();
	}
	cout << endl;
	delete lineData;
	return true;
}

void GwasData::readInput(vector<string> & _vsInputFileList) {
	loadBasicInfo(_vsInputFileList);
	loadDataParallel(_vsInputFileList);
}

void GwasData::setVariantHyperFrequency(UINT32 _iVariantIndex, const vector<vector<vector<int>>>& _vvviHyperGroupInfo){
	vvData.at(_iVariantIndex).setFrequency(_vvviHyperGroupInfo);
}

void GwasData::writeOutput(string & _sFilename) {
	ofstream filePointer;
	filePointer.open(_sFilename.c_str(), ios_base::out);
	vector<UINT32> vi;
	//header
	filePointer << "ID\tChr\tPos";
	for (size_t j = 0; j < iNumGroups; j++) {
		for (size_t i = 0; i < viNumSamplesPerGroup.at(j); i++) {
			filePointer << "\t" << j;
		}
	}
	filePointer << "\n";
	for (size_t i = 0; i < iNumVariants; i++) {
		filePointer << vvData.at(i).sVariantName << "\t";
		filePointer << vvData.at(i).sChromsomeName << "\t";
		filePointer << vvData.at(i).iVariantPosition;
		for (size_t j = 0; j < iNumGroups; j++) {
			getValue(i, j, vi);
			for (size_t k = 0; k < vi.size(); k++) {
				filePointer << "\t";
				filePointer << vi.at(k);
			}
		}
		filePointer << "\n";
	}
	filePointer.close();

}
