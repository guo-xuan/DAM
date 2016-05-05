/*
 * GwasData.cpp
 *
 *  Created on: May 3, 2016
 *      Author: naux
 */

#include "GwasData.h"

char GwasData::cDelim = '\t';

GwasData::GwasData() {
	iNumVariants = 0;
	iNumGroups = 0;
	iNumSamples = 0;
}

GwasData::~GwasData() {
	size_t len = this->vsChromsomeNames.size();
	for (size_t i = 0; i < len; i++) {
		delete this->vsChromsomeNames.at(i);
	}
	for (size_t i = 0; i < this->iNumGroups; i++) {
		for (size_t j = 0; j < this->iNumVariants; j++) {
			delete[] this->vmDataMatrix.at(i)[j];
		}
		delete[] this->vmDataMatrix.at(i);
	}
}

void GwasData::readInput(vector<string> & _vsInputFileList) {

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

bool GwasData::loadBasicInfo(vector<string> & _vsInputFileList) {
	cout << "Initializing Space" << endl;
	ifstream myFile;
	string text;
	string element;
	vector<string> vs;
	vector<string>::const_iterator i_iter;
	myFile.open(_vsInputFileList.at(0).c_str());
	if (!myFile.is_open()) {
		cout << "Unable to open file: " << _vsInputFileList.at(0) << endl;
		return false;
	}
	//read the first line
	getline(myFile, text);
	split(text, cDelim, vs);
	for (size_t i = 3; i < vs.size(); i++) {
		i_iter = find(viSampleNamesList.begin(), viSampleNamesList.end(), vs.at(i));
		if (i_iter == viSampleNamesList.end()) {
			viSampleNamesList.push_back(vs.at(i));
			viGroupSizeList.push_back(1);
		} else {
			(viGroupSizeList.at(i_iter - viSampleNamesList.begin()))++;
			iNumSamples++;
		}
	}
	iNumGroups = viSampleNamesList.size();
	myFile.close();
	vector<UINT32> vi;
	for (size_t i = 0; i < _vsInputFileList.size(); i++) {
		myFile.open(_vsInputFileList.at(i).c_str());
		if (!myFile.is_open()) {
			cout << "Unable to open file: " << _vsInputFileList.at(i) << endl;
			return false;
		}
		//check the group information
		getline(myFile, text);
		split(text, cDelim, vs);
		vi.clear();
		for (size_t i = 3; i < vs.size(); i++) {
			i_iter = find(viSampleNamesList.begin(), viSampleNamesList.end(), vs.at(i));
			if (i_iter == viSampleNamesList.end()) {
				vi.at(i_iter - viSampleNamesList.begin())++;}
			else {
				cout << "Group info is not consistent in file: " << _vsInputFileList.at(i) << endl;
				myFile.close();
				return false;
			}
		}
		if (vi != viGroupSizeList) {
			cout << "Group info is not consistent in file: " << _vsInputFileList.at(i) << endl;
			myFile.close();
			return false;
		}
		//check the variant information
		while (!myFile.eof()) {
			getline(myFile, text);
			split(text, cDelim, vs);
			if (vs.size() - 3 != iNumSamples) {
				cout << "Variant info is not consistent in file: " << _vsInputFileList.at(i) << endl;
				cout << "Error line: " << text << endl;
				myFile.close();
				return false;
			}
			iNumVariants++;
		}
		myFile.close();
	}

	return true;
}

bool GwasData::variantControl(vector<string> & _vsList){

	return true;
}
