/*
 * GwasData.cpp
 *
 *  Created on: May 3, 2016
 *      Author: naux
 */

#include "GwasData.h"

char GwasData::cDelim = '\t';
UINT32 GwasData::iMissingData = -1;
UINT32 GwasData::iLogMaxInteger = log(pow(2, sizeof(UINT32) * 8));
UINT32 GwasData::iChunkSize = 5000;

GwasData::GwasData() {
	iNumVariants = 0;
	iNumGroups = 0;
	iNumTotalSamplesAcrossGroups = 0;
	miPower = new UINT32*[20];
	for (size_t i = 0; i < 20; i++) {
		miPower[i] = new UINT32[32];
		for (size_t j = 0; j < 32; j++) {
			miPower[i][j] = pow(i, j);
		}
	}
}

GwasData::~GwasData() {
	for (size_t i = 0; i < this->iNumGroups; i++) {
		for (size_t j = 0; j < this->iNumVariants; j++) {
			delete[] this->vmDataMatrix.at(i)[j];
		}
		delete[] this->vmDataMatrix.at(i);
	}
	for (size_t i = 0; i < 20; i++) {
		delete[] miPower[i];
	}
	delete[] miPower;

	for (size_t i = 0; i < vviVariantFrequence.size(); i++) {
		for (size_t j = 0; j < vviVariantFrequence.at(i)->size(); j++) {
			delete[] vviVariantFrequence.at(i)->at(j);
		}
		vviVariantFrequence.at(i)->clear();
		delete vviVariantFrequence.at(i);
	}

	for (size_t i = 0; i < viSampleVectorSizePerGroup.size(); i++) {
		delete[] viSampleVectorSizePerGroup.at(i);
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

UINT32 getNumSamplesPer32Bits(double base) {
	return (UINT32) (GwasData::iLogMaxInteger / log(base));
}

void GwasData::getValue(UINT32 _iVariantIndex, UINT32 _iGroupIndex, vector<UINT32> & _vi) {
	_vi.clear();
	UINT32 iCount = 0, temp = viNumSamples32Bits.at(_iVariantIndex), value = 0, iMark = viNumVariantTypes.at(
			_iVariantIndex);
	for (size_t i = 0; i < viSampleVectorSizePerGroup.at(_iVariantIndex)[_iGroupIndex]; i++) {
		value = vmDataMatrix.at(_iGroupIndex)[_iVariantIndex][i];
		for (size_t j = 0; j < temp; j++) {
			_vi.push_back(value % iMark);
			value /= iMark;
			iCount++;
			if (iCount >= viGroupSizeList.at(_iGroupIndex)) {
				return;
			}
		}
	}
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

bool GwasData::loadBasicInfo(vector<string> & _vsInputFileList) {
	cout << "Initializing Space" << endl;
	ifstream myFile;
	string text;
	string element;
	vector<string> vs;
	vector<UINT32> vi;
	vector<UINT32> vi2;
	vector<string>::const_iterator vs_iter;
	vector<UINT32>::const_iterator vi_iter;
	//open the first file
	myFile.open(_vsInputFileList.at(0).c_str());
	if (!myFile.is_open()) {
		cout << "Unable to open file: " << _vsInputFileList.at(0) << endl;
		return false;
	}
	//read the first line
	getline(myFile, text);
	split(text, cDelim, vs);
	for (size_t i = 3; i < vs.size(); i++) {
		vs_iter = find(viGroupNamesSet.begin(), viGroupNamesSet.end(), vs.at(i));
		if (vs_iter == viGroupNamesSet.end()) {
			viGroupNamesSet.push_back(vs.at(i));
			viGroupSizeList.push_back(1);
			iNumTotalSamplesAcrossGroups++;
		} else {
			(viGroupSizeList.at(vs_iter - viGroupNamesSet.begin()))++;
			iNumTotalSamplesAcrossGroups++;
		}
	}
	iNumGroups = viGroupNamesSet.size();
	myFile.close();
	//check the group info and the number of variant types per variant for all input files
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
		vi.resize(iNumGroups, 0);
		for (size_t i = 3; i < vs.size(); i++) {
			vs_iter = find(viGroupNamesSet.begin(), viGroupNamesSet.end(), vs.at(i));
			if (vs_iter == viGroupNamesSet.end()) {
				cout << "Group info is not consistent in file: " << _vsInputFileList.at(i) << endl;
				myFile.close();
				return false;
			} else {
				vi.at(vs_iter - viGroupNamesSet.begin())++;}
			}
		if (vi != viGroupSizeList) {
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
				vi2.clear();
				vs_iter = find(vsChromsomeNamesSet.begin(), vsChromsomeNamesSet.end(), vs.at(1));
				if (vs_iter == vsChromsomeNamesSet.end()) {
					vsChromsomeNamesSet.push_back(vs.at(1));
				}
				for (size_t j = 0; j < vi.size(); j++) {
					if (vi.at(j) == iMissingData) {
						continue;
					}
					vi_iter = find(vi2.begin(), vi2.end(), vi.at(j));
					if (vi_iter == vi2.end()) {
						vi2.push_back(vi.at(j));
					}
				}
				viNumVariantTypes.push_back(vi2.size());
				iNumVariants++;
			}

		}
		myFile.close();
	}
	//ask for memory for the matrix
	double temp = 0;
	viNumSamples32Bits.resize(iNumVariants);
	viSampleVectorSizePerGroup.resize(iNumVariants);
	for (size_t i = 0; i < iNumVariants; i++) {
		temp = getNumSamplesPer32Bits(viNumVariantTypes.at(i));
		viNumSamples32Bits.at(i) = temp;
		viSampleVectorSizePerGroup.at(i) = new UINT16[iNumGroups]();
		for (size_t j = 0; j < iNumGroups; j++) {
			viSampleVectorSizePerGroup.at(i)[j] = ceil(((double) viGroupSizeList.at(j)) / temp);
		}
	}
	for (size_t i = 0; i < iNumGroups; i++) {
		UINT32** matrix = new UINT32*[iNumVariants];
		for (size_t j = 0; j < iNumVariants; j++) {
			matrix[j] = new UINT32[viSampleVectorSizePerGroup.at(j)[i]]();

		}
		vmDataMatrix.push_back(matrix);
	}

	viVariantChromsomeNameList.resize(iNumVariants);
	viVariantPositionList.resize(iNumVariants);
	vsVariantNameList.resize(iNumVariants);
	//frequency for each variant
	vviVariantFrequence.resize(iNumGroups, NULL);
	for (UINT32 i = 0; i < iNumGroups; i++) {
		vviVariantFrequence.at(i) = new vector<UINT32*>();
		vviVariantFrequence.at(i)->resize(iNumVariants, NULL);
		for (UINT32 j = 0; j < iNumVariants; j++) {
			vviVariantFrequence.at(i)->at(j) = new UINT32[viNumVariantTypes.at(j)]();
		}
	}
	return true;
}

bool GwasData::loadDataParallel(vector<string> & _vsInputFileList) {
	ifstream myFile;
	string text;
	UINT32 chunkCounter;
	vector<string>* lineData = new vector<string>();
	UINT32 iVariantIndex = 0;
	vector<string> vsGroupInfo;
	vector<UINT8> vsGroupIndex;
	vsGroupIndex.resize(iNumTotalSamplesAcrossGroups, 0);
	vector<string>::const_iterator vs_iterGroup;
	for (size_t i = 0; i < _vsInputFileList.size(); i++) {
		myFile.open(_vsInputFileList.at(i).c_str());
		if (!myFile.is_open()) {
			cout << "Unable to open file: " << _vsInputFileList.at(i) << endl;
			return false;
		}
		//get the header line
		getline(myFile, text);
		split(text, cDelim, vsGroupInfo);
		for (size_t j = 3; j < vsGroupInfo.size(); j++) {
			vs_iterGroup = find(viGroupNamesSet.begin(), viGroupNamesSet.end(), vsGroupInfo.at(j));
			vsGroupIndex.at(j - 3) = (vs_iterGroup - viGroupNamesSet.begin());
		}
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
#pragma omp parallel
			{
				vector<string> vs;
				vector<UINT32> vi;
				string sLine;
				vector<string>::const_iterator vs_iter;
				UINT32 iArrayIndex = 0;
				UINT32 iValueIndex = 0;
				vector<UINT32> viGroupCount(iNumGroups, 0);
#pragma omp for schedule(guided)
				for (UINT32 j = 0; j < chunkCounter; j++) {
					split(lineData->at(j), cDelim, 3, vs, vi);
					vsVariantNameList.at(iVariantIndex + j) = vs.at(0);
					vs_iter = find(vsChromsomeNamesSet.begin(), vsChromsomeNamesSet.end(), vs.at(1));
					viVariantChromsomeNameList.at(iVariantIndex + j) = (vs_iter - vsChromsomeNamesSet.begin());
					viVariantPositionList.at(iVariantIndex + j) = stoi(vs.at(2));
					for (UINT32 k = 0; k < vi.size(); k++) {
						iArrayIndex = viGroupCount.at(vsGroupIndex.at(k));
						iValueIndex = iArrayIndex % viNumSamples32Bits.at(iVariantIndex + j);
						iArrayIndex = iArrayIndex / viNumSamples32Bits.at(iVariantIndex + j);
						vmDataMatrix.at(vsGroupIndex.at(k))[iVariantIndex + j][iArrayIndex] += (vi.at(k)
								* miPower[viNumVariantTypes.at(iVariantIndex + j)][iValueIndex]);
						++(viGroupCount.at(vsGroupIndex.at(k)));
						++(vviVariantFrequence.at(vsGroupIndex.at(k))->at(iVariantIndex + j)[vi.at(k)]);
					}
					fill(viGroupCount.begin(), viGroupCount.end(), 0);
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

void GwasData::writeOutput(string & _sFilename) {
	ofstream filePointer;
	filePointer.open(_sFilename.c_str(), ios_base::out);
	vector<UINT32> vi;
	//header
	filePointer << "ID\tChr\tPos";
	for (size_t j = 0; j < iNumGroups; j++) {
		for (size_t i = 0; i < viGroupSizeList.at(j); i++) {
			filePointer << "\t" << j;
		}
	}
	filePointer << "\n";
	for (size_t i = 0; i < iNumVariants; i++) {
		filePointer << vsVariantNameList.at(i) << "\t";
		filePointer << vsChromsomeNamesSet.at(viVariantChromsomeNameList.at(i)) << "\t";
		filePointer << viVariantPositionList.at(i);
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
