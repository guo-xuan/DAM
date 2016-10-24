/*
 * Variant.cpp
 *
 *  Created on: Sep 4, 2016
 *      Author: naux
 */

#include "Variant.h"

UINT32 Variant::iMissingData = -1;
UINT32 Variant::iLogMaxInteger = log(pow(2, sizeof(UINT32) * 8));
static UINT32 ** generate_data() {
	UINT32 ** miPower = new UINT32*[20];
	for(size_t i = 0;i < 20;i++) {
		miPower[i] = new UINT32[32];
		for(size_t j = 0;j < 32;j++) {
			miPower[i][j] = pow(i, j);
		}
	}
	return miPower;
}
UINT32 ** Variant::miPower = generate_data();
vector<UINT32> Variant::viNumSamplePerGroup;

Variant::Variant() {
	iNumTypes = 0;
	iNumSamplesPerInteger = 0;
	iVariantPosition = 0;
}

Variant::~Variant() {

}

void Variant::getDataInGroup(UINT32 _iGroupIndex, vector<UINT32> & _vi) {
	_vi.clear();
	UINT32 iCount = 0, iValue;
	for(size_t i = 0;i < vviData.at(_iGroupIndex).size();i++) {
		iValue = vviData.at(_iGroupIndex).at(i);
		for(size_t j = 0;j < iNumSamplesPerInteger;j++) {
			_vi.push_back(iValue % iNumTypes);
			iValue /= iNumTypes;
			iCount++;
			if(iCount >= viNumSamplePerGroup.at(_iGroupIndex)) {
				return;
			}
		}
	}
}

vector<vector<vector<double>>>* Variant::getHyperGroupFrequency() {
	return (& vvviHyperGroupFrequency);
}

vector<vector<vector<UINT32>>>* Variant::getHyperGroupOccurrence() {
	return (& vvviHyperGroupOccurrence);
}

bool Variant::operator <(const Variant & _o) const {
	if(iChromsomeNameIndex < _o.iChromsomeNameIndex) {
		return true;
	} else if(iChromsomeNameIndex > _o.iChromsomeNameIndex) {
		return false;
	} else {
		return (iVariantPosition < _o.iVariantPosition);
	}
}

void Variant::setData(const string & _sLine, char _delim, const vector<UINT32> & _viNumSamplePerGroup,
		const vector<UINT32> & _viSampleGroupInfo, const map<string, UINT32> & msiChromsomeIndex) {
	vector<UINT32> _vector;
	map<UINT32, UINT32> _mii;
	vector<UINT32> _viGroupCount;
	setData(_sLine, _delim, _viNumSamplePerGroup, _viSampleGroupInfo, _vector, _mii, _viGroupCount, msiChromsomeIndex);
}

void Variant::setData(const string & _sLine, char _delim, const vector<UINT32> & _viNumSamplePerGroup,
		const vector<UINT32> & _viSampleGroupInfo, vector<UINT32> & _vector, map<UINT32, UINT32> & _mii,
		vector<UINT32> & _viGroupCount, const map<string, UINT32> & msiChromsomeIndex) {
	_vector.clear();
	_mii.clear();
	_viGroupCount.resize(_viNumSamplePerGroup.size());
	fill(_viGroupCount.begin(), _viGroupCount.end(), 0);

	size_t i = 0;
	size_t pos = _sLine.find(_delim);
	sVariantName = _sLine.substr(i, pos - i);
	i = ++pos;
	pos = _sLine.find(_delim, pos);
	iChromsomeNameIndex = msiChromsomeIndex.at(_sLine.substr(i, pos - i));
	i = ++pos;
	pos = _sLine.find(_delim, pos);
	iVariantPosition = atoi(_sLine.substr(i, pos - i).c_str());
	_vector.clear();
	int iValue = 0;
	int iTemp = 0;
	i = ++pos;
	pos = _sLine.find(_delim, pos);
	while(pos != string::npos) {
		iValue = atoi(_sLine.substr(i, pos - i).c_str());
		_vector.push_back(iValue);
		if(_mii.find(iValue) == _mii.end()) {
			_mii[iValue] = 1;
		} else {
			_mii[iValue] += 1;
		}
		i = ++pos;
		pos = _sLine.find(_delim, pos);

		if(pos == string::npos) {
			iValue = atoi(_sLine.substr(i, _sLine.length() - i).c_str());
			_vector.push_back(iValue);
			if(_mii.find(iValue) == _mii.end()) {
				_mii[iValue] = 1;
			} else {
				_mii[iValue] += 1;
			}
		}
	}

	iNumTypes = _mii.size();
	iNumSamplesPerInteger = getNumSamplesPerInteger(iNumTypes);
	for(size_t i = 0;i < _viNumSamplePerGroup.size();i++) {
		iTemp = ceil(((double) _viNumSamplePerGroup.at(i)) / ((double) iNumSamplesPerInteger));
		vviData.push_back(vector<UINT32>(iTemp, 0));
		vviOccurrence.push_back(vector<UINT32>(iNumTypes, 0));
	}

	UINT32 iArrayIndex = 0;
	UINT32 iValueIndex = 0;

	for(UINT32 k = 0;k < _vector.size();k++) {
		if(_vector.at(k) == iMissingData) {
			(_viGroupCount.at(_viSampleGroupInfo.at(k))) += 1;
			continue;
		}
		iArrayIndex = _viGroupCount.at(_viSampleGroupInfo.at(k));
		iValueIndex = iArrayIndex % iNumSamplesPerInteger;
		iArrayIndex = iArrayIndex / iNumSamplesPerInteger;
		vviData.at(_viSampleGroupInfo.at(k)).at(iArrayIndex) += (_vector.at(k) * miPower[iNumTypes][iValueIndex]);
		vviOccurrence.at(_viSampleGroupInfo.at(k)).at(_vector.at(k)) += 1;
		(_viGroupCount.at(_viSampleGroupInfo.at(k))) += 1;
	}
}

void Variant::setFrequency(const vector<vector<vector<int>>>& _vvviHyperGroupInfo) {
	double sum = 0;
	for(size_t i=0;i<_vvviHyperGroupInfo.size();i++) {
		vvviHyperGroupOccurrence.push_back(vector<vector<UINT32>>());
		vvviHyperGroupFrequency.push_back(vector<vector<double>>());
		for(size_t j=0;j<_vvviHyperGroupInfo.at(i).size();j++) {
			vvviHyperGroupOccurrence.at(i).push_back(vector<UINT32>());
			vvviHyperGroupFrequency.at(i).push_back(vector<double>());
			// k is the group index
			vvviHyperGroupOccurrence.at(i).at(j).resize(iNumTypes, 0);
			vvviHyperGroupFrequency.at(i).at(j).resize(iNumTypes, 0.0);
			for(size_t k=0;k<_vvviHyperGroupInfo.at(i).at(j).size();k++) {
				for(size_t l=0;l<iNumTypes;l++) {
					vvviHyperGroupOccurrence.at(i).at(j).at(l) += vviOccurrence.at(_vvviHyperGroupInfo.at(i).at(j).at(k)).at(l);
				}
			}
			sum = 0;
			for(size_t l=0;l<iNumTypes;l++) {
				sum += vvviHyperGroupOccurrence.at(i).at(j).at(l);
			}
			for(size_t l=0;l<iNumTypes;l++) {
				vvviHyperGroupFrequency.at(i).at(j).at(l) = ((double)vvviHyperGroupOccurrence.at(i).at(j).at(l))/sum;
			}
		}
	}
}

//---------------------------------------------PRIVATE--------------------------------------------------------------------------------------------

UINT32 Variant::getNumSamplesPerInteger(double base) {
	return (UINT32) (iLogMaxInteger / log(base));
}

