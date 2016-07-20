/*
 * HyperGroup.cpp
 *
 *  Created on: Jun 23, 2016
 *      Author: naux
 */

#include "HyperGroup.h"
#include "Config.h"

HyperGroup::HyperGroup() {

}

HyperGroup::~HyperGroup() {

}

int HyperGroup::BellNumber(int _n) {
	int iBellNumber = 0;
	int * aiBellNumber = new int[_n + 1];
	aiBellNumber[0] = 1;
	for (int i = 1; i < _n + 1; i++) {
		aiBellNumber[i] = 0;
		for (int j = 0; j < i; j++) {
			aiBellNumber[i] += aiBellNumber[j] * BinomialCoefficient(i - 1, j);
		}
	}

	iBellNumber = aiBellNumber[_n];
	delete[] aiBellNumber;
	return iBellNumber;
}

int HyperGroup::BinomialCoefficient(int _n, int _k) {
	// Base Cases
	if (_k == 0 || _k == _n) {
		return 1;
	}

	// Recur
	return BinomialCoefficient(_n - 1, _k - 1) + BinomialCoefficient(_n - 1, _k);
}

void decoding(int _base, int * _aiEncoding, int _len, int _iNumber) {
	for (int i = 0; i < _len; i++) {
		_aiEncoding[_len - i - 1] = (_iNumber / ((int) pow(_base, i))) % (_base);
	}
}

int getNBaseNumber(int _base, int * _aiEncoding, int _len) {
	int iNumber = 0;
	for (int i = 0; i < _len; i++) {
		iNumber += _aiEncoding[_len - i - 1] * pow(_base, i);
	}
	return iNumber;
}

bool isDuplicate(int _n, int * _aiEncoding) {
	bool bGoodDigit = false;
	for (int i = 0; i < _n; i++) {
		if (_aiEncoding[i] == 0) {
			return true;
		}
		if (_aiEncoding[i] == 1) {
			continue;
		}
		bGoodDigit = false;
		for (int j = 0; j < i; j++) {
			if (_aiEncoding[i] - _aiEncoding[j] == 1) {
				bGoodDigit = true;
			}
		}
		if (!bGoodDigit) {
			return true;
		}
	}
	return false;
}

/**
 * idea is from
 * https://compprog.wordpress.com/2007/10/15/generating-the-partitions-of-a-set/
 * group index starts at 0
 */
vector<vector<vector<int>>> HyperGroup::getHyperGroup(int _n) {
	vector<vector<vector<int>>> vvviHyperGroup;
	int iBase = _n + 1;
	int aiEncoding[_n];
	for(int i=0;i<_n;i++) {
		//the subset index; all 1 means all groups in the same subset
		aiEncoding[i] = i+1;
	}
	//translate to _n-base number system
	int iNumber = 0, iUpperBound;
	iUpperBound = getNBaseNumber(iBase, aiEncoding, _n);
	for(int i=0;i<_n;i++) {
		//the subset index; all 1 means all groups in the same subset
		aiEncoding[i] = 1;
	}
	iNumber = getNBaseNumber(iBase, aiEncoding, _n);
	do {
		//decoding
		decoding(iBase, aiEncoding, _n, iNumber);
		//is a duplicate encoding
		if(!isDuplicate(_n, aiEncoding)) {
			vector<vector<int>> vviHyperGroup;
			for(int i=0;i<_n;i++) {
				if(aiEncoding[i]>(int)vviHyperGroup.size()) {
					vviHyperGroup.push_back(vector<int>());
				}
				vviHyperGroup.at(aiEncoding[i]-1).push_back(i);
			}
			vvviHyperGroup.push_back(vviHyperGroup);
		}
		iNumber++;
	}while(iNumber<=iUpperBound);

	return vvviHyperGroup;
}
