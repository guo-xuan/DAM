/*
 * HyperGroup.h
 *
 *  Created on: Jun 23, 2016
 *      Author: naux
 */

#ifndef HYPERGROUP_H_
#define HYPERGROUP_H_

#include "Config.h"

class HyperGroup {
public:
	HyperGroup();
	~HyperGroup();

	static int BellNumber(int _n);
	static int BinomialCoefficient(int _n, int _k);
	static vector<vector<vector<int>>> getHyperGroup(int _n);
};

#endif /* HYPERGROUP_H_ */
