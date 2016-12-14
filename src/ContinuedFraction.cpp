/*
 * ContinuedFraction.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: naux
 */

#include "ContinuedFraction.h"
const double ContinuedFraction::DEFAULT_EPSILON = 10e-9;

ContinuedFraction::ContinuedFraction() {
	// TODO Auto-generated constructor stub

}

ContinuedFraction::~ContinuedFraction() {
	// TODO Auto-generated destructor stub
}

double ContinuedFraction::evaluate(double x, double epsilon, int maxIterations, double a_int) {
	const double small = 1e-50;
	double hPrev = getA(0, x, a_int);

	// use the value of small as epsilon criteria for zero checks
	if(fabs(hPrev) <= small) {
		hPrev = small;
	}

	int n = 1;
	double dPrev = 0.0;
	double cPrev = hPrev;
	double hN = hPrev;

	while(n < maxIterations) {
		const double a = getA(n, x, a_int);
		const double b = getB(n, x, a_int);

		double dN = a + b * dPrev;
		if(fabs(dN) <= small) {
			dN = small;
		}
		double cN = a + b / cPrev;
		if(fabs(cN) <= small) {
			cN = small;
		}

		dN = 1 / dN;
		const double deltaN = cN * dN;
		hN = hPrev * deltaN;

		if(fabs(deltaN - 1.0) < epsilon) {
			break;
		}

		dPrev = dN;
		cPrev = cN;
		hPrev = hN;
		n++;
	}
	return hN;
}
