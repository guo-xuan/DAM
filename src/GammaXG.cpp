/*
 * GammaXG.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: naux
 */

#include "GammaXG.h"

const double GammaXG::GAMMA = 0.577215664901532860606512090082;
const double GammaXG::LANCZOS_G = 607.0 / 128.0;
const double GammaXG::DEFAULT_EPSILON = 10e-15;
const double GammaXG::LANCZOS[] = { 0.99999999999999709182, 57.156235665862923517, -59.597960355475491248,
		14.136097974741747174, -0.49191381609762019978, .33994649984811888699e-4, .46523628927048575665e-4,
		-.98374475304879564677e-4, .15808870322491248884e-3, -.21026444172410488319e-3, .21743961811521264320e-3,
		-.16431810653676389022e-3, .84418223983852743293e-4, -.26190838401581408670e-4, .36899182659531622704e-5, };
const double GammaXG::HALF_LOG_2_PI = 0.5 * log(2.0 * M_PI);
const double GammaXG::SQRT_TWO_PI = 2.506628274631000502;
const double GammaXG::C_LIMIT = 49;
const double GammaXG::S_LIMIT = 1e-5;
const double GammaXG::INV_GAMMA1P_M1_A0 = .611609510448141581788E-08;
const double GammaXG::INV_GAMMA1P_M1_A1 = .624730830116465516210E-08;
const double GammaXG::INV_GAMMA1P_M1_B1 = .203610414066806987300E+00;
const double GammaXG::INV_GAMMA1P_M1_B2 = .266205348428949217746E-01;
const double GammaXG::INV_GAMMA1P_M1_B3 = .493944979382446875238E-03;
const double GammaXG::INV_GAMMA1P_M1_B4 = -.851419432440314906588E-05;
const double GammaXG::INV_GAMMA1P_M1_B5 = -.643045481779353022248E-05;
const double GammaXG::INV_GAMMA1P_M1_B6 = .992641840672773722196E-06;
const double GammaXG::INV_GAMMA1P_M1_B7 = -.607761895722825260739E-07;
const double GammaXG::INV_GAMMA1P_M1_B8 = .195755836614639731882E-09;
const double GammaXG::INV_GAMMA1P_M1_P0 = .6116095104481415817861E-08;
const double GammaXG::INV_GAMMA1P_M1_P1 = .6871674113067198736152E-08;
const double GammaXG::INV_GAMMA1P_M1_P2 = .6820161668496170657918E-09;
const double GammaXG::INV_GAMMA1P_M1_P3 = .4686843322948848031080E-10;
const double GammaXG::INV_GAMMA1P_M1_P4 = .1572833027710446286995E-11;
const double GammaXG::INV_GAMMA1P_M1_P5 = -.1249441572276366213222E-12;
const double GammaXG::INV_GAMMA1P_M1_P6 = .4343529937408594255178E-14;
const double GammaXG::INV_GAMMA1P_M1_Q1 = .3056961078365221025009E+00;
const double GammaXG::INV_GAMMA1P_M1_Q2 = .5464213086042296536016E-01;
const double GammaXG::INV_GAMMA1P_M1_Q3 = .4956830093825887312020E-02;
const double GammaXG::INV_GAMMA1P_M1_Q4 = .2692369466186361192876E-03;
const double GammaXG::INV_GAMMA1P_M1_C = -.422784335098467139393487909917598E+00;
const double GammaXG::INV_GAMMA1P_M1_C0 = .577215664901532860606512090082402E+00;
const double GammaXG::INV_GAMMA1P_M1_C1 = -.655878071520253881077019515145390E+00;
const double GammaXG::INV_GAMMA1P_M1_C2 = -.420026350340952355290039348754298E-01;
const double GammaXG::INV_GAMMA1P_M1_C3 = .166538611382291489501700795102105E+00;
const double GammaXG::INV_GAMMA1P_M1_C4 = -.421977345555443367482083012891874E-01;
const double GammaXG::INV_GAMMA1P_M1_C5 = -.962197152787697356211492167234820E-02;
const double GammaXG::INV_GAMMA1P_M1_C6 = .721894324666309954239501034044657E-02;
const double GammaXG::INV_GAMMA1P_M1_C7 = -.116516759185906511211397108401839E-02;
const double GammaXG::INV_GAMMA1P_M1_C8 = -.215241674114950972815729963053648E-03;
const double GammaXG::INV_GAMMA1P_M1_C9 = .128050282388116186153198626328164E-03;
const double GammaXG::INV_GAMMA1P_M1_C10 = -.201348547807882386556893914210218E-04;
const double GammaXG::INV_GAMMA1P_M1_C11 = -.125049348214267065734535947383309E-05;
const double GammaXG::INV_GAMMA1P_M1_C12 = .113302723198169588237412962033074E-05;
const double GammaXG::INV_GAMMA1P_M1_C13 = -.205633841697760710345015413002057E-06;

GammaXG::GammaXG() {

}

GammaXG::~GammaXG() {

}

double GammaXG::logGamma(double x) {
	double ret;

	if(isnan(x) || (x <= 0.0)) {
		ret = nan;
	} else if(x < 0.5) {
		return logGamma1p(x) - log(x);
	} else if(x <= 2.5) {
		return logGamma1p((x - 0.5) - 0.5);
	} else if(x <= 8.0) {
		int n = (int) floor(x - 1.5);
		double prod = 1.0;
		for(int i = 1;i <= n;i++) {
			prod *= x - i;
		}
		return logGamma1p(x - (n + 1)) + log(prod);
	} else {
		double sum = lanczos(x);
		double tmp = x + LANCZOS_G + .5;
		ret = ((x + .5) * log(tmp)) - tmp + HALF_LOG_2_PI + log(sum / x);
	}

	return ret;
}

double GammaXG::regularizedGammaP(double a, double x) {
	return regularizedGammaP(a, x, DEFAULT_EPSILON, LONG_MAX);
}

double GammaXG::regularizedGammaP(double a, double x, double epsilon, int maxIterations) {
	double ret = 0;

	if(isnan(a) || isnan(x) || (a <= 0.0) || (x < 0.0)) {
		ret = nan;
	} else if(x == 0.0) {
		ret = 0.0;
	} else if(x >= a + 1) {
		// use regularizedGammaQ because it should converge faster in this
		// case.
		ret = 1.0 - regularizedGammaQ(a, x, epsilon, maxIterations);
	} else {
		// calculate series
		double n = 0.0; // current element index
		double an = 1.0 / a; // n-th element in the series
		double sum = an; // partial sum
		while(abs(an / sum) > epsilon && n < maxIterations && sum < INFINITY) {
			// compute next element in the series
			n = n + 1.0;
			an = an * (x / (a + n));

			// update partial sum
			sum = sum + an;
		}
		if(n >= maxIterations) {

		} else if(isinf(sum)) {
			ret = 1.0;
		} else {
			ret = exp(-x + (a * log(x)) - logGamma(a)) * sum;
		}
	}

	return ret;
}

double GammaXG::regularizedGammaQ(double a, double x) {
	return regularizedGammaQ(a, x, DEFAULT_EPSILON, LONG_MAX);
}

double GammaXG::regularizedGammaQ(const double a, double x, double epsilon, int maxIterations) {
	double ret;

	if(isnan(a) || isnan(x) || (a <= 0.0) || (x < 0.0)) {
		ret = nan;
	} else if(x == 0.0) {
		ret = 1.0;
	} else if(x < a + 1.0) {
		// use regularizedGammaP because it should converge faster in this
		// case.
		ret = 1.0 - regularizedGammaP(a, x, epsilon, maxIterations);
	} else {
		// create continued fraction
		ContinuedFraction cf = ContinuedFraction()
		{double getA(int n, double x) {
			return ((2.0 * n) + 1.0) - a + x;
		}

		double getB(int n, double x) {
			return n * (a - n);
		}
	};

	ret = 1.0 / cf.evaluate(x, epsilon, maxIterations);
	double loggammaa = -logGamma(a);
	double summm = (a * log(x));
	summm = summm - x - loggammaa;
	summm = exp(summm);
	summm = summm * ret;
	ret = exp(-x + (a * log(x)) - logGamma(a)) * ret;
}

return ret;
}
