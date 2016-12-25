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
		ret = NAN;
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

double GammaXG::regularizedGammaP(double a, double x, double epsilon, long maxIterations) {
	double ret = 0;

	if(isnan(a) || isnan(x) || (a <= 0.0) || (x < 0.0)) {
		ret = NAN;
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
		while(fabs(an / sum) > epsilon && n < maxIterations && sum < INFINITY) {
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

double GammaXG::regularizedGammaQ(const double a, double x, double epsilon, long maxIterations) {
	double ret;

	if(isnan(a) || isnan(x) || (a <= 0.0) || (x < 0.0)) {
		ret = NAN;
	} else if(x == 0.0) {
		ret = 1.0;
	} else if(x < a + 1.0) {
		// use regularizedGammaP because it should converge faster in this
		// case.
		ret = 1.0 - regularizedGammaP(a, x, epsilon, maxIterations);
	} else {
		// create continued fraction
		class : public ContinuedFraction {
		protected:
			double getA(int n, double x, double a_int) {
				return ((2.0 * n) + 1.0) - a_int + x;
			};

			double getB(int n, double x, double a_int) {
				return n * (a_int - n);
			};
		} cf;

		ret = 1.0 / cf.evaluate(x, epsilon, maxIterations, a);
		double loggammaa = -logGamma(a);
		double summm = (a * log(x));
		summm = summm - x - loggammaa;
		summm = exp(summm);
		summm = summm * ret;
		ret = exp(-x + (a * log(x)) - logGamma(a)) * ret;
	}

	return ret;
}

double GammaXG::digamma(double x) {
	if(x > 0 && x <= S_LIMIT) {
		// use method 5 from Bernardo AS103
		// accurate to O(x)
		return -GAMMA - 1 / x;
	}

	if(x >= C_LIMIT) {
		// use method 4 (accurate to O(1/x^8)
		double inv = 1 / (x * x);
		//            1       1        1         1
		// log(x) -  --- - ------ + ------- - -------
		//           2 x   12 x^2   120 x^4   252 x^6
		return log(x) - 0.5 / x - inv * ((1.0 / 12) + inv * (1.0 / 120 - inv / 252));
	}

	return digamma(x + 1) - 1 / x;
}

double GammaXG::trigamma(double x) {
	if(x > 0 && x <= S_LIMIT) {
		return 1 / (x * x);
	}

	if(x >= C_LIMIT) {
		double inv = 1 / (x * x);
		//  1    1      1       1       1
		//  - + ---- + ---- - ----- + -----
		//  x      2      3       5       7
		//      2 x    6 x    30 x    42 x
		return 1 / x + inv / 2 + inv / x * (1.0 / 6 - inv * (1.0 / 30 + inv / 42));
	}

	return trigamma(x + 1) + 1 / (x * x);
}

double GammaXG::lanczos(const double x) {
	double sum = 0.0;
	for(int i = (sizeof(LANCZOS) / sizeof(*LANCZOS)) - 1;i > 0;--i) {
		sum = sum + (LANCZOS[i] / (x + i));
	}
	return sum + LANCZOS[0];
}

double GammaXG::invGamma1pm1(const double x) {

	double ret;
	const double t = x <= 0.5 ? x : (x - 0.5) - 0.5;
	if(t < 0.0) {
		const double a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
		double b = INV_GAMMA1P_M1_B8;
		b = INV_GAMMA1P_M1_B7 + t * b;
		b = INV_GAMMA1P_M1_B6 + t * b;
		b = INV_GAMMA1P_M1_B5 + t * b;
		b = INV_GAMMA1P_M1_B4 + t * b;
		b = INV_GAMMA1P_M1_B3 + t * b;
		b = INV_GAMMA1P_M1_B2 + t * b;
		b = INV_GAMMA1P_M1_B1 + t * b;
		b = 1.0 + t * b;

		double c = INV_GAMMA1P_M1_C13 + t * (a / b);
		c = INV_GAMMA1P_M1_C12 + t * c;
		c = INV_GAMMA1P_M1_C11 + t * c;
		c = INV_GAMMA1P_M1_C10 + t * c;
		c = INV_GAMMA1P_M1_C9 + t * c;
		c = INV_GAMMA1P_M1_C8 + t * c;
		c = INV_GAMMA1P_M1_C7 + t * c;
		c = INV_GAMMA1P_M1_C6 + t * c;
		c = INV_GAMMA1P_M1_C5 + t * c;
		c = INV_GAMMA1P_M1_C4 + t * c;
		c = INV_GAMMA1P_M1_C3 + t * c;
		c = INV_GAMMA1P_M1_C2 + t * c;
		c = INV_GAMMA1P_M1_C1 + t * c;
		c = INV_GAMMA1P_M1_C + t * c;
		if(x > 0.5) {
			ret = t * c / x;
		} else {
			ret = x * ((c + 0.5) + 0.5);
		}
	} else {
		double p = INV_GAMMA1P_M1_P6;
		p = INV_GAMMA1P_M1_P5 + t * p;
		p = INV_GAMMA1P_M1_P4 + t * p;
		p = INV_GAMMA1P_M1_P3 + t * p;
		p = INV_GAMMA1P_M1_P2 + t * p;
		p = INV_GAMMA1P_M1_P1 + t * p;
		p = INV_GAMMA1P_M1_P0 + t * p;

		double q = INV_GAMMA1P_M1_Q4;
		q = INV_GAMMA1P_M1_Q3 + t * q;
		q = INV_GAMMA1P_M1_Q2 + t * q;
		q = INV_GAMMA1P_M1_Q1 + t * q;
		q = 1.0 + t * q;

		double c = INV_GAMMA1P_M1_C13 + (p / q) * t;
		c = INV_GAMMA1P_M1_C12 + t * c;
		c = INV_GAMMA1P_M1_C11 + t * c;
		c = INV_GAMMA1P_M1_C10 + t * c;
		c = INV_GAMMA1P_M1_C9 + t * c;
		c = INV_GAMMA1P_M1_C8 + t * c;
		c = INV_GAMMA1P_M1_C7 + t * c;
		c = INV_GAMMA1P_M1_C6 + t * c;
		c = INV_GAMMA1P_M1_C5 + t * c;
		c = INV_GAMMA1P_M1_C4 + t * c;
		c = INV_GAMMA1P_M1_C3 + t * c;
		c = INV_GAMMA1P_M1_C2 + t * c;
		c = INV_GAMMA1P_M1_C1 + t * c;
		c = INV_GAMMA1P_M1_C0 + t * c;

		if(x > 0.5) {
			ret = (t / x) * ((c - 0.5) - 0.5);
		} else {
			ret = x * c;
		}
	}

	return ret;
}

double GammaXG::logGamma1p(const double x) {
	return -log1p(invGamma1pm1(x));
}

double GammaXG::gamma(const double x) {

	if((x == rint(x)) && (x <= 0.0)) {
		return NAN;
	}

	double ret;
	const double absX = fabs(x);
	if(absX <= 20.0) {
		if(x >= 1.0) {
			/*
			 * From the recurrence relation
			 * Gamma(x) = (x - 1) * ... * (x - n) * Gamma(x - n),
			 * then
			 * Gamma(t) = 1 / [1 + invGamma1pm1(t - 1)],
			 * where t = x - n. This means that t must satisfy
			 * -0.5 <= t - 1 <= 1.5.
			 */
			double prod = 1.0;
			double t = x;
			while(t > 2.5) {
				t = t - 1.0;
				prod *= t;
			}
			ret = prod / (1.0 + invGamma1pm1(t - 1.0));
		} else {
			/*
			 * From the recurrence relation
			 * Gamma(x) = Gamma(x + n + 1) / [x * (x + 1) * ... * (x + n)]
			 * then
			 * Gamma(x + n + 1) = 1 / [1 + invGamma1pm1(x + n)],
			 * which requires -0.5 <= x + n <= 1.5.
			 */
			double prod = x;
			double t = x;
			while(t < -0.5) {
				t = t + 1.0;
				prod *= t;
			}
			ret = 1.0 / (prod * (1.0 + invGamma1pm1(t)));
		}
	} else {
		const double y = absX + LANCZOS_G + 0.5;
		const double gammaAbs = SQRT_TWO_PI / x * pow(y, absX + 0.5) * exp(-y) * lanczos(absX);
		if(x > 0.0) {
			ret = gammaAbs;
		} else {
			/*
			 * From the reflection formula
			 * Gamma(x) * Gamma(1 - x) * sin(pi * x) = pi,
			 * and the recurrence relation
			 * Gamma(1 - x) = -x * Gamma(-x),
			 * it is found
			 * Gamma(x) = -pi / [x * sin(pi * x) * Gamma(-x)].
			 */
			ret = -M_PI / (x * sin(M_PI * x) * gammaAbs);
		}
	}
	return ret;
}
