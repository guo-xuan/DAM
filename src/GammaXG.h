/*
 * GammaXG.h
 *
 *  Created on: Dec 13, 2016
 *      Author: naux
 */

#ifndef GAMMAXG_H_
#define GAMMAXG_H_

#include <math.h>
#include <climits>
#include <limits>

#include "ContinuedFraction.h"

class GammaXG {
public:
	GammaXG();
	~GammaXG();

	static const double GAMMA;

	/**
	 * The value of the {@code g} constant in the Lanczos approximation, see
	 * {@link #lanczos(double)}.
	 * @since 3.1
	 */
	static const double LANCZOS_G;

	/**
	 * <p>
	 * Returns the value of log&nbsp;&Gamma;(x) for x&nbsp;&gt;&nbsp;0.
	 * </p>
	 * <p>
	 * For x &le; 8, the implementation is based on the double precision
	 * implementation in the <em>NSWC Library of Mathematics Subroutines</em>,
	 * {@code DGAMLN}. For x &gt; 8, the implementation is based on
	 * </p>
	 * <ul>
	 * <li><a href="http://mathworld.wolfram.com/GammaFunction.html">Gamma
	 *     Function</a>, equation (28).</li>
	 * <li><a href="http://mathworld.wolfram.com/LanczosApproximation.html">
	 *     Lanczos Approximation</a>, equations (1) through (5).</li>
	 * <li><a href="http://my.fit.edu/~gabdo/gamma.txt">Paul Godfrey, A note on
	 *     the computation of the convergent Lanczos complex Gamma
	 *     approximation</a></li>
	 * </ul>
	 *
	 * @param x Argument.
	 * @return the value of {@code log(Gamma(x))}, {@code Double.NaN} if
	 * {@code x <= 0.0}.
	 */
	static double logGamma(double x);

	/**
	 * Returns the regularized gamma function P(a, x).
	 *
	 * @param a Parameter.
	 * @param x Value.
	 * @return the regularized gamma function P(a, x).
	 * @throws MaxCountExceededException if the algorithm fails to converge.
	 */
	static double regularizedGammaP(double a, double x);

	/**
	 * Returns the regularized gamma function P(a, x).
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/RegularizedGammaFunction.html">
	 *   Regularized Gamma Function</a>, equation (1)
	 *  </li>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/IncompleteGammaFunction.html">
	 *   Incomplete Gamma Function</a>, equation (4).
	 *  </li>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html">
	 *   Confluent Hypergeometric Function of the First Kind</a>, equation (1).
	 *  </li>
	 * </ul>
	 *
	 * @param a the a parameter.
	 * @param x the value.
	 * @param epsilon When the absolute value of the nth item in the
	 * series is less than epsilon the approximation ceases to calculate
	 * further elements in the series.
	 * @param maxIterations Maximum number of "iterations" to complete.
	 * @return the regularized gamma function P(a, x)
	 * @throws MaxCountExceededException if the algorithm fails to converge.
	 */
	static double regularizedGammaP(double a, double x, double epsilon, int maxIterations);

	/**
	 * Returns the regularized gamma function Q(a, x) = 1 - P(a, x).
	 *
	 * @param a the a parameter.
	 * @param x the value.
	 * @return the regularized gamma function Q(a, x)
	 * @throws MaxCountExceededException if the algorithm fails to converge.
	 */
	static double regularizedGammaQ(double a, double x);

	/**
	 * Returns the regularized gamma function Q(a, x) = 1 - P(a, x).
	 *
	 * The implementation of this method is based on:
	 * <ul>
	 *  <li>
	 *   <a href="http://mathworld.wolfram.com/RegularizedGammaFunction.html">
	 *   Regularized Gamma Function</a>, equation (1).
	 *  </li>
	 *  <li>
	 *   <a href="http://functions.wolfram.com/GammaBetaErf/GammaRegularized/10/0003/">
	 *   Regularized incomplete gamma function: Continued fraction representations
	 *   (formula 06.08.10.0003)</a>
	 *  </li>
	 * </ul>
	 *
	 * @param a the a parameter.
	 * @param x the value.
	 * @param epsilon When the absolute value of the nth item in the
	 * series is less than epsilon the approximation ceases to calculate
	 * further elements in the series.
	 * @param maxIterations Maximum number of "iterations" to complete.
	 * @return the regularized gamma function P(a, x)
	 * @throws MaxCountExceededException if the algorithm fails to converge.
	 */
	static double regularizedGammaQ(const double a, double x, double epsilon, int maxIterations);

private:
	/** Maximum allowed numerical error. */
	static const double DEFAULT_EPSILON;

	/** Lanczos coefficients */
	static const double LANCZOS[];

	/** Avoid repeated computation of log of 2 PI in logGamma */
	static const double HALF_LOG_2_PI;

	/** The constant value of &radic;(2&pi;). */
	static const double SQRT_TWO_PI;

	// limits for switching algorithm in digamma
	/** C limit. */
	static const double C_LIMIT;

	/** S limit. */
	static const double S_LIMIT;

	/*
	 * Constants for the computation of double invGamma1pm1(double).
	 * Copied from DGAM1 in the NSWC library.
	 */

	/** The constant {@code A0} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_A0;

	/** The constant {@code A1} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_A1;

	/** The constant {@code B1} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B1;

	/** The constant {@code B2} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B2;

	/** The constant {@code B3} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B3;

	/** The constant {@code B4} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B4;

	/** The constant {@code B5} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B5;

	/** The constant {@code B6} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B6;

	/** The constant {@code B7} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B7;

	/** The constant {@code B8} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_B8;

	/** The constant {@code P0} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P0;

	/** The constant {@code P1} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P1;

	/** The constant {@code P2} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P2;

	/** The constant {@code P3} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P3;

	/** The constant {@code P4} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P4;

	/** The constant {@code P5} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P5;

	/** The constant {@code P6} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_P6;

	/** The constant {@code Q1} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_Q1;

	/** The constant {@code Q2} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_Q2;

	/** The constant {@code Q3} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_Q3;

	/** The constant {@code Q4} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_Q4;

	/** The constant {@code C} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C;

	/** The constant {@code C0} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C0;

	/** The constant {@code C1} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C1;

	/** The constant {@code C2} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C2;

	/** The constant {@code C3} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C3;

	/** The constant {@code C4} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C4;

	/** The constant {@code C5} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C5;

	/** The constant {@code C6} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C6;

	/** The constant {@code C7} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C7;

	/** The constant {@code C8} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C8;

	/** The constant {@code C9} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C9;

	/** The constant {@code C10} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C10;

	/** The constant {@code C11} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C11;

	/** The constant {@code C12} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C12;

	/** The constant {@code C13} defined in {@code DGAM1}. */
	static const double INV_GAMMA1P_M1_C13;
};

#endif /* GAMMAXG_H_ */
