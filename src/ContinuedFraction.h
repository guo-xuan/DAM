/*
 * ContinuedFraction.h
 *
 *  Created on: Dec 13, 2016
 *      Author: naux
 */

#ifndef CONTINUEDFRACTION_H_
#define CONTINUEDFRACTION_H_

#include <climits>
#include <limits>
#include <math.h>

class ContinuedFraction {
public:
	ContinuedFraction();
	virtual ~ContinuedFraction();

	/**
	 * Evaluates the continued fraction at the value x.
	 * <p>
	 * The implementation of this method is based on the modified Lentz algorithm as described
	 * on page 18 ff. in:
	 * <ul>
	 *   <li>
	 *   I. J. Thompson,  A. R. Barnett. "Coulomb and Bessel Functions of Complex Arguments and Order."
	 *   <a target="_blank" href="http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf">
	 *   http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf</a>
	 *   </li>
	 * </ul>
	 * <b>Note:</b> the implementation uses the terms a<sub>i</sub> and b<sub>i</sub> as defined in
	 * <a href="http://mathworld.wolfram.com/ContinuedFraction.html">Continued Fraction @ MathWorld</a>.
	 * </p>
	 *
	 * @param x the evaluation point.
	 * @param epsilon maximum error allowed.
	 * @param maxIterations maximum number of convergents
	 * @return the value of the continued fraction evaluated at x.
	 * @throws ConvergenceException if the algorithm fails to converge.
	 * @throws MaxCountExceededException if maximal number of iterations is reached
	 */
	double evaluate(double x, double epsilon, int maxIterations, double a_int);

private:
	/** Maximum allowed numerical error. */
	static const double DEFAULT_EPSILON;

protected:
	/**
	 * Access the n-th a coefficient of the continued fraction.  Since a can be
	 * a function of the evaluation point, x, that is passed in as well.
	 * @param n the coefficient index to retrieve.
	 * @param x the evaluation point.
	 * @return the n-th a coefficient.
	 */
	virtual double getA(int n, double x, double a_int) = 0;

	/**
	 * Access the n-th b coefficient of the continued fraction.  Since b can be
	 * a function of the evaluation point, x, that is passed in as well.
	 * @param n the coefficient index to retrieve.
	 * @param x the evaluation point.
	 * @return the n-th b coefficient.
	 */
	virtual double getB(int n, double x, double a_int) = 0;
};

#endif /* CONTINUEDFRACTION_H_ */
