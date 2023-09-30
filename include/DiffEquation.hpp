#ifndef DIFF_EQUATION_H
#define DIFF_EQUATION_H


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>




//--------------------------------------------------DiffEquation-------------------------------------------------------------------

/**
 * @brief class DiffEquation - class for represent a second order differential equation:
 *                              ..   .
 *                             Ax + Bx + Cx + D = 0 
 *                             as a system of two first order equations:
 *                              _  .
 *                             |   x = u
 *                            <     .
 *                             |_  Au + Bu + Cx = 0
 */
class DiffEquation
{
protected:
	double A_, B_, C_, D_;

public:
	DiffEquation(double A, double B, double C, double D) : A_(A), B_(B), C_(C), D_(D) {};
	// Frequency
	double W() const { return sqrt(C_ / A_); };
	// Attenuation
	double G() const { return B_ / A_; };

};

//------------------------------------------------HarmonicEquation----------------------------------------------------------------

/**
 * @brief class DiffEquation - equation of harmonic oscillations with frequency W.
 *  
 */
class HarmonicEquation : public DiffEquation
{
public:
	HarmonicEquation(double W) : DiffEquation(1.0, 0, W * W, 0) {};

};

#endif // DIFF_EQUATION_H
