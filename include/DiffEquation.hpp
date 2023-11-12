#ifndef DIFF_EQUATION_H
#define DIFF_EQUATION_H


#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <linalg.h>




template <typename T, unsigned Dim>
using Coordinates = linalg::vec<T, Dim>;

//--------------------------------------------------DiffEquation-------------------------------------------------------------------

/**
 * @brief class DiffEquation - class for represent Dim order differential equation. 
 *                             It is reduced to a system of equations using expressions for the highest derivative.
 *                             The DiffEquation can return a derivative vector in any state.
 *              
 */
template <typename T, unsigned Dim>
class DiffEquation
{
public:
	DiffEquation() {};
	virtual Coordinates<T, Dim> getDerivative(Coordinates<T, Dim> State) const { return Coordinates<T, Dim>(); }
	virtual Coordinates<T, Dim - 1> getConstants(Coordinates<T, Dim> StartCoords) const { return Coordinates<T, Dim - 1>(); }
	virtual Coordinates<T, Dim> getState(T Time, Coordinates<T, Dim - 1> Constants) const { return Coordinates<T, Dim>(); }
};

//------------------------------------------------HarmonicEquation----------------------------------------------------------------

/**
 * @brief class HarmonicEquation - equation of harmonic oscillations with frequency W.
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..   .
 *              HarmonicEquation - x + Ax + Bx + C = 0 
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = -Au - Bx - C
 *               
 * 		   getDerivative(x, u) == [u, -Au - Bx - C]
 *              
 */
template <typename T>
class HarmonicEquation : public DiffEquation<T, 3>
{
	T A_, B_, C_;

public:
	HarmonicEquation(T W) : DiffEquation<T, 3>(), A_(0), B_(W * W), C_(0) {};
	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -A_ * V - B_ * X - C_};
	}

	Coordinates<T, 2> getConstants(Coordinates<T, 3> StartCoords) const override
	{
		T W = sqrt(B_);
		T T0 = StartCoords[0], X0 = StartCoords[1], U0 = StartCoords[2];
	    T C1 = (X0 * W * sin(W * T0) + U0 * cos(W * T0)) / W;
	    T C2 = (X0 * W * cos(W * T0) - U0 * sin(W * T0)) / W;
		return Coordinates<T, 2>{C1, C2};
	}

	Coordinates<T, 3> getState(T Time, Coordinates<T, 2> Constants) const override
	{
		T W = sqrt(B_);
		T C1 = Constants[0], C2 = Constants[1];
		T X = C1 * sin(W * Time) + C2 * cos(W * Time);
		T V = C1 * W * cos(W * Time) - C2 * W * sin(W * Time);
		return Coordinates<T, 3>{Time, X, V}; 
	}
	// Frequency
	T W() const { return sqrt(B_); };
	// Attenuation
	T G() const { return C_; };
};

//------------------------------------------------PhysOscillEquation----------------------------------------------------------------

/**
 * @brief class PhysOscillEquation - equation of physical harmonic oscillations with frequency W.
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..   
 *              HarmonicEquation - x  + W^2 sin x + C = 0 
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = -W^2 sin x - C
 *               
 * 		   getDerivative(x, u) == [u, -W^2 sin x - C]
 *              
 */
template <typename T>
class PhysOscillEquation : public DiffEquation<T, 3>
{
	T W_, C_;

public:
	PhysOscillEquation(T W, T C = 0) : DiffEquation<T, 3>(), W_(W), C_(C) {};
	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -W_ * W_ * sin(X) - C_};
	}

	// Frequency
	T W() const { return W_; };
	// Attenuation
	T G() const { return C_; };
};

//------------------------------------------------HarmonicEquationWithFriction------------------------------------------------------

/**
 * @brief class HarmonicEquationWithFriction - HarmonicEquation with friction G.
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..    .
 *              HarmonicEquation - x + 2Gx + W^2 x = 0 
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = -2Gu - W^2 x
 *               
 * 		   getDerivative(x, u) == [u, -2Gu - W^2 x]
 *              
 */
template <typename T>
class HarmonicEquationWithFriction : public DiffEquation<T, 3>
{
	T W_, G_;

public:
	HarmonicEquationWithFriction(T W, T G) : DiffEquation<T, 3>(), W_(W), G_(G) {};
	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -2 * G_ * V - W_ * W_ * X};
	}

	// Coordinates<T, 2> getConstants(Coordinates<T, 3> StartCoords) const override
	// {
	// 	T T0 = StartCoords[0], X0 = StartCoords[1], U0 = StartCoords[2];
	//     T C1 = (X0 * W * sin(W * T0) + U0 * cos(W * T0)) / W;
	//     T C2 = (X0 * W * cos(W * T0) - U0 * sin(W * T0)) / W;
	// 	return Coordinates<T, 2>{C1, C2};
	// }

	// Coordinates<T, 3> getState(T Time, Coordinates<T, 2> Constants) const override
	// {
	// 	T C1 = Constants[0], C2 = Constants[1];
	// 	T X = C1 * sin(W_ * Time) + C2 * cos(W_ * Time);
	// 	T V = C1 * W_ * cos(W_ * Time) - C2 * W_ * sin(W_ * Time);
	// 	return Coordinates<T, 3>{Time, X, V}; 
	// }
	// Frequency
	T W() const { return W_; };
	// Attenuation
	T G() const { return G_; };
};


#endif // DIFF_EQUATION_H
