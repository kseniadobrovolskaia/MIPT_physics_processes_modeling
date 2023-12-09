#ifndef DIFF_EQUATION_H
#define DIFF_EQUATION_H



#include "DrivenForce.hpp"
#include "magic_enum.hpp"


enum class Models
{
	Math,
	Phys,
	MathWithFric,
	MathWithDriv
};


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
	virtual const std::basic_string_view<char> getName() const { return "BaseModel"; }
};

//------------------------------------------------HarmonicEquation----------------------------------------------------------------

/**
 * @brief class HarmonicEquation - equation of harmonic oscillations with frequency W.
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..    .                .
 *              HarmonicEquation - x  + Bx = 0
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = - Bx
 *               
 * 		   getDerivative(x, u) == [u, - Bx]
 *              
 */
template <typename T>
class HarmonicEquation : public DiffEquation<T, 3>
{
	T B_;

public:
	HarmonicEquation(T W) : DiffEquation<T, 3>(), B_(W * W) {};
	const std::basic_string_view<char> getName() const override { return magic_enum::enum_name(Models::Math); }

	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, - B_ * X};
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
};

//------------------------------------------------PhysOscillEquation----------------------------------------------------------------

/**
 * @brief class PhysOscillEquation - equation of physical harmonic oscillations with frequency W.
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..   
 *              HarmonicEquation - x  + W^2 sin x = 0 
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = -W^2 sin x
 *               
 * 		   getDerivative(x, u) == [u, -W^2 sin x]
 *              
 */
template <typename T>
class PhysOscillEquation : public DiffEquation<T, 3>
{
	T W_;

public:
	PhysOscillEquation(T W) : DiffEquation<T, 3>(), W_(W) {};
	const std::basic_string_view<char> getName() const override { return magic_enum::enum_name(Models::Phys); }

	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -W_ * W_ * sin(X)};
	}

	// Frequency
	T W() const { return W_; };
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
	const std::basic_string_view<char> getName() const override { return magic_enum::enum_name(Models::MathWithFric); }

	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -2 * G_ * V - W_ * W_ * X};
	}

	Coordinates<T, 2> getConstants(Coordinates<T, 3> StartCoords) const override
	{
		if (G_ > W_)
			return getConstantsFirst(StartCoords);
		else if (G_ < W_)
			return getConstantsSecond(StartCoords);
		return getConstantsThird(StartCoords);
	}

	Coordinates<T, 2> getConstantsFirst(Coordinates<T, 3> StartCoords) const
	{
		T X0 = StartCoords[1], U0 = StartCoords[2];
		T Alpha = sqrt(G_ * G_ - W_ * W_);

	    T C1 = (U0 + (Alpha + G_) * X0)  / (2 * Alpha);
	    T C2 = (-U0 + (Alpha - G_) * X0) / (2 * Alpha);
		return Coordinates<T, 2>{C1, C2};
	}

	Coordinates<T, 2> getConstantsSecond(Coordinates<T, 3> StartCoords) const
	{
		T X0 = StartCoords[1], U0 = StartCoords[2];
		T W0 = sqrt(-G_ * G_ + W_ * W_);

	    T C1 = X0;
	    T C2 = (U0 + G_* X0) / W0;
		return Coordinates<T, 2>{C1, C2};
	}

	Coordinates<T, 2> getConstantsThird(Coordinates<T, 3> StartCoords) const
	{
		T X0 = StartCoords[1], U0 = StartCoords[2];
		
	    T C1 = X0;
	    T C2 = U0 + G_ * X0;
		return Coordinates<T, 2>{C1, C2};
	}

	Coordinates<T, 3> getState(T Time, Coordinates<T, 2> Constants) const override
	{
		if (G_ > W_)
			return getStateFirst(Time, Constants);
		else if (G_ < W_)
			return getStateSecond(Time, Constants);
		return getStateThird(Time, Constants);
	}

	Coordinates<T, 3> getStateFirst(T Time, Coordinates<T, 2> Constants) const
	{
		T C1 = Constants[0], C2 = Constants[1];
		T Alpha = sqrt(G_ * G_ - W_ * W_);

		T X = exp(-G_ * Time) * (C1 * exp(Alpha * Time) + C2 * exp(-Alpha * Time));
		T V = -G_ * X + exp(-G_ * Time) * (Alpha * C1 * exp(Alpha * Time) - Alpha * C2 * exp(-Alpha * Time));
		return Coordinates<T, 3>{Time, X, V}; 
	}

	Coordinates<T, 3> getStateSecond(T Time, Coordinates<T, 2> Constants) const
	{
		T C1 = Constants[0], C2 = Constants[1];
		T W0 = sqrt(-G_ * G_ + W_ * W_);

		T X = exp(-G_ * Time) * (C1 * cos(W0 * Time) + C2 * sin(W0 * Time));
		T V = -G_ * X + W0 * exp(-G_ * Time) * (-C1 * sin(W0 * Time) + C2 * cos(W0 * Time));
		return Coordinates<T, 3>{Time, X, V}; 
	}

	Coordinates<T, 3> getStateThird(T Time, Coordinates<T, 2> Constants) const
	{
		T C1 = Constants[0], C2 = Constants[1];

		T X = exp(-G_ * Time) * (C1 + C2 * Time);
		T V = -G_ * X + exp(-G_ * Time) * C2;
		return Coordinates<T, 3>{Time, X, V}; 
	}

	// Frequency
	T W() const { return W_; };
	// Attenuation
	T G() const { return G_; };
};

//------------------------------------------------DrivenOscillatorEquation----------------------------------------------------------------

/**
 * @brief class DrivenOscillatorEquation - equation of driven oscillations with frequency W, friction G and driving force F(t, X, V).
 *                             
 *              Coordinates<T, 3> = {Time, X, V}
 *                                 ..     .                   .
 *              HarmonicEquation - x  + 2Gx + W^2 x = F(t, x, x)
 *                     
 *                              _  .
 *                             |   x = u
 *                            <    .
 *                             |_  u = -2Gu - W^2 x + F(t, x, u)
 *               
 * 		   getDerivative(x, u) == [u, -2Gu - W^2 x + F(t, x, u)]
 *              
 */
template <typename T>
class DrivenOscillatorEquation : public DiffEquation<T, 3>
{
	T W_, G_;
	const DrivenForce<T> &F_;


public:
	DrivenOscillatorEquation(T W, T G, const DrivenForce<T> &F) : DiffEquation<T, 3>(), W_(W), G_(G), F_(F) {};
	const std::basic_string_view<char> getName() const override { return magic_enum::enum_name(Models::MathWithDriv); }

	Coordinates<T, 3> getDerivative(Coordinates<T, 3> State) const override
	{
		T X = State[1];
		T V = State[2];
		return Coordinates<T, 3>{1, V, -2 * G_ * V - W_ * W_ * X + F_(State)};
	}

	Coordinates<T, 3> getState(T Time, Coordinates<T, 2> Constants) const override
	{
		T X = getX(Time, Constants);
		T V = getU(Time, Constants);
		return Coordinates<T, 3>{Time, X, V}; 
	}

    T getX(T Time, Coordinates<T, 2> Constants) const
	{
		T C1 = Constants[0], C2 = Constants[1];
		T W0 = F_.getW(), F = F_.getF(), W = W_, G = G_, t = Time;

   		T X = C1*exp(t*(-G + sqrt(-G + W)*sqrt(G + W))) 
   		+ C2*exp(-t*(G + sqrt(-G + W)*sqrt(G + W))) + 2*F*G*W0*sin(W0*t)/(4*pow(G, 2)*pow(W0, 2) 
   		+ pow(W, 4) - 2*pow(W, 2)*pow(W0, 2) + pow(W0, 4))
   		+ F*pow(W, 2)*cos(W0*t)/(4*pow(G, 2)*pow(W0, 2) 
   		+ pow(W, 4) - 2*pow(W, 2)*pow(W0, 2) + pow(W0, 4)) 
   		- F*pow(W0, 2)*cos(W0*t)/(4*pow(G, 2)*pow(W0, 2) + pow(W, 4) 
   		- 2*pow(W, 2)*pow(W0, 2) + pow(W0, 4));

   		return X;
	}

    T getU(T Time, Coordinates<T, 2> Constants) const
	{
		T C1 = Constants[0], C2 = Constants[1];
		T W0 = F_.getW(), F = F_.getF(), W = W_, G = G_, t = Time;
		
		T U = C1*(-G + sqrt(-G + W)*sqrt(G + W))*exp(t*(-G + sqrt(-G + W)*sqrt(G + W))) 
		+ C2*(-G - sqrt(-G + W)*sqrt(G + W))*exp(-t*(G + sqrt(-G + W)*sqrt(G + W))) 
		+ 2*F*G*pow(W0, 2)*cos(W0*t)/(4*pow(G, 2)*pow(W0, 2) + pow(W, 4) - 2*pow(W, 2)*pow(W0, 2) 
		+ pow(W0, 4)) - F*pow(W, 2)*W0*sin(W0*t)/(4*pow(G, 2)*pow(W0, 2) + pow(W, 4) 
		- 2*pow(W, 2)*pow(W0, 2) + pow(W0, 4)) + F*pow(W0, 3)*sin(W0*t)/(4*pow(G, 2)*pow(W0, 2) 
		+ pow(W, 4) - 2*pow(W, 2)*pow(W0, 2) + pow(W0, 4));

		return U;
	}

	Coordinates<T, 2> getConstants(Coordinates<T, 3> StartCoords) const override
	{
		T X0 = StartCoords[1], U0 = StartCoords[2];
		T W0 = F_.getW(), F = F_.getF(), W = W_, G = G_;

		T P1 = (-G + sqrt(-G + W)*sqrt(G + W)), P2 = (-G - sqrt(-G + W)*sqrt(G + W));
		T P3 = F*W*W/(4*G*G*W0*W0 + pow(W, 4) - 2*W*W*W0*W0 + pow(W0, 4)) 
		- F*W0*W0/(4*G*G*W0*W0 + pow(W, 4) - 2*W*W*W0*W0 + pow(W0, 4));
		T P4 = 2*F*G*W0*W0/(4*G*G*W0*W0 + pow(W, 4) - 2*W*W*W0*W0 + pow(W0, 4));

		T C2 = (U0 - P1*X0 + P3*P1 - P4) / (P2 - P1);
	    T C1 = X0 - C2 - P3;
	    
		return Coordinates<T, 2>{C1, C2};
	}

	// Frequency
	T W() const { return W_; };
	// Attenuation
	T G() const { return G_; };
	// Driving force
	DrivenForce<T> F() const { return F_; };
};


#endif // DIFF_EQUATION_H
