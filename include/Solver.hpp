#ifndef SOLVER_H
#define SOLVER_H


#include "DiffEquation.hpp"




template <typename T>
struct TimeRange
{
	T Start;
	T Stop;
	T DeltaT;

	TimeRange() {};

	TimeRange(T S, T St, T DT = 0.01) : Start(S), Stop(St), DeltaT(DT)
	{
		if (Start < 0)
			throw std::logic_error("Start time can't be negative");
		if (Stop < 0)
			throw std::logic_error("Stop time can't be negative");
		if (DeltaT < 0)
			throw std::logic_error("DeltaT can't be negative");
		if (Start > Stop)
			throw std::logic_error("Start time can't be bigger than Stop time");
		if (DeltaT > (Stop - Start))
			throw std::logic_error("DeltaT can't be bigger than (Stop - Start)");
	}
};

template <typename T, unsigned Dim>
using SequenceOfStates = std::vector<Coordinates<T, Dim>>;

template <typename T, unsigned Dim>
using SequenceOfConstants = Coordinates<T, Dim - 1>;

//----------------------------------------------------Solver----------------------------------------------------------------------

/**
 * @brief class Solver - abstract class for solving a Dim-order differential equation Equation_,
 *                       with a concrete method that each child class defines
 */
template <typename T, unsigned Dim>
class Solver
{
protected:
	const DiffEquation<T, Dim> &Equation_;
	SequenceOfStates<T, Dim> Trajectory_;

public:

	Solver(const DiffEquation<T, Dim> &Equation) : Equation_(Equation) {};
	virtual void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) { return; };
	virtual bool isCalculated() const { return !Trajectory_.empty(); }
	virtual void writeSolution(std::ofstream &FileWithSolution) const
	{
		std::for_each(Solver<T, Dim>::Trajectory_.begin(), Solver<T, Dim>::Trajectory_.end(), [&](const Coordinates<T, Dim>& K) 
				 	{ 
				 		FileWithSolution.write((const char *)(&K), sizeof(K)); 
				 	});
	}

	virtual void writeEnergy(std::ofstream &FileWithSolution) const
	{
		T W = static_cast<const HarmonicEquation<T>&>(Solver<T, Dim>::Equation_).W();
		std::for_each(Solver<T, Dim>::Trajectory_.begin(), Solver<T, Dim>::Trajectory_.end(), [&](const Coordinates<T, Dim>& K) 
				 	{
				 		T Energy = 0, X = K[1], V = K[2];
						Energy = V * V / 2 + W * W * X * X / 2;
						FileWithSolution.write((const char *)(&K[0]), sizeof(K[0]));
						FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
				 	});
	}

	virtual Coordinates<T, Dim> getCoords(unsigned Step) const
	{
		if (Step >= Solver<T, Dim>::Trajectory_.size())
			throw std::logic_error("Сoordinates were not calculated at this time step");
		return Solver<T, Dim>::Trajectory_[Step];
	}

	virtual const SequenceOfStates<T, Dim> &getTrajectory() const { return Trajectory_; }
	virtual const DiffEquation<T, Dim>     &getEquation()   const { return Equation_;   }

};

//------------------------------------------------AnalyticalSolver----------------------------------------------------------------

/**
 * @brief class AnalyticalSolver - solves the Equation_ analytically
 */
template <typename T, unsigned Dim>
class AnalyticalSolver : public Solver<T, Dim>
{
	SequenceOfConstants<T, Dim> Constants_;

public:
	AnalyticalSolver(DiffEquation<T, Dim> &Equation) : Solver<T, Dim>(Equation) {};
	void setConstants(Coordinates<T, Dim> StartCoords)
	{
		Constants_ = Solver<T, Dim>::Equation_.getConstants(StartCoords);
	}

	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) override
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		for (T Time = Start; Time < Stop; Time += DeltaT)
			Solver<T, Dim>::Trajectory_.emplace_back(Solver<T, Dim>::Equation_.getState(Time, Constants_));
	}
};
	
//---------------------------------------------------EilerSolver----------------------------------------------------------------

/**
 * @brief class EilerSolver - solves the Equation_ using Euler's numerical iterative method
 */
template <typename T, unsigned Dim>
class EilerSolver : public Solver<T, Dim>
{
	T DeltaT_;

public:
	EilerSolver(DiffEquation<T, Dim> &Equation, T DeltaT = 0.01) : Solver<T, Dim>(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) override
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		Coordinates<T, Dim> K1, K0 = getStart(StartCoords, DeltaT);
		for (T Time = Start; Time < Stop; Time += DeltaT)
		{
			K1    = K0   + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			Solver<T, Dim>::Trajectory_.push_back(K0);
			K0    = K1;
		}
	}

	Coordinates<T, Dim> getStart(Coordinates<T, Dim> StartCoords, T DeltaT) const
	{
		Coordinates<T, Dim> K1, K0 = StartCoords;
		for (T Time = 0; Time < StartCoords[0]; Time += DeltaT)
		{
			K1 = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K0 = K1;
		}
		return K0;
	}
};

//---------------------------------------------------HeunSolver-------------------------------------------------------------------

/**
 * @brief class HeunSolver - solves the Equation_ using specified Euler's method
 *                            called Heun's scheme. Iterative numerical solution in two stages (predictive-corrector)
 *                            based on the trapezoid method.
 */
template <typename T, unsigned Dim>
class HeunSolver : public Solver<T, Dim>
{
	T DeltaT_;

public:
	HeunSolver(DiffEquation<T, Dim> &Equation, T DeltaT = 0.01) : Solver<T, Dim>(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) override
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		Coordinates<T, Dim> K1, K2, K0 = getStart(StartCoords, DeltaT);
		for (T Time = Start; Time < Stop; Time += DeltaT)
		{
			K1    = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K2    = K0 + DeltaT / 2 * (Solver<T, Dim>::Equation_.getDerivative(K0) + Solver<T, Dim>::Equation_.getDerivative(K1));
			Solver<T, Dim>::Trajectory_.push_back(K0);
			K0   = K2;
		}
	}

	Coordinates<T, Dim> getStart(Coordinates<T, Dim> StartCoords, T DeltaT) const
	{
		Coordinates<T, Dim> K1, K2, K0 = StartCoords;
		for (T Time = 0; Time < StartCoords[0]; Time += DeltaT)
		{
			K1 = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K2 = K0 + DeltaT / 2 * (Solver<T, Dim>::Equation_.getDerivative(K0) + Solver<T, Dim>::Equation_.getDerivative(K1));
			K0 = K2;
		}
		return K0;
	}
};

//---------------------------------------------------RungeKuttaSolver-------------------------------------------------------------------

/**
 * @brief class RungeKuttaSolver - solves the Equation_ using Runge-Kutta method,
 *                                 that helps to numerically integrate the equations 
 *                                 by calculating a new value in four steps.
 */
template <typename T, unsigned Dim>
class RungeKuttaSolver : public Solver<T, Dim>
{
	T DeltaT_;

public:
	RungeKuttaSolver(DiffEquation<T, Dim> &Equation, T DeltaT = 0.01) : Solver<T, Dim>(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) override
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		Coordinates<T, Dim> K1, K2, D1, D2, D3, D4, K0 = getStart(StartCoords, DeltaT);
		for (T Time = Start; Time < Stop; Time += DeltaT)
		{
			D1 = Solver<T, Dim>::Equation_.getDerivative(K0);
			D2 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT / 2 * D1);
			D3 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT / 2 * D2);
			D4 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT * D3);
			K1 = K0 + DeltaT / 6 * (D1 + 2 * D2 + 2 * D3 + D4);
			Solver<T, Dim>::Trajectory_.push_back(K0);
			K0 = K1;
		}
	}

	Coordinates<T, Dim> getStart(Coordinates<T, Dim> StartCoords, T DeltaT) const
	{
		Coordinates<T, Dim> K1, K2, D1, D2, D3, D4, K0 = StartCoords;
		for (T Time = 0; Time < StartCoords[0]; Time += DeltaT)
		{
			D1 = Solver<T, Dim>::Equation_.getDerivative(K0);
			D2 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT / 2 * D1);
			D3 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT / 2 * D2);
			D4 = Solver<T, Dim>::Equation_.getDerivative(K0 + DeltaT * D3);
			K1 = K0 + DeltaT / 6 * (D1 + 2 * D2 + 2 * D3 + D4);
			K0 = K1;
		}
		return K0;
	}
};


#endif // SOLVER_H
