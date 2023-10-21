#ifndef SOLVER_H
#define SOLVER_H


#include "DiffEquation.hpp"




template <typename T>
struct TimeRange
{
	T Start;
	T Stop;
	T DeltaT;

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
	const HarmonicEquation<T> &Equation_;
	SequenceOfStates<T, Dim> Trajectory_;

public:

	Solver(const HarmonicEquation<T> &Equation) : Equation_(Equation) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range) { return; };
	void writeSolution(std::ofstream &FileWithSolution) const
	{
		std::for_each(Solver<T, Dim>::Trajectory_.begin(), Solver<T, Dim>::Trajectory_.end(), [&](const Coordinates<T, Dim>& K) 
				 	{ 
				 		FileWithSolution.write((const char *)(&K), sizeof(K)); 
				 	});
	}

	void writeEnergy(std::ofstream &FileWithSolution) const
	{
		T W = Solver<T, Dim>::Equation_.W();
		std::for_each(Solver<T, Dim>::Trajectory_.begin(), Solver<T, Dim>::Trajectory_.end(), [&](const Coordinates<T, Dim>& K) 
				 	{
				 		T Energy = 0, X = K[1], V = K[2];
						Energy = V * V / 2 + W * W * X * X / 2;
						FileWithSolution.write((const char *)(&K[0]), sizeof(K[0]));
						FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
				 	});
	}

	Coordinates<T, Dim> getCoords(unsigned Step) const
	{
		if (Step >= Solver<T, Dim>::Trajectory_.size())
			throw std::logic_error("Сoordinates were not calculated at this time step");
		return Solver<T, Dim>::Trajectory_[Step];
	}

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
	AnalyticalSolver(HarmonicEquation<T> &Equation) : Solver<T, Dim>(Equation) {};
	void setConstants(Coordinates<T, Dim> StartCoords)
	{
		Constants_ = Solver<T, Dim>::Equation_.getConstants(StartCoords);
	}

	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range)
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
	EilerSolver(HarmonicEquation<T> &Equation, T DeltaT = 0.01) : Solver<T, Dim>(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range)
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		Coordinates<T, Dim> K1, K0 = getStart(StartCoords, DeltaT);
		for (T Time = Start; Time < Stop; Time += DeltaT)
		{
			K1    = K0   + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K1[0] = Time + DeltaT;
			Solver<T, Dim>::Trajectory_.push_back(K0);
			K0    = K1;
		}
	}

	Coordinates<T, Dim> getStart(Coordinates<T, Dim> StartCoords, T DeltaT) const
	{
		Coordinates<T, Dim> K0 = StartCoords, K1;
		for (T Time = 0; Time < StartCoords[0]; Time += DeltaT)
		{
			K1 = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K0 = K1;
		}
		K0[0] = StartCoords[0];
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
	HeunSolver(HarmonicEquation<T> &Equation, T DeltaT = 0.01) : Solver<T, Dim>(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates<T, Dim> StartCoords, TimeRange<T> Range)
	{
		T Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
		Solver<T, Dim>::Trajectory_.clear();
		Solver<T, Dim>::Trajectory_.reserve((Stop - Start) / DeltaT + 1);
		Coordinates<T, Dim> K1, K2, K0 = getStart(StartCoords, DeltaT);
		for (T Time = Start; Time < Stop; Time += DeltaT)
		{
			K1    = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K2[0] = Time + DeltaT;
			K2    = K0 + DeltaT / 2 * (Solver<T, Dim>::Equation_.getDerivative(K0) + Solver<T, Dim>::Equation_.getDerivative(K1));
			Solver<T, Dim>::Trajectory_.push_back(K0);
			K0   = K2;
		}
	}

	Coordinates<T, Dim> getStart(Coordinates<T, Dim> StartCoords, T DeltaT) const
	{
		Coordinates<T, Dim> K0 = StartCoords, K1, K2;
		for (T Time = 0; Time < StartCoords[0]; Time += DeltaT)
		{
			K1 = K0 + DeltaT * Solver<T, Dim>::Equation_.getDerivative(K0);
			K2 = K0 + DeltaT / 2 * (Solver<T, Dim>::Equation_.getDerivative(K0) + Solver<T, Dim>::Equation_.getDerivative(K1));
			K0 = K2;
		}
		K0[0] = StartCoords[0];
		return K0;
	}
};


#endif // SOLVER_H
