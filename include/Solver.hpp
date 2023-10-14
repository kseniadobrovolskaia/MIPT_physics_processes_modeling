#ifndef SOLVER_H
#define SOLVER_H


#include "DiffEquation.hpp"



struct Coordinates
{
	double Time;
	double X;
	double U;
};

struct TimeRange
{
	double Start;
	double Stop;
	double DeltaT;

	TimeRange(double S, double St, double DT = 0.01) : Start(S), Stop(St), DeltaT(DT)
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

using SequenceOfStates = std::vector<Coordinates>;

//----------------------------------------------------Solver----------------------------------------------------------------------

/**
 * @brief class Solver - abstract class for solving a second-order differential equation Equation_,
 *                       with a concrete method that each child class defines
 */
class Solver
{
protected:
	const DiffEquation &Equation_;
	SequenceOfStates Trajectory_;

public:

	Solver(const DiffEquation &Equation) : Equation_(Equation) {};
	virtual void calculateTrajectory(Coordinates StartCoords, TimeRange Range) = 0;
	virtual void writeSolution(std::ofstream &FileWithSolution) const;
	virtual void writeEnergy(std::ofstream &FileWithSolution) const;
	virtual Coordinates getCoords(unsigned Step) const;
};

//------------------------------------------------AnalyticalSolver----------------------------------------------------------------

/**
 * @brief class AnalyticalSolver - solves the Equation_ analytically
 */
class AnalyticalSolver : public Solver
{
	double A_, B_;

public:
	AnalyticalSolver(DiffEquation &Equation) : Solver(Equation) {};
	void setCoefficients(Coordinates StartCoords);
	void calculateTrajectory(Coordinates StartCoords, TimeRange Range) override;
};
	
//---------------------------------------------------EilerSolver----------------------------------------------------------------

/**
 * @brief class EilerSolver - solves the Equation_ using Euler's numerical iterative method
 */
class EilerSolver : public Solver
{
	double DeltaT_;

public:
	EilerSolver(DiffEquation &Equation, double DeltaT = 0.01) : Solver(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates StartCoords, TimeRange Range) override;
	Coordinates getStart(Coordinates StartCoords, double DeltaT) const;
};

//---------------------------------------------------HeunSolver----------------------------------------------------------------

/**
 * @brief class HeunSolver - solves the Equation_ using specified Euler's method
 *                            called Heun's scheme. Iterative numerical solution in two stages (predictive-corrector)
 *                            based on the trapezoid method.


 */
class HeunSolver : public Solver
{
	double DeltaT_;

public:
	HeunSolver(DiffEquation &Equation, double DeltaT = 0.01) : Solver(Equation), DeltaT_(DeltaT) {};
	void calculateTrajectory(Coordinates StartCoords, TimeRange Range) override;
	Coordinates getStart(Coordinates StartCoords, double DeltaT) const;
};


#endif // SOLVER_H
