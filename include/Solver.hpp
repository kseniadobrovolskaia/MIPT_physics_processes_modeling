#ifndef SOLVER_H
#define SOLVER_H


#include "DiffEquation.hpp"



struct Coordinates
{
	double X;
	double U;
};

//----------------------------------------------------Solver----------------------------------------------------------------------

/**
 * @brief class Solver - abstract class for solving a second-order differential equation Equation_,
 *                       with a concrete method that each child class defines
 */
class Solver
{
protected:
	const DiffEquation &Equation_;

public:

	Solver(const DiffEquation &Equation) : Equation_(Equation) {};
	virtual void setStartConditions(double X0, double U0) = 0;
	virtual void writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const;
	virtual Coordinates getXV(double Time) const = 0;
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
	void setStartConditions(double X0, double U0) override;
	Coordinates getXV(double Time) const override;
};

#endif // SOLVER_H
