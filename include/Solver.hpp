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
	virtual void writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const;
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
	
//---------------------------------------------------EilerSolver----------------------------------------------------------------

/**
 * @brief class EilerSolver - solves the Equation_ using Euler's numerical iterative method
 */
class EilerSolver : public Solver
{
	Coordinates StartXU;
	double DeltaT_;
	
public:
	EilerSolver(DiffEquation &Equation, double DeltaT = 0.01) : Solver(Equation), DeltaT_(DeltaT) {};
	void setStartConditions(double X0, double U0) override;
	void writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const override;
	void writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const override;
	Coordinates getStart(double Start, double DeltaT) const;
	Coordinates getXV(double Time) const override;
};

//---------------------------------------------------HeunSolver----------------------------------------------------------------

/**
 * @brief class HeunSolver - solves the Equation_ using specified Euler's method
 *                            called Heun's scheme. Iterative numerical solution in two stages (predictive-corrector)
 *                            based on the trapezoid method.


 */
class HeunSolver : public Solver
{
	Coordinates StartXU;
	double DeltaT_;
	
public:
	HeunSolver(DiffEquation &Equation, double DeltaT = 0.01) : Solver(Equation), DeltaT_(DeltaT) {};
	void setStartConditions(double X0, double U0) override;
	void writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const override;
	void writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const override;
	Coordinates getStart(double Start, double DeltaT) const;
	Coordinates getXV(double Time) const override;
};


#endif // SOLVER_H
