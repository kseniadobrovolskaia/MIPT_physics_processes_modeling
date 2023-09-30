#include "Solver.hpp"




//----------------------------------------------------Solver----------------------------------------------------------------------

void Solver::writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double Time = Start;
	while (Time < Stop)
	{
		auto Coords = this->getXV(Time);
		FileWithSolution.write((const char *)(&Time), sizeof(Time));
		FileWithSolution.write((const char *)(&Coords), sizeof(Coords));
		Time += DeltaT;
	}
}

//------------------------------------------------AnalyticalSolver----------------------------------------------------------------

void AnalyticalSolver::setStartConditions(double X0, double U0)
{
	double W = Equation_.W();
    A_ = U0 / W;
    B_ = X0;
}

Coordinates AnalyticalSolver::getXV(double Time) const
{
	double W = Equation_.W();
	return Coordinates{A_ * sin(W * Time) + B_ * cos(W + Time), A_ * W * cos(W * Time) - B_ * W * sin(W * Time)};
}

