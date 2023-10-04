#include "Solver.hpp"




//----------------------------------------------------Solver----------------------------------------------------------------------

void Solver::writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	for (double T = Start; T < Stop; T += DeltaT)
	{
		auto Coords = this->getXV(T);
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&Coords), sizeof(Coords));
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
	return Coordinates{A_ * sin(W * Time) + B_ * cos(W * Time), A_ * W * cos(W * Time) - B_ * W * sin(W * Time)};
}

//------------------------------------------------EilerSolver----------------------------------------------------------------

void EilerSolver::writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double W = Equation_.W();
	Coordinates K1, K = getStart(Start, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X = K.X + DeltaT * K.U;
		K1.U = K.U - DeltaT * W * W * K.X;
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&K1), sizeof(K1));
		K = K1;
	}
}

Coordinates EilerSolver::getStart(double Start, double DeltaT) const
{
	Coordinates K = StartXU, K1;
	double W = Equation_.W();	
	for (double T = 0; T < Start; T += DeltaT)
	{
		K1.X = K.X + DeltaT * K.U;
		K1.U = K.U - DeltaT * W * W * K.X;
		K = K1;
	}
	return K1;
}

void EilerSolver::setStartConditions(double X0, double U0)
{
	StartXU.X = X0;
	StartXU.U = U0;
}

Coordinates EilerSolver::getXV(double Time) const
{
	return getStart(Time, DeltaT_);
}

