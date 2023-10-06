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

void Solver::writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double Energy = 0;
	double W = Equation_.W();
	for (double T = Start; T < Stop; T += DeltaT)
	{
		auto Coords = this->getXV(T);
		Energy = Coords.U * Coords.U / 2 + W * W * Coords.X * Coords.X / 2;
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
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
	Coordinates K1, K0 = getStart(Start, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X = K0.X + DeltaT * K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&K0), sizeof(K0));
		K0   = K1;
	}
}

void EilerSolver::writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double Energy = 0;
	double W = Equation_.W();
	Coordinates K1, K0 = getStart(Start, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X = K0.X + DeltaT * K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;
		Energy = K0.U * K0.U / 2 + W * W * K0.X * K0.X / 2;
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
		K0   = K1;
	}
}

Coordinates EilerSolver::getStart(double Start, double DeltaT) const
{
	double W = Equation_.W();
	Coordinates K0 = StartXU, K1;	
	for (double T = 0; T < Start; T += DeltaT)
	{
		K1.X = K0.X + DeltaT * K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;
		K0   = K1;
	}
	return K0;
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

//------------------------------------------------HeunSolver----------------------------------------------------------------

void HeunSolver::writeSolution(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double W = Equation_.W();
	Coordinates K1, K2, K0 = getStart(Start, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X = K0.X + DeltaT *         K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;

		K2.X = K0.X + DeltaT / 2 *         (K0.U + K1.U);
		K2.U = K0.U - DeltaT / 2 * W * W * (K0.X + K1.X);

		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&K0), sizeof(K0));

		K0   = K2;
	}
}

void HeunSolver::writeEnergy(double Start, double Stop, double DeltaT, std::ofstream &FileWithSolution) const
{
	double Energy = 0;
	double W = Equation_.W();
	Coordinates K1, K2, K0 = getStart(Start, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X = K0.X + DeltaT *         K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;

		K2.X = K0.X + DeltaT / 2 *         (K0.U + K1.U);
		K2.U = K0.U - DeltaT / 2 * W * W * (K0.X + K1.X);

		Energy = K0.U * K0.U / 2 + W * W * K0.X * K0.X / 2;
		FileWithSolution.write((const char *)(&T), sizeof(T));
		FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
		K0   = K2;
	}
}

Coordinates HeunSolver::getStart(double Start, double DeltaT) const
{
	double W = Equation_.W();
	Coordinates K0 = StartXU, K1, K2;
	for (double T = 0; T < Start; T += DeltaT)
	{
		K1.X = K0.X + DeltaT *         K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;

		K2.X = K0.X + DeltaT / 2 *         (K0.U + K1.U);
		K2.U = K0.U - DeltaT / 2 * W * W * (K0.X + K1.X);
		K0   = K2;
	}
	return K0;
}

void HeunSolver::setStartConditions(double X0, double U0)
{
	StartXU.X = X0;
	StartXU.U = U0;
}

Coordinates HeunSolver::getXV(double Time) const
{
	return getStart(Time, DeltaT_);
}

