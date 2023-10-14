#include "Solver.hpp"




//----------------------------------------------------Solver----------------------------------------------------------------------

void Solver::writeSolution(std::ofstream &FileWithSolution) const
{
	std::for_each(Trajectory_.begin(), Trajectory_.end(), [&](const Coordinates& K) 
			 	{ 
			 		FileWithSolution.write((const char *)(&K), sizeof(K)); 
			 	});
}

void Solver::writeEnergy(std::ofstream &FileWithSolution) const
{
	double W = Equation_.W();
	std::for_each(Trajectory_.begin(), Trajectory_.end(), [&](const Coordinates& K) 
			 	{
			 		double Energy = 0;
					Energy = K.U * K.U / 2 + W * W * K.X * K.X / 2;
					FileWithSolution.write((const char *)(&K.Time), sizeof(K.Time));
					FileWithSolution.write((const char *)(&Energy), sizeof(Energy));
			 	});
}

Coordinates Solver::getCoords(unsigned Step) const
{
	if (Step >= Trajectory_.size())
		throw std::logic_error("Сoordinates were not calculated at this time step");
	return Trajectory_[Step];
}

//------------------------------------------------AnalyticalSolver----------------------------------------------------------------

void AnalyticalSolver::setCoefficients(Coordinates StartCoords)
{
	double W = Equation_.W();
	double T0 = StartCoords.Time, X0 = StartCoords.X, U0 = StartCoords.U;
    A_ = (X0 * W * sin(W * T0) + U0 * cos(W * T0)) / W;
    B_ = (X0 * W * cos(W * T0) - U0 * sin(W * T0)) / W;
}

void AnalyticalSolver::calculateTrajectory(Coordinates StartCoords, TimeRange Range)
{
	double W = Equation_.W();
	double Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
	Trajectory_.clear();
	Trajectory_.reserve((Stop - Start) / DeltaT + 1);
	for (double T = Start; T < Stop; T += DeltaT)
		Trajectory_.emplace_back(T, A_ * sin(W * T) + B_ * cos(W * T), A_ * W * cos(W * T) - B_ * W * sin(W * T));
}

//------------------------------------------------EilerSolver----------------------------------------------------------------

void EilerSolver::calculateTrajectory(Coordinates StartCoords, TimeRange Range)
{
	double W = Equation_.W();
	double Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
	Trajectory_.clear();
	Trajectory_.reserve((Stop - Start) / DeltaT + 1);
	Coordinates K1, K0 = getStart(StartCoords, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.Time = T + DeltaT;
		K1.X    = K0.X + DeltaT * K0.U;
		K1.U    = K0.U - DeltaT * W * W * K0.X;
		Trajectory_.push_back(K0);
		K0   = K1;
	}
}

Coordinates EilerSolver::getStart(Coordinates StartCoords, double DeltaT) const
{
	double W = Equation_.W();
	Coordinates K0 = StartCoords, K1;
	for (double T = 0; T < StartCoords.Time; T += DeltaT)
	{
		K1.X = K0.X + DeltaT * K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;
		K0   = K1;
	}
	K0.Time = StartCoords.Time;
	return K0;
}

//------------------------------------------------HeunSolver----------------------------------------------------------------

void HeunSolver::calculateTrajectory(Coordinates StartCoords, TimeRange Range)
{
	double W = Equation_.W();
	double Start = Range.Start, Stop = Range.Stop, DeltaT = Range.DeltaT;
	Trajectory_.clear();
	Trajectory_.reserve((Stop - Start) / DeltaT + 1);
	Coordinates K1, K2, K0 = getStart(StartCoords, DeltaT);
	for (double T = Start; T < Stop; T += DeltaT)
	{
		K1.X    = K0.X + DeltaT *         K0.U;
		K1.U    = K0.U - DeltaT * W * W * K0.X;
		K2.Time = T + DeltaT;
		K2.X    = K0.X + DeltaT / 2 *         (K0.U + K1.U);
		K2.U    = K0.U - DeltaT / 2 * W * W * (K0.X + K1.X);
		Trajectory_.push_back(K0);
		K0   = K2;
	}
}

Coordinates HeunSolver::getStart(Coordinates StartCoords, double DeltaT) const
{
	double W = Equation_.W();
	Coordinates K0 = StartCoords, K1, K2;
	for (double T = 0; T < StartCoords.Time; T += DeltaT)
	{
		K1.X = K0.X + DeltaT *         K0.U;
		K1.U = K0.U - DeltaT * W * W * K0.X;
		K2.X = K0.X + DeltaT / 2 *         (K0.U + K1.U);
		K2.U = K0.U - DeltaT / 2 * W * W * (K0.X + K1.X);
		K0   = K2;
	}
	K0.Time = StartCoords.Time;
	return K0;
}

