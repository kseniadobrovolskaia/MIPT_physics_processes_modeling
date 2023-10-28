#ifndef SOLVER_WITH_NAME_H
#define SOLVER_WITH_NAME_H


#include "Solver.hpp"




template <typename T, unsigned Dim>
struct SolverWithName
{
	const std::string SolverName_;
	const std::string EquationName_;
	Solver<T, Dim> &Solver_;

	SolverWithName(const std::string SolverName, const std::string EquationName, Solver<T, Dim> &Solver) :
	SolverName_(SolverName), EquationName_(EquationName), Solver_(Solver) {};

	void writeSolution(Coordinates<T, Dim> StartCoords, TimeRange<T> Range)
	{
		if (SolverName_ == "Analitic")
			static_cast<AnalyticalSolver<T, Dim>&>(Solver_).setConstants(StartCoords);
		Solver_.calculateTrajectory(StartCoords, Range);

		std::string FileSolutionName = SolverName_ + EquationName_ + ".bin";
		std::ofstream FileSolution(FileSolutionName, std::ios::binary);
		Solver_.writeSolution(FileSolution);
	}

	void writeEnergy() const
	{
		if (!Solver_.isCalculated())
			throw std::logic_error("You must first writeSolution before calculating the energy");
		std::string FileEnergyName = SolverName_ + EquationName_ + "Energy.bin";
		std::ofstream FileEnergy(FileEnergyName, std::ios::binary);
		Solver_.writeEnergy(FileEnergy);
	}

	void writeSolutionAndEnergy(Coordinates<T, Dim> &StartCoords, TimeRange<T> &Range)
	{
		if ((EquationName_ == "Phys") && (SolverName_ == "Analitic"))
			return;
		writeSolution(StartCoords, Range);
		writeEnergy();
	}
};


#endif // SOLVER_WITH_NAME_H
