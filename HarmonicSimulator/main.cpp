#include "SolverWithName.hpp"
#include "DrivenForce.hpp"
#include "json.hpp"



using TypeForCoords = double;
#define Dim 3


std::string getConfigName(const int argc, const char *argv[]);
void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &G, double &F, double &W0, Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range, std::string &Model, std::string &Solver);
void writeSolutionAndEnergyForMethod(Solvers Solver, DiffEquation<TypeForCoords, Dim> &Equation, 
	                                     Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range);




int main(const int argc, const char *argv[])
{
	double W, G, F, W0;
	TimeRange<TypeForCoords> Range;
	Coordinates<TypeForCoords, Dim> StartCoords;
	std::string ModelStr, SolverStr;

	const std::string ConfigFileName = getConfigName(argc, argv);
	getStartConditionsFromConfigFile(ConfigFileName, W, G, F, W0, StartCoords, Range, ModelStr, SolverStr);

	//---------------Create_Equations------------------------------------

	auto Model = magic_enum::enum_cast<Models>(ModelStr);
	if (!Model.has_value())
	{
		std::cout << "We dont know this Model: " << ModelStr << "\n";
		return 0;
	}
	auto Solver = magic_enum::enum_cast<Solvers>(SolverStr);
	if (!Solver.has_value())
	{
		std::cout << "We dont know this Solver: " << SolverStr << "\n";
		return 0;
	}

	switch (Model.value())
	{
		case Models::Math:
		{
			HarmonicEquation<TypeForCoords> MathOscilliator(W);
			writeSolutionAndEnergyForMethod(Solver.value(), MathOscilliator, StartCoords, Range);
		}
		case Models::Phys:
		{
			PhysOscillEquation<TypeForCoords> PhysOscilliator(W);
			writeSolutionAndEnergyForMethod(Solver.value(), PhysOscilliator, StartCoords, Range);
		}
		case Models::MathWithFric:
		{
			HarmonicEquationWithFriction<TypeForCoords> MathWithFriction(W, G);
			writeSolutionAndEnergyForMethod(Solver.value(), MathWithFriction, StartCoords, Range);
		}
		case Models::MathWithDriv:
		{
			auto DrivenForceLambda = [](TypeForCoords F, TypeForCoords W0, Coordinates<TypeForCoords, 3> State) -> TypeForCoords { return F * cos(W0 * State[0]); };
			auto Force = DrivenForce<TypeForCoords>(F, W0, DrivenForceLambda);
			DrivenOscillatorEquation<TypeForCoords> MathWithDriven(W, G, Force);
			writeSolutionAndEnergyForMethod(Solver.value(), MathWithDriven, StartCoords, Range);
		}
	}

	return 0;
}




std::string getConfigName(const int argc, const char *argv[])
{
	if (argc < 1)
		throw std::logic_error("There are not Config File Name in args");
	return argv[1];
}

void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &G, double &F, double &W0, Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range, std::string &Model, std::string &Solver)
{
	std::ifstream ConfigFile(ConfigFileName);
	nlohmann::json Config = nlohmann::json::parse(ConfigFile);

	Model = Config["Model"];
	Solver = Config["Solver"];
	W = Config["W"];
	G = Config["G"];
	F = Config["F"];
	W0 = Config["W0"];
	StartCoords[0] = Config["T0"];
	StartCoords[1] = Config["X0"];
	StartCoords[2] = Config["V0"];
	Range.Start = Config["Start"];
	Range.Stop = Config["Stop"];
	Range.DeltaT = Config["Step"];
}

void writeSolutionAndEnergyForMethod(const Solvers Solver, DiffEquation<TypeForCoords, Dim> &Equation, 
	                                     Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range)
{
	const std::string EquationName(Equation.getName());
	const std::string SolverName(magic_enum::enum_name(Solver));

	switch (Solver)
	{
		case Solvers::Analitic:
		{
			AnalyticalSolver<TypeForCoords, Dim> Analitic(Equation);
			SolverWithName<TypeForCoords, Dim> AnaliticWithMath(SolverName, EquationName, Analitic);
			AnaliticWithMath.writeSolutionAndEnergy(StartCoords, Range);
		}
		case Solvers::Eiler:
		{
			EilerSolver<TypeForCoords, Dim> Eiler(Equation, Range.DeltaT);
			SolverWithName<TypeForCoords, Dim> EilerWithMath(SolverName, EquationName, Eiler);
			EilerWithMath.writeSolutionAndEnergy(StartCoords, Range);
		}
		case Solvers::Heun:
		{
			HeunSolver<TypeForCoords, Dim> Heun(Equation, Range.DeltaT);
			SolverWithName<TypeForCoords, Dim> HeunWithMath(SolverName, EquationName, Heun);
			HeunWithMath.writeSolutionAndEnergy(StartCoords, Range);
		}
		case Solvers::RungeKutta:
		{
			RungeKuttaSolver<TypeForCoords, Dim> RungeKutta(Equation, Range.DeltaT);
			SolverWithName<TypeForCoords, Dim> RungeKuttaWithMath(SolverName, EquationName, RungeKutta);
			RungeKuttaWithMath.writeSolutionAndEnergy(StartCoords, Range);
		}
	}
}
