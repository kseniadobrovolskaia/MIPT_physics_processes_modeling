#include "SolverWithName.hpp"
#include "json.hpp"




using TypeForCoords = double;
#define Dim 3


std::string getConfigName(const int argc, const char *argv[]);
void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range);
void writeSolutionAndEnergyForAllMethods(const std::string EquationName, DiffEquation<TypeForCoords, Dim> &Equation, 
	                                     Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range);


const std::string AnaliticName = "Analitic";
const std::string EilerName = "Eiler";
const std::string HeunName = "Heun";
const std::string RungeKuttaName = "RungeKutta";

const std::string MathOscilliatorName = "Math";
const std::string PhysOscilliatorName = "Phys";
const std::string MathWithFrictionOscilliatorName = "MathWithFric";




int main(const int argc, const char *argv[])
{
	double W;
	TimeRange<TypeForCoords> Range;
	Coordinates<TypeForCoords, Dim> StartCoords;

	const std::string ConfigFileName = getConfigName(argc, argv);
	getStartConditionsFromConfigFile(ConfigFileName, W, StartCoords, Range);

	//---------------Create_Equations------------------------------------

	HarmonicEquation<TypeForCoords> MathOscilliator(W);
	writeSolutionAndEnergyForAllMethods(MathOscilliatorName, MathOscilliator, StartCoords, Range);

	PhysOscillEquation<TypeForCoords> PhysOscilliator(W);
	writeSolutionAndEnergyForAllMethods(PhysOscilliatorName, PhysOscilliator, StartCoords, Range);

	HarmonicEquationWithFriction<TypeForCoords> MathWithFriction(W, 0.2);
	writeSolutionAndEnergyForAllMethods(MathWithFrictionOscilliatorName, MathWithFriction, StartCoords, Range);
	return 0;
}




std::string getConfigName(const int argc, const char *argv[])
{
	if (argc < 1)
		throw std::logic_error("There are not Config File Name in args");
	return argv[1];
}

void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range)
{
	std::ifstream ConfigFile(ConfigFileName);
	nlohmann::json Config = nlohmann::json::parse(ConfigFile);

	W = Config["W"];
	StartCoords[0] = Config["T0"];
	StartCoords[1] = Config["X0"];
	StartCoords[2] = Config["V0"];
	Range.Start = Config["Start"];
	Range.Stop = Config["Stop"];
	Range.DeltaT = Config["Step"];
}

void writeSolutionAndEnergyForAllMethods(const std::string EquationName, DiffEquation<TypeForCoords, Dim> &Equation, 
	                                     Coordinates<TypeForCoords, Dim> &StartCoords, TimeRange<TypeForCoords> &Range)
{
	AnalyticalSolver<TypeForCoords, Dim> Analitic(Equation);
	SolverWithName<TypeForCoords, Dim> AnaliticWithMath(AnaliticName, EquationName, Analitic);
	AnaliticWithMath.writeSolutionAndEnergy(StartCoords, Range);

	EilerSolver<TypeForCoords, Dim> Eiler(Equation, Range.DeltaT);
	SolverWithName<TypeForCoords, Dim> EilerWithMath(EilerName, EquationName, Eiler);
	EilerWithMath.writeSolutionAndEnergy(StartCoords, Range);

	HeunSolver<TypeForCoords, Dim> Heun(Equation, Range.DeltaT);
	SolverWithName<TypeForCoords, Dim> HeunWithMath(HeunName, EquationName, Heun);
	HeunWithMath.writeSolutionAndEnergy(StartCoords, Range);

	RungeKuttaSolver<TypeForCoords, Dim> RungeKutta(Equation, Range.DeltaT);
	SolverWithName<TypeForCoords, Dim> RungeKuttaWithMath(RungeKuttaName, EquationName, RungeKutta);
	RungeKuttaWithMath.writeSolutionAndEnergy(StartCoords, Range);
}
