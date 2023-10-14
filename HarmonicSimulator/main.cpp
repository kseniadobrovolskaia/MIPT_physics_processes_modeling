#include "Solver.hpp"
#include "json.hpp"



std::string getConfigName(const int argc, const char *argv[]);
void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &T0, double &X0, double &U0, double &Start, double &Stop, double &DeltaT);


int main(const int argc, const char *argv[])
{
	double W, T0, X0, U0, Start, Stop, DeltaT;
	const char *FileNameAnalitic = "Analitic.bin";
	const char *FileNameEiler = "Eiler.bin";
	const char *FileNameHeun = "Heun.bin";

	const char *FileNameEnergyAnalitic = "AnaliticEnergy.bin";
 	const char *FileNameEnergyEiler = "EilerEnergy.bin";
 	const char *FileNameEnergyHeun = "HeunEnergy.bin";

	std::ofstream FileAnalitic(FileNameAnalitic, std::ios::binary);
	std::ofstream FileEiler(FileNameEiler, std::ios::binary);
	std::ofstream FileHeun(FileNameHeun, std::ios::binary);

	std::ofstream FileEnergyAnalitic(FileNameEnergyAnalitic, std::ios::binary);
	std::ofstream FileEnergyEiler(FileNameEnergyEiler, std::ios::binary);
	std::ofstream FileEnergyHeun(FileNameEnergyHeun, std::ios::binary);

	const std::string ConfigFileName = getConfigName(argc, argv);
	getStartConditionsFromConfigFile(ConfigFileName, W, T0, X0, U0, Start, Stop, DeltaT);

	TimeRange Range(Start, Stop, DeltaT);
	Coordinates StartCoords(T0, X0, U0);

	HarmonicEquation Oscilliator(W);

	AnalyticalSolver Analitic(Oscilliator);
	Analitic.setCoefficients(StartCoords);
	Analitic.calculateTrajectory(StartCoords, Range);

	EilerSolver Eiler(Oscilliator, DeltaT);
	Eiler.calculateTrajectory(StartCoords, Range);

	HeunSolver Heun(Oscilliator, DeltaT);
	Heun.calculateTrajectory(StartCoords, Range);

	Analitic.writeSolution(FileAnalitic);
	Eiler.writeSolution(FileEiler);
	Heun.writeSolution(FileHeun);

	Analitic.writeEnergy(FileEnergyAnalitic);
	Eiler.writeEnergy(FileEnergyEiler);
	Heun.writeEnergy(FileEnergyHeun);

	return 0;
}

std::string getConfigName(const int argc, const char *argv[])
{
	if (argc < 1)
		throw std::logic_error("There are not Config File Name in args");
	return argv[1];
}

void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &T0, double &X0, double &U0, double &Start, double &Stop, double &DeltaT)
{
	std::ifstream ConfigFile(ConfigFileName);
	nlohmann::json Config = nlohmann::json::parse(ConfigFile);

	W = Config["W"];
	T0 = Config["T0"];
	X0 = Config["X0"];
	U0 = Config["V0"];
	Start = Config["Start"];
	Stop = Config["Stop"];
	DeltaT = Config["Step"];
}
