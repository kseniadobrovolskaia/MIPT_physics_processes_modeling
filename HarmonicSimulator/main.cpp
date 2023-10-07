#include "Solver.hpp"
#include "json.hpp"



std::string getConfigName(const int argc, const char *argv[]);
void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &X0, double &U0, double &Start, double &Stop, double &DeltaT);


int main(const int argc, const char *argv[])
{
	double W, X0, U0, Start, Stop, DeltaT;
	const char *FileNameAnalitic = "Analitic.bin";
	const char *FileNameEiler = "Eiler.bin";
	const char *FileNameHeun = "Heun.bin";
	std::ofstream FileWithAnalitic(FileNameAnalitic, std::ios::binary);
	std::ofstream FileWithEiler(FileNameEiler, std::ios::binary);
	std::ofstream FileWithHeun(FileNameHeun, std::ios::binary);

	const std::string ConfigFileName = getConfigName(argc, argv);
	getStartConditionsFromConfigFile(ConfigFileName, W, X0, U0, Start, Stop, DeltaT);

	HarmonicEquation Oscilliator(W);

	AnalyticalSolver Analitic(Oscilliator);
	Analitic.setStartConditions(X0, U0);

	EilerSolver Eiler(Oscilliator, DeltaT);
	Eiler.setStartConditions(X0, U0);

	HeunSolver Heun(Oscilliator, DeltaT);
	Heun.setStartConditions(X0, U0);

	Analitic.writeSolution(Start, Stop, DeltaT, FileWithAnalitic);
	Eiler.writeSolution(Start, Stop, DeltaT, FileWithEiler);
	Heun.writeSolution(Start, Stop, DeltaT, FileWithHeun);

	return 0;
}


std::string getConfigName(const int argc, const char *argv[])
{
	if (argc < 1)
	{
		throw std::logic_error("There are not Config File Name in args");
	}
	return argv[1];
}


void getStartConditionsFromConfigFile(const std::string ConfigFileName, double &W, double &X0, double &U0, double &Start, double &Stop, double &DeltaT)
{
	std::ifstream ConfigFile(ConfigFileName);
	nlohmann::json Config = nlohmann::json::parse(ConfigFile);

	W = Config["W"];
	X0 = Config["X0"];
	U0 = Config["V0"];
	Start = Config["Start"];
	Stop = Config["Stop"];
	DeltaT = Config["Step"];
}
