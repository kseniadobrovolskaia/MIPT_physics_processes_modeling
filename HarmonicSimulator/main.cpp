#include "Solver.hpp"


int main()
{
	double W, X0, U0, Start, Stop, DeltaT = 0.01;
	const char *FileNameAnalitic = "Analitic.bin";
	const char *FileNameEiler = "Eiler.bin";
	const char *FileNameHeun = "Heun.bin";
	std::ofstream FileWithAnalitic(FileNameAnalitic, std::ios::binary);
	std::ofstream FileWithEiler(FileNameEiler, std::ios::binary);
	std::ofstream FileWithHeun(FileNameHeun, std::ios::binary);
	
	std::cout << "Дайте частоту : ";
	std::cin >> W;
	std::cout << "Дайте начальную координату и скорость: X0 = ";
	std::cin >> X0;
	std::cout << "U0 = ";
	std::cin >> U0;

	HarmonicEquation Oscilliator(W);

	AnalyticalSolver Analitic(Oscilliator);
	Analitic.setStartConditions(X0, U0);

	EilerSolver Eiler(Oscilliator, DeltaT);
	Eiler.setStartConditions(X0, U0);

	HeunSolver Heun(Oscilliator, DeltaT);
	Heun.setStartConditions(X0, U0);

	std::cout << "Промежуток времени: Tstart = ";
	std::cin >> Start;
	std::cout << "Tstop = ";
	std::cin >> Stop;

	Analitic.writeSolution(Start, Stop, DeltaT, FileWithAnalitic);
	Eiler.writeSolution(Start, Stop, DeltaT, FileWithEiler);
	Heun.writeSolution(Start, Stop, DeltaT, FileWithHeun);

	return 0;
}
