#include "Solver.hpp"


int main()
{
	double W, X0, U0, Start, Stop, DeltaT;
	const char *FileName = "Solution.bin";
	std::ofstream FileWithSolution(FileName, std::ios::binary);
	std::cout << "Дайте частоту : ";
	std::cin >> W;
	std::cout << "Дайте начальную координату и скорость: X0 = ";
	std::cin >> X0;
	std::cout << "U0 = ";
	std::cin >> U0;

	HarmonicEquation Oscilliator(W);
	AnalyticalSolver Analitic(Oscilliator);
	Analitic.setStartConditions(X0, U0);

	std::cout << "Промежуток времени: Tstart = ";
	std::cin >> Start;
	std::cout << "Tstop = ";
	std::cin >> Stop;
	DeltaT = 0.1;
	auto [X, Y] = Analitic.getXV(Start);
	auto [X2, Y2] = Analitic.getXV(Stop);

	Analitic.writeSolution(Start, Stop, DeltaT, FileWithSolution);

	std::cout << "Решение в Tstart : " << X << " " << Y << "\n";
	std::cout << "Решение в Tstop : " << X2 << " " << Y2 << "\n";
	std::cout << "Решение на промежутке записано в файл " << FileName << "\n";
	return 0;
}
