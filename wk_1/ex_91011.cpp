#include <iostream>
#include <array>
#include <vector>
#include <fstream>

double garmonic(unsigned N)
{
	return 1.0 / N;
}

int main()
{
	const unsigned N = 10;
	unsigned Num = 0;
	std::array<double, N> Array;
	std::cout << std::scientific;
	std::cout << "Task 9: ";

	while (N - Num++)
	{
		Array[Num - 1] = garmonic(Num);
		std::cout << Array[Num - 1] << " ";
	}
	
	std::vector<float> Vector;
	Vector.reserve(N);

	std::cout << "\nTask 10: use vector\n";
	for (auto It = Array.rbegin(); It != Array.rend(); It++)
		Vector.push_back(*It);
	
	std::ofstream BinFile;
	BinFile.open("BinFile.txt");
	std::cout << "Task 11: created file BinFile.txt\n";

	BinFile.write((char *)&Vector[0], Vector.size() * sizeof(Vector[0]));

	BinFile.close();
	return 0;
}


