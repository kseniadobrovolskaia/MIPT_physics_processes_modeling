#include <iostream>
#include <fstream>
#include <array>


unsigned long long Fib(unsigned N) {
    unsigned long long First = 0ull, Second = 1ull;
    int Idx;
    if (N == 0)
	   return 0ull;
    for (Idx = 2; Idx <= N; ++Idx)
    {
        unsigned long long Tmp = Second;
        Second = Second + First;
        First = Tmp;
    }
    return Second;
}


int main()
{
	std::ofstream File;
	File.open("File.txt");
	std::cout << "Task 7: created file File.txt\n";
	const unsigned Cnt = 10;
	unsigned Num = 0;
	std::array<unsigned, Cnt> Array;

	while (Cnt - Num++)
	{
		Array[Num - 1] = Fib(Num);
		File << Num << ": " << Array[Num - 1]  << "\n";
	}
	File << "\n";
	File.close();

	std::cout << "Task 8: ";
	for (auto Elem : Array)
		std::cout << Elem << " ";
	std::cout << "\n";

	return 0;
}
