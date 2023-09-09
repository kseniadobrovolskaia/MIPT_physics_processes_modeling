#include <iostream>
#include <fstream>

void print_Nature_Numbers(std::ostream &Os, unsigned N)
{
        N++;
	unsigned Num = 1;
        do
                Os << Num++ << " ";
        while (N - Num);
        Os << "\n";
}

int main()
{
	std::cout << "Task 3: ";
	print_Nature_Numbers(std::cout, 30);

	std::cout << "Task 4: created file File.txt\n";
	std::ofstream File;
	File.open("File.txt");
	
	print_Nature_Numbers(File, 30);

	File.close();

	std::cout << "Task 5: input number n: ";
	unsigned n;
	std::cin >> n;
	if (std::cin.fail())
		throw std::logic_error("Bad number n. Abort.");

        std::ofstream File1;
        File1.open("File1.txt");

	print_Nature_Numbers(File1, n);

	std::cout << "created file File1.txt\n";

	return 0;
}
