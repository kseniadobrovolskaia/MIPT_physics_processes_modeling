#include <iostream>


void print_Nature_Numbers(std::ostream &Os, unsigned N);


int main(int argc, const char *argv[])
{
	std::string Number = argv[1];
        unsigned N;

	if (!(N = std::stoi(Number)))
		throw std::logic_error("Bad number n. Abort.");	
	print_Nature_Numbers(std::cout, N);
	return 0;
}
