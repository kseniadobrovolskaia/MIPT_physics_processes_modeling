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

