#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "armadillo"
using namespace std;
using namespace arma;


int main()
{
    vec R = zeros(10);
    cout << R << endl;
    return 0;

    Solver S();

    for (int i = 0; i < 100; i++)
    {
         S.solve(mc);
    }

    cout << "equilib reached" << endl;

}
