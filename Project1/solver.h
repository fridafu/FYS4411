#ifndef SOLVER_H
#define SOLVER_H
#include "armadillo"
#include <random>
#include <iostream>
#include <fstream>

using namespace arma;
using namespace std;

class Solver{
public:
    Solver()
    double beta;
    double hbar;
    int N;
    vec R;
    double E_L;
    double alpha;
    double a_h0;
    vec g;
    vec f;
    int mc; //num MC cycles
    double rho; //position update parameter


private:
}
