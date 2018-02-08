#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <random>
#include "armadillo"
using namespace std;
using namespace arma;

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
    double r;
    double r2;
    double A;
    vec f;
    int mc; //num MC cycles
    double rho; //position update parameter


private:
}
