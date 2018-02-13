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
    Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim);
    double beta;
    double hbar;
    int N; //number of particles
    vec R;
    double E_L;
    double alpha;
    double a_h0;
    vec g;
    double r;
    double r2;
    double A;
    double omega;
    double m;
    vec f;
    double phi;
    int mc; //num MC cycles
    int dim;
    double rho; //position update parameter
    // functions in class
    double PDF(vec R, double alpha_);
    void solve(int mc, int N);// int mc=10, int N=1
    double wavefunc(vec R, double alpha_);
    double Elocal(double omega); // later also R
private:
};
#endif
