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
    Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h);
    double beta;
    double hbar;
    int N; //number of particles
    mat R;
    double E_L;
    double alpha;
    double a_h0;
    double g;
    double r;
    double r2;
    double A;
    double omega;
    double m;
    double f;
    double psi;
    double h;
    double h2;
    double wavefuncplus;
    double wavefuncminus;
    int mc; //num MC cycles
    int dim;
    double rho; //position update parameter
    // functions in class
    double PDF(mat R, double alpha_);
    void solve();// int mc=10, int N=1
    void solve_num();
    double wavefunc(mat R, double alpha_);
    double energy_local(double omega); // later also R
    mat init_pos();
    double energy_num(mat R, double alphanow);
private:
};
#endif
