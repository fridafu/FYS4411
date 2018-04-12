#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <random>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

class Solver{
public:
    clock_t start, end;
    Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    double beta;
    double hbar;
    int N; //number of particles
    double E_L;
    double alpha;
    double a_h0;
    double B;
    double omega;
    double m;
    double h;
    double h2;
    double dt;
    int mc; //num MC cycles
    int dim;
    double rho; //position update parameter

    // functions in class
    double wavefunc(const mat &R, double alpha_);
    double d_wavefunc(const mat &R, double alpha_);
    mat init_pos_gaus();
    mat distance_part(const mat &R);
    double energy_num(const mat &R, double alpha);
    mat F(const mat &R_, double alpha_);
    double energy_real(const mat &R, double alpha); //
private:
};
#endif
