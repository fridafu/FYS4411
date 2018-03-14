#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <fstream>
#include <random>
#include "armadillo"
#include <cmath>
using namespace std;
using namespace arma;

class Solver{
public:
    Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    double beta;
    double hbar;
    int N; //number of particles
    double E_L;
    double alpha;
    double a_h0;
    //double g;
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
    double PDF(mat &R, double alpha_);
    void solve(std::ofstream &myfile);// int mc=10, int N=1
    void solve_num(std::ofstream &myfile);
    void langevin(std::ofstream &myfile);
    mat F(mat &R_);
    double wavefunc(mat &R, double alpha_);
    double energy_local(); // later also R
    mat init_pos();
    mat init_pos_interact();
    double energy_real(mat &R); // not analytical solution
    double energy_num(mat &R, double alphanow);
    clock_t start, end;
private:
};
#endif
