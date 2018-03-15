#ifndef BRUTEFORCE_H
#define BRUTEFORCE_H
#include "solver.h"

class Bruteforce : public Solver {

public:
    Bruteforce(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    double energy_local();
    mat init_pos();
    double energy_real(mat &R); // not analytical solution
    void solve(std::ofstream &myfile);// int mc=10, int N=1
    void solve_num(std::ofstream &myfile);
};


#endif
