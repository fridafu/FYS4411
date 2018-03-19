#ifndef INTERACT.H
#define INTERACT.H
#include "solver.h"

class Interact : public Solver{
public:
    Interact(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    mat init_pos_interact();
    mat distancematrix;
    mat too_close(mat &Rtull);
    double wavefunc_interact(mat &R, double alphanow, mat &distancematrix);
    double energy_interact(mat &R, double alphanow);
    void solve_interact(std::ofstream &myfile);
};
#endif
