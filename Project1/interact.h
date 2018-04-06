#ifndef INTERACT.H
#define INTERACT.H
#include "solver.h"

class Interact : public Solver{
public:
    double a;
    Interact(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt);
    mat init_pos_interact();
    mat too_close(mat Rtull);
    double wavefunc_interact(mat R, double alphanow, mat distancematrix);
    double energy_interact(mat R, double alphanow);
    vec solve_interact(std::ofstream &myfile, ofstream &myfile5, double alphanoe);
    double d_wavefunc_interact(mat R, double alphanow, mat distancematrix);

    mat lapphi(mat R, double alpha_);
    mat nablaphi(mat R, double alpha_);
    mat nablaphinablaF(mat R, mat distR, double alpha_);
    mat nablaf(mat R, mat distR);
    mat doublesum(mat R, mat distanceR);
    mat suma2(mat distanceR, mat R);
    mat quantumF(mat R, double alpha_, mat rij);
    mat nablafsquared(mat init_distance, mat R);
    mat newnablaf(mat init_distance, mat R);

};
#endif
