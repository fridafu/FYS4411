#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "armadillo"
using namespace std;
using namespace arma;


int main(){
    double a_h0 = 1.;
    double alpha = 1./((2*a_h0)*(2*a_h0));
    double rho = 0.001;
    // monte carlo cycles
    int mc = 10;
    // N particles
    int numpart = 10;
    int howmanyDs = 3;
    // initialize Solver class
    Solver S(1, 1, 1, 1, 1, 0.001, mc, numpart, howmanyDs);
    // Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim)

    for (int i = 0; i < 100; i++){
        S.solve();
    }
    cout << "equilib reached" << endl;

}
