#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "armadillo"
using namespace std;
using namespace arma;


int main(){
    clock_t start, end;
    double alpha = 0.5;
    double rho = 0.7;
    // monte carlo cycles
    int mc = 500000;
    // N particles
    int numpart = 10;
    int howmanyDs = 2;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 1e-8;
    // initialize Solver class
    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h);
    //Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h);


    start=clock();
    S.solve();
    end=clock();

    cout<<scientific<<"Local Energy CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    start=clock();
    S.solve_num();
    end=clock();
    cout<<scientific<<"Numerical Energy CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;

    start=clock();
    S.langevin();
    end=clock();
    cout<<scientific<<"Importance Energy CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
}
