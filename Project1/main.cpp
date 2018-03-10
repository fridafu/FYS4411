#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "armadillo"
using namespace std;
using namespace arma;


int main(){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;

    // monte carlo cycles
    int mc = 1000000;
    // N particles
    int numpart = 100; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 1e-8;
    // initialize Solver class
    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    //Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h);
    ofstream myfile;

    /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    myfile.open("RtilR1wavefunconlylangvin.dat");
    //S.solve(myfile);
    //S.solve_num(myfile);
    S.langevin(myfile);
    myfile.close();

}
