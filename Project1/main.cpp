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
    int mc = 100000; // monte carlo cycles
    int numpart = 100; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 1e-3;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    ofstream myfile;
    /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    myfile.open("2316.dat");
    //S.solve(myfile);
    S.solve_num(myfile);
    //S.langevin(myfile);
    myfile.close();
}
