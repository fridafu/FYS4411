#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "bruteforce.h"
#include "bruteforce.cpp"
#include "impsamp.h"
#include "impsamp.cpp"
#include "interact.h"
#include "interact.cpp"
#include "armadillo"
using namespace std;
using namespace arma;

int main(){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    double h = 0.001;
    int mc = 1000; // monte carlo cycles
    int numpart = 100; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

    //ofstream myfile;
    //myfile.open("classy.dat");     /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    //B->solve(myfile);
    //B->solve_num(myfile);
    //Imp->langevin(myfile);
    //Int->solve_interact(myfile);
    //myfile.close();
    Imp->best_alpha();
    delete B;
    delete Imp;
    delete Int;
}
