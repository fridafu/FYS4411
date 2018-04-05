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
#include "catch.hpp"

using namespace std;
using namespace arma;

int main(){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001; // [0.001, 0.01]
    double h = 0.001;
    int mc = 1048576; // monte carlo cycles

    int numpart = 10; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!

    int howmanyDs = 1;
    double beta = 1;//sqrt(8);
    double hbar = 1;
    double mass = 1;
    double omega = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;


    myfile.open("blocktest.dat");     /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    myfile2.open("Ea.dat");
    myfile3.open("En.dat");
    myfile4.open("Eimp.dat");
    B->solve(myfile,myfile2);
    B->solve_num(myfile,myfile3);
    Imp->langevin(myfile,myfile4,alpha);
    //Int->solve_interact(myfile, alpha);
    myfile.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    //Imp->best_alpha();
    delete B;
    delete Imp;
    delete Int;
}
