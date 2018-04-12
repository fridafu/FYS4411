#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "interact.h"
#include "catch.hpp"

using namespace std;
using namespace arma;

int main(){
    vec alphavec = linspace<vec>(0.3, 0.7, 9);
    double rho = 0.1;
    double dt = 0.001; // [0.001, 0.01]
    double h = 0.0001;
    int numpart = 100; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = 1e5 / numpart; // monte carlo cycles
    int howmanyDs = 3;
    double beta = sqrt(8);
    double hbar = 1;
    double mass = 1;
    double omega = 1;

    //for(int elem=0; elem<size(alphavec,0); elem++){
    double alpha = 0.5; //alphavec(elem);
    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

    ofstream myfile;
    ofstream myfile2;
    ofstream myfile3;
    ofstream myfile4;
    ofstream myfile5;

    string filename = "b_n100_d1";

    myfile.open(filename + "_alpha" + std::to_string(alpha) + ".dat");     /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    //myfile2.open("Ea" + filename + "_alpha"  + std::to_string(alpha) + ".dat");
    //myfile3.open("En" + filename + "_alpha" + std::to_string(alpha) + ".dat");
    //myfile4.open("Eimp" + filename + "_alpha" + std::to_string(alpha) + ".dat");
    myfile5.open("Eint" + filename + "_alpha" + std::to_string(alpha) + ".dat");
    //B->solve(myfile,myfile2);
    //B->solve_num(myfile,myfile3);
    //Imp->langevin(myfile,myfile4,alpha);
    Int->solve_interact(myfile, myfile5, alpha);
    myfile.close();
    //myfile2.close();
    //myfile3.close();
    //myfile4.close();
    myfile5.close();

    //Imp->best_alpha();
    delete Int;
    delete Imp;
    delete B;
    //}
}
