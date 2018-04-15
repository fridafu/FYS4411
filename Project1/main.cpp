#include <iostream>
#include "solver.h"
#include "bruteforce.h"
#include "impsamp.h"
#include "interact.h"
//#include "catch.hpp"

using namespace std;
using namespace arma;

int main(){
    vec alphavec = linspace<vec>(0.3, 0.7, 9);
    vec dtvec = logspace<vec>(-7,1,17);
    cout << dtvec << endl;
    double rho = 0.1;
    //double dt = 0.01; // [0.001, 0.01]
    double h = 0.0001;
    int numpart = 10; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int mc = (1048576 + 1000) / numpart; // monte carlo cycles
    int howmanyDs = 3;
    double beta = 1;//sqrt(8);
    double hbar = 1;
    double mass = 1;
    double omega = 1;

    ofstream myfile;

    string folder = "alpha0_55/";
    myfile.open(folder + "dt_table_test.dat");

    for(int elem=0; elem<size(dtvec,0); elem++){
        double alpha = 0.55; //alphavec(elem);
        double dt = dtvec(elem); // 0.0000001;//
        cout << "Doing " << scientific << to_string(dt) << endl;
        Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
        Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
        Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
        Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

        ofstream myfile2;
        ofstream myfile3;
        ofstream myfile4;
        ofstream myfile5;
        string filename =folder + "c_n" + std::to_string(numpart) + "_d" + std::to_string(howmanyDs) + "_mc" + std::to_string(mc);
        //myfile.open(filename + "dt" + std::to_string(dt) + ".dat");     //CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //myfile2.open("Ea_" + filename + "_alpha"  + std::to_string(alpha) + ".dat");
        //myfile3.open("En_" + filename + "_alpha" + std::to_string(alpha) + ".dat");
        //myfile4.open("Eimp_" + filename + "dt" + std::to_string(dt) + ".dat");
        myfile4.open(folder+std::to_string(log10(dt)) + ".dat");
        //myfile5.open("Eint_" + filename + "_alpha" + std::to_string(alpha) + ".dat");
        //B->solve(myfile,myfile2);
        //B->solve_num(myfile,myfile3);
        Imp->langevin(myfile,myfile4,alpha);
        //Int->solve_interact(myfile, myfile5, alpha);
        //myfile.close();
        //myfile2.close();
        //myfile3.close();
        myfile4.close();
        //myfile5.close();
        //Imp->best_alpha();
        delete Int;
        delete Imp;
        delete B;
    }
    myfile.close();
}
