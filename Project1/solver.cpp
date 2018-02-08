#include "armadillo"
#include "solver.h"

using namespace arma;


void Solver::solve(mc=10,N=1){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha
    int num_alpha = 0;
    int dim = 1;
    while(num_alpha < len(alpha)){
        // initialize random positions
        int i; int j;
        //i = rand() % N;
        //j = rand() % N;

        //initialize expectation values
        int expEL = 0;
        int expEL2 = 0;
        vec R = zeros((N,dim));
        vec Rnew = R;
        // call for function that initialize position when we are ready

        // iterate over MC cycles
        for(i=0;i<mc;i++){
           //set up PDF |phi|^2
           pdf = PDF(R);
           vec Rnew;

           //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
                for(j=0;j<N;j++){
                    Rnew = zeros(N); // to suggest new position for boson
                    r = doubleRNG(gen) - 0.5;
                    Rnew(j) = R(j) + r*rho;
                    pdf = PDF(R);
                    // how to determine if we accept or reject new position



           }
        }
        num_alpha += 1
    }
}

void wavefunc(R,alpha){
    //bool interact = y/n
    //g = exp(-alpha*(x(i)^2 + y(i)^2 + beta*z(i)^2));
    g = exp(-alpha*dot(R,R));
    f = 1; //no interaction here!!
    phi = g*f;
}

void PDF(vec R){
    return abs(wavefunc(R))^2;
}

Solver::Solver(double beta, double hbar, double m, double omega, double a_h0, double alpha, double rho){
    beta = 1;
    hbar = 1;
    m = 1;
    omega = 1;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = 1./(2*a_h0^2);
    //W = zeros(N,N)
    rho = 0.001;
}
