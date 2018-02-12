#include "armadillo"
#include "solver.h"

using namespace arma;


void Solver::solve(int mc, int N){// need mc, N

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha
    int num_alpha = 0;

    vec alpha = ones(1);
    double pdf;
    double current_alpha;
    double energySum = 0;
    double energySquaredSum = 0;
    //nt len_alpha = size(alpha);
    while(num_alpha < size(alpha,0)){
        current_alpha = alpha(num_alpha);
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
           pdf = PDF(R,current_alpha);
           vec Rnew;

           //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
           for(j=0;j<N;j++){
                Rnew = zeros(N); // to suggest new position for boson
                r = doubleRNG(gen) - 0.5;
                Rnew(j) = R(j) + r*rho;
                pdf = PDF(R,current_alpha); // now we calculate it two times... not necessary
                // how to determine if we accept or reject new position
                A = (wavefunc(Rnew,current_alpha))/wavefunc(R,current_alpha);
                if(A > 1 || A > doubleRNG(gen)){
                    //accept new position
                    R(j) = Rnew(j);
                }
                //else{
                    //compare probability with a random number between 0 and 1
                    //r2 = doubleRNG(gen);
                    //if(A > r2){
                    //    //accept new position
                    //    R(j) = Rnew(j);
                    //}
                //}


           }
        double deltaE = Elocal(omega);
        energySum += deltaE;
        energySquaredSum += deltaE*deltaE;
        }
        num_alpha += 1;
    }
    double energy = energySum/(mc * N);
    double energySquared = energySquaredSum/(mc * N);
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
}

double Solver::wavefunc(vec R, double alpha_){// need R, alpha
    //bool interact = y/n
    //g = exp(-alpha*(x(i)^2 + y(i)^2 + beta*z(i)^2));
    g = exp(-alpha_*dot(R,R));
    f = ones(1); //no interaction here!!
    phi = dot(g,f);
    return phi;
}


double Solver::PDF(vec R, double alpha_){
    return pow(abs(wavefunc(R, alpha_)),2);

}
double Solver::Elocal(double omega){
    // cout << "in Elocal"<< endl;
    cout << "hbar " << hbar << endl;// * omega * N
    cout << "omega " << omega << endl;
    cout << "N " << N << endl;
    return 0.5 * hbar * omega * N;
}

Solver::Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim){
    beta = s_beta;
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = 1./((2*a_h0)*(2*a_h0));
    //W = zeros(N,N)
    rho = s_rho;// 0.001;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
}
