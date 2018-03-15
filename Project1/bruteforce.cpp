#include "armadillo"
#include "solver.h"
#include "bruteforce.h"

using namespace arma;

Bruteforce::Bruteforce(double s_beta,
                       double s_hbar,
                       double mass,
                       double s_omega,
                       double s_alpha,
                       double s_rho,
                       int s_mc,
                       int s_N,
                       int s_dim,
                       double s_h,
                       double s_dt)
:
    Solver(s_beta, s_hbar, mass,s_omega, s_alpha, s_rho, s_mc, s_N, s_dim, s_h, s_dt)
{}

double Bruteforce::energy_local(){
    return 0.5 * hbar * omega * N * dim;
}

void Bruteforce::solve( std::ofstream &myfile){
    double energy = energy_local();

    myfile << "dim = " << dim << ", N = " << N << ", dt = " << dt << ", alpha = " << alpha << " and mc = " << mc << endl << endl;
    myfile << scientific << "Theoretical Energy = " << energy << endl << endl;
    myfile << "Brute force:" << endl;

    start=clock();
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double newE = 0;
    mat testR = init_pos();
    mat dis = distance_part(testR);
    cout << testR << endl;
    cout << dis << endl;
    while(num_alpha < size(alpha_vec,0)){
        current_alpha = alpha;

        mat R = init_pos(); // initialize random positions
        mat Rnew = R;
        int i; int j; int q;

        double accept = 0;

        for(i=0;i<mc;i++){ // iterate over MC cycles
            //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
            for(j=0;j<N;j++){
                for(q=0;q<dim;q++){
                    Rnew(j,q) = R(j,q) + (doubleRNG(genMT64) - 0.5)*rho;
                }

                double A = (wavefunc(Rnew,current_alpha))/wavefunc(R, current_alpha);
                A *= A;

                // test if new position is more probable than random number between 0 and 1.
                if((A > 1) || (A > doubleRNG(genMT64))){
                    R(j) = Rnew(j); //accept new position
                    accept += 1;
                } else {
                    Rnew(j) = R(j);
                }
                newE += energy_real(R, current_alpha); // calculate change in energy
           }
        }

        num_alpha += 1;
        myfile << scientific << "Acceptance = " << accept/(mc*N) << endl;
        cout << "Brute force finished! Hang in there <3" << endl;
    }

    /*
    double energy = energySum/(mc * N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);
    */
    end=clock();
    myfile << scientific << "Calculated energy =  " << newE/(mc*N) << endl;
    myfile<< scientific <<"Brute Force CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
}

mat Bruteforce::init_pos(){
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(0,1);
    int k; int l;
    mat position = zeros(N,dim);
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            //position(k,l) = (rando() - 0.5)*rho;
            position(k,l) = (doubleRNG(genMT64) - 0.5)*rho;
        }
    }
    return position;
}
void Bruteforce::solve_num( std::ofstream &myfile){
    random_device rd;
    mt19937_64 genMT64(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    myfile << endl << "Numerical derivation of kinetic energy:" << endl;
    start=clock();

    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double sumKE = 0;

    while(num_alpha < size(alpha_vec,0)){
        current_alpha = alpha;//alpha_vec(num_alpha);
        // initialize random positions
        mat R2 = init_pos();
        mat R2new = R2;
        int i; int j; int q;

        //initialize expectation values
        double accept = 0;

        // iterate over MC cycles
        for(i=0;i<mc;i++){
            //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
            for(j=0;j<N;j++){
                for(q=0;q<dim;q++){
                    //R2new(j,q) = R2(j,q) + (rando() - 0.5)*rho;
                    //cout << R2new(j,q) << endl;
                    R2new(j,q) = R2(j,q) + (doubleRNG(genMT64) - 0.5)*rho;

                }

                double A = (wavefunc(R2new,current_alpha))/wavefunc(R2,current_alpha);

                A *= A;
                // test if new position is more probable than random number between 0 and 1.
                if((A > 1) || (A > doubleRNG(genMT64))) {

                    R2(j) = R2new(j); //accept new position
                    accept += 1;
                } else {
                    R2new(j) = R2(j);
                }
                // calculate change in energy
                double drit = energy_num(R2, current_alpha);
                // cout << R2 << endl;

                //sumKE += energy_num(R2, current_alpha);
                sumKE += drit;
            }
        }
        num_alpha += 1;
        myfile << "Acceptance =" << accept/(mc*N) << endl;
    }

    myfile << scientific << "Numerical Energy = " << sumKE/(N*mc) << endl;
    end=clock();
    myfile<<scientific<<"Num dev CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    cout << "Numerical energy finished!!!" << endl;
}
