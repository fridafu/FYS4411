#include "armadillo"
#include "solver.h"

using namespace arma;


void Solver::solve(){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double energySum = 0;
    double energySquaredSum = 0;
    double energynewSum = 0;

    while(num_alpha < size(alpha_vec,0)){
        current_alpha = alpha;//alpha_vec(num_alpha);
        // initialize random positions
        mat R = init_pos();
        mat Rnew = R;
        int i; int j; int q;

        //initialize expectation values
        mat Rplus = zeros(N,dim);
        mat Rminus = zeros(N,dim);
        double accept = 0;

        // iterate over MC cycles
        for(i=0;i<mc;i++){
            //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
            for(j=0;j<N;j++){

                for(q=0;q<dim;q++){
                    r = doubleRNG(gen) - 0.5;
                    Rnew(j,q) = R(j,q) + r*rho;
                }

                A = (wavefunc(Rnew,current_alpha))/wavefunc(R,current_alpha);
                A *= A;
                // test if new position is more probable or if larger than random number doubleRNG(gen) in (0,1)
                if((A>1)||(A > doubleRNG(gen))){
                    //accept new position
                    R(j) = Rnew(j);
                    accept += 1;
                }
                // calculate change in energy
                double deltaE = energy_local();
                energynewSum += energy_real(R);
                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;
                }
        }
        num_alpha += 1;
        cout << "accept " << accept/(mc*N) << endl;
    }


    double energy = energySum/(mc * N);
    double newenergy = energynewSum/(mc*N);
    double totalenergy = energySum/mc;
    double energySquared = energySquaredSum/(mc * N);

    cout << "E_tot = " << totalenergy << endl;
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    cout << "New try on Energy: " << newenergy << endl;

}

double Solver::wavefunc(mat R, double alpha_){// need R, alpha
    //bool interact = y/n
    int i;
    g = 0;
    if(dim==1){
        for(i=0;i<N;i++){
            g += -alpha_*dot(R.row(i),R.row(i)); // take Product of Pi(g(Ri)
        }
    }
    else{
        for(i=0;i<N;i++){
            g += -alpha_*dot(R.row(i),R.row(i));
        }
    }
    f = 1; //no interaction here!!
    psi = exp(g)*f;
    return psi;
}

mat Solver::init_pos(){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG2(-1,1);
    int k; int l;
    mat position = zeros(N,dim);
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            position(k,l) = doubleRNG2(gen)*rho;
        }
    }
    return position;
}

double Solver::PDF(mat R, double alpha_){
    return pow(abs(wavefunc(R, alpha_)),2);
}

double Solver::energy_local(){
    return 0.5 * hbar * omega * N * dim;
}

double Solver::energy_real(mat R){


    // Calculate kinetic energy by numerical derivation
    /*
    else{
        for(j = 0; j < N; j++) {
            r2 = 0;

            for(q = 0; q < dim; q++) {
                r2 += R(j,q)*R(j,q);

                energy += alpha;
        }
            energy += -2*alpha*alpha*r2  + 0.5*omega*omega*r2;
        //calculate potential energy

    }
    }
    */
    int i;
    double r2 = 0;
    for(i=0;i < N;i++){
        r2 = dot(R.row(i), R.row(i));
    }

    return (0.5*omega*omega-2*alpha*alpha)*r2 + N*alpha*dim;
}
void Solver::solve_num(){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double sumKE = 0;

    while(num_alpha < size(alpha_vec,0)){
        current_alpha = alpha;//alpha_vec(num_alpha);
        // initialize random positions
        mat R = init_pos();
        mat Rnew = R;
        int i; int j; int q;

        //initialize expectation values
        mat Rplus = zeros(N,dim);
        mat Rminus = zeros(N,dim);
        double accept = 0;

        // iterate over MC cycles
        for(i=0;i<mc;i++){
            //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
            for(j=0;j<N;j++){
                for(q=0;q<dim;q++){
                    r = doubleRNG(gen) - 0.5;
                    Rnew(j,q) = R(j,q) + r*rho;
                }

                A = (wavefunc(Rnew,current_alpha))/wavefunc(R,current_alpha);
                A *= A;
                // test if new position is more probable or if larger than random number doubleRNG(gen) in (0,1)
                if((A > 1) || (A > doubleRNG(gen))){
                    //accept new position
                    R(j) = Rnew(j);
                    accept += 1;
                }
                // calculate change in energy
                double deltakinE = energy_num(R, current_alpha);
                sumKE += deltakinE;
                // calculate total energy
                }
        }
        num_alpha += 1;
        cout << "accept " << accept/(mc*N) << endl;
    }

    cout << "sumKE+Vext (should be equal to Energy) " << sumKE/(N*mc) << endl;
}

double Solver::energy_num(mat R, double alphanow){
    double wavefuncnow = wavefunc(R, alphanow);

    mat Rplus;
    mat Rminus;
    double Ek = 0;
    double Vext = 0;
    double r2 = 0;
    // Calculate kinetic energy by numerical derivation
    Rplus = Rminus = R;
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            Rplus(j,q) += h;
            Rminus(j,q) -= h;
            r2 += R(j,q)*R(j,q);
        }
        //calculate potential energy
        Vext += 0.5*m*omega*omega*r2;
    }
    wavefuncplus = wavefunc(Rplus, alphanow);
    wavefuncminus = wavefunc(Rminus, alphanow);
    Ek -= (wavefuncplus+wavefuncminus - 2*wavefuncnow);
    Ek = 0.5 * Ek * h2 / wavefuncnow;
    return Ek + Vext;
}


Solver::Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h){
    beta = s_beta;
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = s_alpha;
    rho = s_rho;// 0.001;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
    h = s_h;
    h2 = 1.0/(h*h);
}
