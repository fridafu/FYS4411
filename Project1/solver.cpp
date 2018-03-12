#include "armadillo"
#include "solver.h"

using namespace arma;


void Solver::solve( std::ofstream &myfile){
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
                newE += energy_real(R); // calculate change in energy
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

double Solver::wavefunc(mat &R, double alpha_){// need R, alpha
    //bool interact = y/n
    int i; int j;
    double g = 0;
    if(dim==1){
        for(i=0;i<N;i++){
            g += R(i)*R(i); // take Product of Pi(g(Ri)
        }
    } else{
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                g += R(i,j)*R(i,j);//g += dot(R.row(i),R.row(i));
            }
        }
    }
    double f = 1; //no interaction here!!
    double psi = exp(-alpha_*g)*f;
    //cout << psi << endl;
    return psi;
}

mat Solver::init_pos(){
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

double Solver::PDF(mat &R, double alpha_){
    return pow(abs(wavefunc(R, alpha_)),2);
}

double Solver::energy_local(){
    return 0.5 * hbar * omega * N * dim;
}

//double Solver::rando(){
//    return doubleRNG(genMT64);
//}

void Solver::solve_num( std::ofstream &myfile){
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
    cout << "Numerical energy finished! One to go!!!" << endl;
}

mat Solver::F(mat &R_){
    return -4*R_*alpha;
}

double Solver::energy_real(mat &R){ //done optimization
    int i; int j;
    double energy = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            energy += R(i,j)*R(i,j);
        }
    }
    return (0.5*omega*omega - 2*alpha*alpha)*energy + alpha*dim*N;
}

double Solver::energy_num(mat &R, double alphanow){

    double wavefuncnow = wavefunc(R, alphanow);
    double Ek = 0;
    double Vext = 0;
    double r2 = 0;
    // Calculate kinetic energy by numerical derivation
    mat Rplus = R + h;
    mat Rminus = R - h;

    double c = 0.5*m*omega*omega;
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            r2 += R(j,q)*R(j,q);
        }
        Vext += c*r2; //calculate potential energy
    }
    double wavefuncplus = wavefunc(Rplus, alphanow);
    double wavefuncminus = wavefunc(Rminus, alphanow);
    Ek -= (wavefuncplus+wavefuncminus - 2*wavefuncnow);
    Ek = 0.5 * Ek * h2 / wavefuncnow;
    return Ek + Vext;
}

void Solver::langevin( std::ofstream &myfile){
    myfile << endl << "Importance Sampling:" << endl;
    start=clock();
    double D = 0.5; //diffusion coefficient

    double Ddt = D*dt;
    double Ddt05 = Ddt*0.5;

    random_device rd;
    mt19937_64 genMT64(rd());
    normal_distribution<double> gaussianRNG(0.,0.5);


    // loop over alpha when we try out
    int num_alpha = 0;
    vec alpha_vec = ones(1);
    double current_alpha;
    double sumKE = 0;
    double sdt = sqrt(dt);
    double alpha4 = alpha*(-4);
    while(num_alpha < size(alpha_vec,0)){
        current_alpha = alpha;//alpha_vec(num_alpha);
        // initialize random positions
        mat R3 = init_pos();
        mat R3new = R3;
        int i; int j; int q;
        mat Fq = F(R3);
        //initialize expectation values
        mat R3plus = zeros(N,dim);
        mat R3minus = zeros(N,dim);
        double accept = 0;
        mat Fqnew = Fq;
        double greens;
        // iterate over MC cycles

        for(i=0;i<mc;i++){
            //propose a new position Rnew(boson_j) by moving one boson from position R(boson_j) one at the time
            for(j=0;j<N;j++){
                greens = 0;

                for(q=0;q<dim;q++){
                    R3new(j,q) = R3(j,q) + Ddt*Fq(j,q) + gaussianRNG(genMT64)*sdt;
                    Fqnew(j,q) = alpha4*R3new(j,q);
                    greens += 0.5*(Fq(j,q) + Fqnew(j,q))*(Ddt05*(Fq(j,q)-Fqnew(j,q))-R3new(j,q)+R3(j,q));

                }
                greens = exp(greens);
                double A = greens*(wavefunc(R3new,current_alpha))/wavefunc(R3,current_alpha);
                A *= A;
                // test if new position is more probable than random number between 0 and 1.
                if((A > 1) || (A > gaussianRNG(genMT64))){
                    R3(j) = R3new(j); //accept new position
                    Fq(j) = Fqnew(j);
                    accept += 1;
                }else {
                    R3new(j) = R3(j);
                    Fqnew(j) = Fq(j);
                }
                // calculate change in energy
                double deltakinE = energy_num(R3, current_alpha);
                sumKE += deltakinE;
                }
        }
        num_alpha += 1;
        myfile << scientific << "Acceptance = " << accept/(mc*N) << endl;
    }
    myfile <<scientific << "Kinetic Energy = " << sumKE/(N*mc) << endl;
    end=clock();
    myfile<<scientific<<"Importance sampling CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    cout << "Langevin and all are finished! Yay." << endl;
}

double Solver::absdistance(vec R1, vec R2){
    double r = 0;
    for(int k = 0; k < dim; k++){
        r += pow(R1(k)-R2(k),2);
    }
    return sqrt(r);
}

mat Solver::distance_part(mat &R){
    mat rij = ones(N,N);
    rij.fill(NAN);
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){

            rij(i,j) = norm(R(i)-R(j));
            //rij(i,j) = absdistance(R(i),R(j));
        }
    }
    return rij;
}
/*
mat Solver::init_L2(){
    mat Solver::init_pos_interact(){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> doubleRNG2(0,1);
    int k; int l;
    mat position1 = zeros(N,dim);
    mat position2 = zeros(N,dim);
    double a = 0.043;
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            position1(k,l) = (doubleRNG2(gen) - 0.5)*rho;

            for(int d=0;d<N;d++){
                for(int e=0;d<N;d++){
                    if(norm(position1(e) - position1(d)) > a){

                    } else{

                    }
            }
        }
    }
    return position;
}
*/
Solver::Solver(double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N, int s_dim, double s_h, double s_dt){
    beta = s_beta;
    hbar = s_hbar;
    m = mass;
    omega = s_omega;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = s_alpha;
    rho = s_rho;
    mc = s_mc;
    N = s_N;
    dim = s_dim;
    h = s_h;
    h2 = 1.0/(h*h);
    dt = s_dt;

    //random_device rd;
    //mt19937_64 genMT64(rd);
    //rd.seed(time(NULL));
    //doubleRNG = uniform_real_distribution<double>(0,1);
}
