#include "solver.h"

Solver::Solver(double s_beta,
               double s_hbar,
               double mass,
               double s_omega,
               double s_alpha,
               double s_rho,
               int s_mc,
               int s_N,
               int s_dim,
               double s_h,
               double s_dt){
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


double Solver::wavefunc(mat &R, double alpha_){
    //bool interact = y/n ??
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
    return psi;
}

double Solver::d_wavefunc(mat &R, double alpha_){
    //bool interact = y/n ??
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
    return psi*(-g);
}
mat Solver::init_pos_gaus(){
    random_device rd;
    mt19937_64 genMT64(rd());
    normal_distribution<double> gaussianRNG(0.,1);
    double sdt = sqrt(dt);
    int k; int l;
    mat position = zeros(N,dim);
    for(k=0;k<N;k++){
        for(l=0;l<dim;l++){
            position(k,l) = gaussianRNG(genMT64)*sdt;
        }
    }
    return position;
}
mat Solver::distance_part(mat R){
    mat rij = zeros(N,N);
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            rij(i,j) = norm(R.row(i) - R.row(j));
            rij(j,i) = rij(i,j);

            //rij(i,j) = absdistance(R(i),R(j));
        }
    }
    return rij;
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

mat Solver::F(mat &R_, double alpha_){
    return -4*R_*alpha_;
}

double Solver::energy_real(mat &R, double alpha){ //done optimization
    int i; int j;
    double energy = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            energy += R(i,j)*R(i,j);
        }
    }
    double en = (0.5*omega*omega - 2*alpha*alpha)*energy + alpha*dim*N;
    return en;

}
