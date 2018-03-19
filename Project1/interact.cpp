#include "interact.h"
#include "solver.h"

Interact::Interact(double s_beta,
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
{a = 0.0043;}
mat Interact::init_pos_interact(){
    mat position = init_pos_gaus();
    mat comfort_zone = too_close(position);
    //cout << comfort_zone << endl;
    cout << min(min(distance_part(comfort_zone))) << endl;
    //mat zeropos = zeros(N, dim);
    //mat doesthiswork = too_close(zeropos);
    //cout << "check" << doesthiswork << endl;
    //cout << "is all above a " << distance_part(doesthiswork) << endl;
    return comfort_zone;
    //return doesthiswork;
}
mat Interact::too_close(mat &Rtull){
    random_device rd;
    mt19937_64 genMT64(rd());
    normal_distribution<double> gaussianRNG(0,0.5);
    mat closeness = distance_part(Rtull);
    double sdt = sqrt(dt);
    int counter = 0;
    double a = 0.043;
    int get_away_you_stink = 1;
    while(get_away_you_stink != 0){
        counter += 1;
        get_away_you_stink = 0;
        int o; int p; int dude;
        for(o=0;o<N;o++){

            for(p=o+1;p<N;p++){
                //cout << closeness(o,p) << endl;
                if(closeness(o,p) < a){
                    get_away_you_stink = -1;
                    for(dude=0;dude<dim;dude++){
                        Rtull(o,dude) = gaussianRNG(genMT64)*sdt;
                        //closeness(o,p) = norm(R.row(o)-R.row(p));
                        //cout << Rtull(o,dude)<<endl;
                    }
                }
            }
            closeness = distance_part(Rtull);

        }
    }
    return Rtull;
}
double Interact::wavefunc_interact(mat &R, double alpha_, mat &distanceRij){
    int i; int j;
    double g = 0;
    double f = 1;
    if(dim==1){
        for(i=0;i<N;i++){
            g += R(i)*R(i); // take Product of Pi(g(Ri)
        };
    } else{
        for(i=0;i<N;i++){
            for(j=0;j<dim;j++){
                g += R(i,j)*R(i,j);//g += dot(R.row(i),R.row(i));
                f*=(1 - 1/distanceRij(i,j));
            }
        }
    }

    //double f = (1-(a/distance_part(R));
    double psi = exp(-alpha_*g)*f;
    return psi;
}

void Interact::solve_interact( std::ofstream &myfile){
    myfile << endl << "Calculation with interaction: " << endl;
    start=clock();
    double D = 0.5; //diffusion coefficient
    double Ddt = D*dt;
    double Ddt05 = Ddt*0.5;

    random_device rd;
    mt19937_64 genMT64(rd());
    normal_distribution<double> gaussianRNG(0.,0.5);
    uniform_real_distribution<double> doubleRNG(0,1);

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
        mat R3 = init_pos_interact();
        mat R3new = R3;
        mat distancematrix = distance_part(R3);
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
                double A = greens*(wavefunc_interact(R3new,current_alpha, distancematrix))/wavefunc_interact(R3,current_alpha, distancematrix);
                A *= A;
                // test if new position is more probable than random number between 0 and 1.
                if((A > 1) || (A > doubleRNG(genMT64))){
                    R3(j) = R3new(j); //accept new position
                    Fq(j) = Fqnew(j);
                    accept += 1;
                }else {
                    R3new(j) = R3(j);
                    Fqnew(j) = Fq(j);
                }
                // calculate change in energy
                double deltakinE = energy_interact(R3, current_alpha); // YOU CAN USE energy_num HERE AS WELL. IS THIS RIGHT???
                sumKE += deltakinE;
                }
        }
        num_alpha += 1;
        myfile << scientific << "Acceptance = " << accept/(mc*N) << endl;
    }
    myfile <<scientific << "Kinetic Energy = " << sumKE/(N*mc) << endl;
    end=clock();
    myfile<<scientific<<"Interaction CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    cout << "Interaction and all are finished! Yay." << endl;
}

mat Interact::lapphi(mat &R){
    mat lphi = zeros(N);
    for(int k = 0; k < N; k++){
        if(dim == 3){
            R.col(2) = beta*R.col(2);
            lphi(k) = -2*alpha*(2 + beta -2*alpha*dot(R.row(k),R.row(k))); // write more effecient, calculate 2alpha and beta^2 as own variables
        } else{
            lphi(k) = -2*alpha*((3-dim) + beta -2*alpha*dot(R.row(k),R.row(k)));
        }
    }
    return lphi;
}


mat Interact::nablaphi(mat &R){
    mat newR = R;
    newR.col(2) = R.col(2)*beta;
    return -2*alpha*newR;
}

mat Interact::nablaphinablaF(mat &R, mat &distR){
    mat nphi = -2*alpha*R;
    if(dim == 3){
        nphi = nphi(2)*beta;
    }
    double rkj = 0;
    mat sum = zeros(N);
    int k; int j;
    mat rk;
    mat nphik;
    for(k = 0; k < N; k++){
        nphik = nphi(k);
        rk = R(k);
        for(j = 0; j < N; j++){
            if(k != j){
                rkj = distR(k,j);
                sum(k) += dot(nphik,rk-R(j))*a/(rkj*rkj*(rkj - a));
            }
        }
    }
    return sum;
}

mat Interact::nablaf(mat &R, mat &distR){
    double rkj = 0;
    mat sum = zeros(N,dim);
    int k; int j;
    mat rk;
    for(k = 0; k < N; k++){
        rk = R.row(k);
        for(j = 0; j < N; j++){
            if(k != j){
                rkj = distR(k,j);
                sum.row(k) += (rk-R.row(j))*a/(rkj*rkj*(rkj - a));
            }
        }
    }
    return sum;
}
/*
mat Interact::doublesum(mat &R, mat &distanceR){ // maybe not necessary as this can be stored already from earlier calculations
    mat sumijnotk = nablaf(R,distanceR);
    return dot(sumijnotk,sumijnotk);
}*/

mat Interact::suma2(mat &distanceR){
    int k; int j;
    mat suma = zeros(N);
    double rkj = 0;
    double rkja2 = 0;
    double a2 = a*a;
    for(k = 0; k < N; k++){
        for(j = 0; j < N; j++){
            if(k != j){
                rkj = distanceR(k,j);
                rkja2 = (rkj - a)*(rkj - a);
                suma(k) -= a2/(rkj*rkj*rkja2);
            }
        }
    }
    return suma;
}

double Interact::energy_interact(mat &R, double alphanow){
    mat energy = zeros(N);
    double r2 = 0;
    int i; int j;
    double Vext = 0;
    mat rkj = distance_part(R);
    mat lphi = lapphi(R);
    mat nphi = nablaphi(R);
    mat nf = nablaf(R,rkj);
    mat suma = suma2(rkj);
    double totsum = 0;
    for(int k = 0; k < N; k++){
        energy(k) = lphi(k) + dot(nphi.row(k),nf.row(k)) + dot(nf.row(k),nf.row(k)) + suma(k);
        totsum += energy(k);
    }
    double c = 0.5*m*omega*omega;
    //double Ek = (c - 2*alpha*alpha)*energy + alpha*dim*N; // NYTT UTTRYKK HER
    for(int j = 0; j < N; j++) {
        r2 = 0;
        for(int q = 0; q < dim; q++) {
            r2 += R(j,q)*R(j,q);
        }
        Vext += c*r2; //calculate potential energy
    }
    return totsum + Vext;
}

