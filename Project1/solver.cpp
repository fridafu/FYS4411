#include solver.h


void Solver::solve(mc=10,N=1){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha
    int num_alpha = 0;

    while(num_alpha < len(alpha)){
        // initialize random positions

        i = rand() % N;
        j = rand() % N;

        //initialize expectation values
        int expEL = 0;
        int expEL2 = 0;

        // iterate over MC cycles
        for(i=0;i<mc;i++){
           //set up PDF |phi|^2
           pdf = PDF(R);
           R(0) = 0;
           //propose a new position R by moving one boson at the time
               for(j=0;j<N;j++){
                   r = doubleRNG(gen) - 0.5;
                   R(j+1) = R(j) + r*rho;
                   pdf = PDF(R);
           }
        }
        num_alpha += 1
    }
}

void wavefunc(R,alpha){
    //bool interact = y/n
    //g = exp(-alpha*(x(i)^2 + y(i)^2 + beta*z(i)^2));
    g = exp(-alpha*R)
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
