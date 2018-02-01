#include solver.h


void Solver::solve(n=10,N=1){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> doubleRNG(0,1);

    // loop over alpha
    int num_alpha = 0;

    while(num_alpha < len(alpha)){
        // initialize random positions
        i = rand() % N;
        j = rand() % N;
        r = doubleRNG(gen);
        //initialize expectation values
        int expEL = 0;
        int expEL2 = 0;

        // iterate over MC cycles
        for(i=0;i<n;i++){
           //set up PDF |phi|^2
           pdf = PDF(r)
           //propose a new position R by moving one boson at the time

        }


        num_alpha += 1
    }
}

void wavefunc(R,alpha){
    bool interact = //y/n;
    beta = 1;
    hbar = 1;
    m = 1;
    omega = 1;
    a_h0 = sqrt(hbar/(m*omega));
    alpha = 1./(2*a_h0^2);
    g = exp(-alpha(x(i)^2 + y(i)^2 + beta*z(i)^2));
    f = 1;
    phi = g*f;
}

void PDF(vec R){
    PDF = abs(wavefunc(R))^2
}
