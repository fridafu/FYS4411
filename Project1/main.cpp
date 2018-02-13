#include <iostream>
#include "solver.h"
#include "solver.cpp"
#include "armadillo"
using namespace std;
using namespace arma;


int main()
{
    // vec R = zeros(10);
    // cout << R << endl;
    // return 0;
    double a_h0 = 1.;
    double alpha = 1./((2*a_h0)*(2*a_h0));
    //W = zeros(N,N)
    double rho = 0.001;
    int mc = 10;
    Solver S(1, 1, 1, 1, 1, 0.001, 10, 1,1);//double s_beta, double s_hbar, double mass, double s_omega, double s_alpha, double s_rho, int s_mc, int s_N

    int N = 1;
    //vec alphatest = ones(10);
    //cout<<size(alphatest,0) << endl;
    for (int i = 0; i < 100; i++){
        S.solve(mc, N);
    }
    cout << "equilib reached" << endl;

}
