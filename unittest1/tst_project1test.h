#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "../Project1/solver.h"
#include "../Project1/bruteforce.h"

using namespace testing;


TEST(project1test, wavefunction)
{
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    int mc = 10000; // monte carlo cycles
    int numpart = 500;
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 0.001;
    double m = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
//  mat R = ones(numpart,howmanyDs)*rho;//S.init_pos();
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    mat R = B->init_pos();

    double g = 1;
    for(int i = 0; i < numpart; i++){
        g *= exp(-alpha*dot(R.row(i),R.row(i)));
    }
    double psi = S.wavefunc(R,alpha);

    EXPECT_NEAR(psi,g,1e-9);

    mat Rminus = R - h;
    mat Rplus = R + h;
    double gminus = 1;
    double gplus = 1;
    double Vext;
    for(int j = 0; j < numpart; j++){
        gminus *= exp(-alpha*dot(Rplus.row(j),Rplus.row(j)));
        gplus *= exp(-alpha*dot(Rminus.row(j),Rminus.row(j)));
        Vext += 0.5*omega*omega*m*dot(R.row(j),R.row(j));
    }

    double h2 = 1/(h*h);
    double energy =(-0.5*(gplus - 2*g + gminus)*h2/g) + Vext;
    double Esolver = S.energy_num(R,alpha);
    double energy2 = -0.5*(S.wavefunc(Rplus,alpha) - 2*psi + S.wavefunc(Rminus,alpha))/(h*h*psi) + Vext;

    EXPECT_NEAR(Esolver, energy, 1e-4) << "numerical energy not calculated properly";
    EXPECT_NEAR(Esolver, energy2, 1e-4) << "numerical energy not caclulated properly due to variations in wavefunc";
    delete B;
}

TEST(project1test, interaction){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    int mc = 10000; // monte carlo cycles
    int numpart = 500;
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 0.001;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Bruteforce* B = new Bruteforce(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    mat Rij = B->init_pos();
    mat distance = S.distance_part(Rij);

    int g = 0;
    for(int i = 0; i < numpart; i++){
        for(int j = 0; j < numpart; j++){
            g += isnan(distance(i,j));
        }
    }

    int nannumber = numpart*(numpart+1)/2;
    EXPECT_EQ(g,nannumber) << "lower triangular is not nan";

    int o = 0;
    for(int k = 0; k < numpart; k++){
        for(int l = k+1; l < numpart; l++){
            o += isnan(distance(k,l));
        }
    }
    EXPECT_EQ(o,0) << "distance is nan, should be finite";

    mat R1 = ones(numpart,howmanyDs);
    mat distR1 = S.distance_part(R1);

    double shouldbezero = 0;
    for(int r = 0; r < numpart; r++){
        for(int p = r + 1; p < numpart; p++){
            shouldbezero += distR1(r,p);
        }
    }
    EXPECT_EQ(shouldbezero, 0) << "distance between particles should be 0";

    mat R2 = zeros(numpart, howmanyDs);
    R2(numpart-1,howmanyDs-1) = 1;
    R2(numpart - 1,1) = 2;
    R2(numpart - 1, 0) = 1;

    mat R2dist = S.distance_part(R2);
    mat upperR2 = trimatu(R2dist,1);

    double sumR2 = accu(upperR2);
    double shouldbe = (numpart-1)*sqrt(6);

    EXPECT_NEAR(sumR2, shouldbe, 1e-4) << "distance between |Ri-Rj| not calculated correctly";
    delete B;
}
