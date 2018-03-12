#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "../Project1/solver.h"

using namespace testing;


TEST(project1test, wavefunction)
{
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    int mc = 10000; // monte carlo cycles
    int numpart = 500; //CHANGE THE NAME!!!!!!!!!!!!!!!!!!!!!!!!!
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 0.001;
    double m = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    //mat R = ones(numpart,howmanyDs)*rho;//S.init_pos();
    mat R = S.init_pos();
    //cout << R << endl;
    double g = 1;
    for(int i = 0; i < numpart; i++){
        g *= exp(-alpha*dot(R.row(i),R.row(i)));
    }
    double psi = S.wavefunc(R,alpha);
    cout << g << endl;
    cout << psi << endl;
    //EXPECT_EQ(1, 1);
    //ASSERT_THAT(0, Eq(0));
    EXPECT_NEAR(psi,g,1e-9);

    //ASSERT_DOUBLE_EQ(psi,g);

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
    cout << "g+ " << gplus << endl;
    cout << "g- " << gminus << endl;
    cout << "g " << g << endl;
    double h2 = 1/(h*h);
    double energy =(-0.5*(gplus - 2*g + gminus)*h2/g) + Vext;
    double Esolver = S.energy_num(R,alpha);
    double energy2 = -0.5*(S.wavefunc(Rplus,alpha) - 2*psi + S.wavefunc(Rminus,alpha))/(h*h*psi);
    cout << energy << endl;
    cout << Esolver << endl;
    cout << energy2 << endl;

    EXPECT_NEAR(Esolver, energy,0.5);
    EXPECT_NEAR(Esolver, energy2, 0.5);
}
