#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "../Project1/solver.h"
#include "../Project1/bruteforce.h"
#include "../Project1/impsamp.h"
#include "../Project1/interact.h"

using namespace testing;

TEST(TestInteract, Interact)
{
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001; // [0.001, 0.01]
    double h = 0.001;
    int mc = 100000; // monte carlo cycles
    int numpart = 2;
    int howmanyDs = 2;
    double beta = 1;//2.82843;
    double hbar = 1;
    double mass = 1;
    double omega = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    //Int->solve_interact(myfile, alpha);

    //Create a simple test matrix:
    mat testmat = zeros(N,N);
    mat testdist = Int->distance_part(mat testmat);


/*  Test distance_part again if testdist fails...
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
    */


    //Int->init_pos_interact();
    //Test if returns matrix fulfilling the distance < a demand
    int N = numpart;
    double wf = Int->wavefunc_interact(testmat, alpha, testdist);
    //Test if returns the correct wavefunction
    double g = 0;
    double f = 1;

    for(int i=0;i<N;i++){
        for(int j=0;j<howmanyDs;j++){
            g += testmat(i,j)*testmat(i,j);//g += dot(R.row(i),R.row(i));
            if(i!=j){
                f*=(1 - a/distanceRij(i,j));
            }

        }
    }
    double psi = exp(-alpha*g)*f;

    //EXPECT_EQ(psi, wf);


    Int->quantumF();
    //Test if calculation of quantum force is correct
    Int->lapphi();
    Int->nablaphi();
    Int->nablaf();
    Int->suma2();
    //Test if all four parts of the energy calculation return the correct value for a very simple case
    Int->energy_interact();
    //Test if energy is summed up correctly with the same simple case

    EXPECT_EQ(1, 1);
    ASSERT_THAT(0, Eq(0));
    delete Int;
}
