#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include "../Project1/solver.h"
#include "../Project1/bruteforce.h"
#include "../Project1/interact.h"

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
//    EXPECT_EQ(g,nannumber) << "lower triangular is not nan";

    int o = 0;
    for(int k = 0; k < numpart; k++){
        for(int l = k+1; l < numpart; l++){
            o += isnan(distance(k,l));
        }
    }
//    EXPECT_EQ(o,0) << "distance is nan, should be finite";

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

TEST(project1test,energyinteraction){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    int mc = 10000; // monte carlo cycles
    int numpart = 2;
    int howmanyDs = 2;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 0.001;
    double a = 0.0043;
    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    // test distance function for two dimensions
    ofstream myfile;
    mat pos = S.init_pos_gaus();
    mat Rtest = pos;
//    Rtest.row(0) = Rtest.row(0) + 1;

    mat Rtdist = S.distance_part(Rtest);

    EXPECT_NEAR(sqrt(2),Rtdist(0,1),1e-4) << "wrong in dist_part(0,1)";
    EXPECT_NEAR(sqrt(2),Rtdist(1,0),1e-4) << "wrong in dist_part(1,0)";
    mat lapphi = Int->lapphi(Rtest,alpha);
    mat lapphitest = zeros(howmanyDs);
    for(int i = 0; i < numpart;i++){
        lapphitest(i) = -2*alpha*((howmanyDs-1) - 2*alpha*dot(Rtest.row(i),Rtest.row(i)));
        cout << i << endl;
        EXPECT_NEAR(lapphitest(i),lapphi(i), 1e-5) << "something wring with lapphi";
    }
    mat nablaphi = Int->nablaphi(Rtest,alpha);
    mat nablaphitest = zeros(numpart,howmanyDs);
    double temp = 0;
    for(int j = 0; j < numpart; j++){
//        nablaphitest(j) = -4*alpha*Rtest.row(j);
        for(int p = 0; p < howmanyDs; p++){
            temp = -2*alpha*Rtest(j,p);
            cout << "(" << j <<", " << p << ")" << endl;
            EXPECT_NEAR(temp,nablaphi(j,p),1e-5) << "nablaphi not working!";
        }
    }

    mat nablaf = Int->nablaf(Rtest,Rtdist);
    mat nablaftest = zeros(numpart,howmanyDs);
    double rkj = 0;
//    mat sum = zeros(N,dim);

    mat rk;
    for(int k = 0; k < numpart; k++){
        rk = Rtest.row(k);
        for(int j = 0; j < numpart; j++){
            if(k != j){
                rkj = Rtdist(k,j);
                nablaftest.row(k) += (rk-Rtest.row(j))*a/(rkj*rkj*(rkj - a));
            }
        }
    }
    cout << "testing nablaf" << endl;
    for(int q = 0; q < numpart; q++){
        for(int l = 0; l < howmanyDs; l++){
            EXPECT_NEAR(nablaftest(q,l), nablaf(q,l), 1e-5) << "nablaf not right";
        }
    }
    cout << "testing nablaf^2" << endl;
    mat nablaf2 = Int->nablafsquared(Rtdist,Rtest);
    mat nablaf2test = zeros(numpart);
    for(int u = 0; u < numpart; u++){
        for(int y = 0; y < numpart; y++){
            if(u != y){
                nablaf2test(u) += dot(nablaf.row(u),nablaf.row(y));
            }
        }
        EXPECT_NEAR(nablaf2test(u), nablaf2(u),1e-5) << "nablaf^2 not right";

    }
    cout << "testing suma2" << endl;
    mat suma2 = Int->suma2(Rtdist, Rtest);
    mat suma2test = zeros(numpart);
    for(int t = 0; t < numpart; t++){
        for(int w = 0; w < numpart; w++){
            if(t != w){
                suma2test(t) += ((howmanyDs-1)*a*Rtdist(t,w) - (howmanyDs - 1)*a*a + a*a - 2*a*Rtdist(t,w))/(Rtdist(t,w)*Rtdist(t,w)*(Rtdist(t,w)-a)*(Rtdist(t,w)-a));
            }
        }
        EXPECT_NEAR(suma2test(t), suma2(t),1e-5) << "suma2 not right!";
    }
    double energy = Int->energy_interact(Rtest, alpha);
    double energytest = 0;
    double Vext = 0;
    double r2;
    for(int j = 0; j < numpart; j++) {
        r2 = 0;
        for(int q = 0; q < howmanyDs; q++) {
            r2 += Rtest(j,q)*Rtest(j,q);
        }
        Vext += 0.5*r2; //calculate potential energy
    }
    for(int h = 0; h < numpart; h++){
        energytest += lapphitest(h) + 2*dot(nablaphitest.row(h),nablaftest.row(h)) + nablaf2test(h) + suma2test(h);
    }
    cout << "testing energy" << endl;
    EXPECT_NEAR(energy,-0.5*energytest + Vext,1e-2) << "energy not right";
//    double twoa = 1.4142 - a;
//    double testenergy = 4*alpha- +4*alpha*a/twoa + a*a/(twoa*twoa) - 2*a/(sqrt(2)*twoa*twoa);
//    EXPECT_NEAR(testenergy,energy,0.1) << "energy fo 2x2 not right!";

//    myfile.open("testtest.dat");     /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//    Int->solve_interact(myfile, alpha);
//    myfile.close();
    delete Int;

}
