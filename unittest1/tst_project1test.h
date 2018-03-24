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


TEST(project1test, energy){
    double alpha = 0.5;
    double rho = 0.1;
    double dt = 0.001;
    int mc = 10000; // monte carlo cycles
    int numpart = 10;
    int howmanyDs = 3;
    double beta = 1;
    double hbar = 1;
    double mass = 1;
    double omega = 1;
    double h = 0.001;
    double m = 1;

    Solver S(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt); // initialize Solver class
    //Impsamp* Imp = new Impsamp(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);
    Interact* Int = new Interact(beta, hbar, mass, omega, alpha, rho, mc, numpart, howmanyDs, h, dt);

    ofstream myfile;

    myfile.open("All.dat");     /*CHANGE MY NAME!!!!!!!!!!!!!  DONT YOU DARE NOT CHANGE ME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    //Int->solve_interact(myfile, alpha);
    myfile.close();

    /*
    mat Rtest = ones(numpart);
    double sumlphi = 0;
    mat
    for(int k = 0; k < numpart; k++){

//        sumlphi += -2*alpha*(2 + beta -2*alpha*(2 + beta*beta));
        sumlphi += -4*alpha -2*alpha*beta + 8*alpha*alpha + 4*alpha*beta*beta;
    }
    cout << "sum" << sumlphi << endl;
    Rtest.print();
    mat lphi = Int->lapphi(Rtest,alpha);
    lphi.print();
    */
//    mat R = Int->init_pos_interact();

    mat R = ones(numpart,howmanyDs);
    for(int p = 0; p < numpart; p++){
        R.row(p) = 2*0.0043*p*R.row(p);
    }

    mat distRkj = Int->distance_part(R);
    double EL = 0;
//    double lphi = 0;
    double x;
    double y;
    double z;
    double a = 0.0043;
    mat nphi = -2*alpha*R;
    nphi.col(2) = beta*nphi.col(2);
 //   mat nf;
    double nf2;
    mat lphi = zeros(numpart);
    double sum0 = 0;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;

    for(int k = 0; k < numpart; k++){
        x = R(0);
        y = R(1);
        z = R(2);

        lphi(k) = -2*alpha*(howmanyDs-1 + beta -2*alpha*(x*x + y*y + beta*beta*z*z));
        sum0 += lphi(k);
        for(int j = 0; j < numpart; j++){
            if(k != j){
                mat rkrj = R.row(k) - R.row(j);
                double rkj = distRkj(k,j);
                mat nf = rkrj*a/(rkj*rkj*(rkj-a));
                sum1 += 2*dot(nphi.row(k),nf);
                nf2 = dot(nf,nf);
    //            cout << nf2;
                sum2 += nf2;
                sum3 += -a*a/(rkj*rkj*(rkj-a)*(rkj-a));
            }

        }/*
        for(int j = k+1; j < numpart; j++){
            mat rkrj = R.row(k) - R.row(j);
            double rkj = distRkj(k,j);
            mat nf = rkrj*a/(rkj*rkj*(rkj-a));
            sum1 += 2*dot(nphi.row(k),nf);
            nf2 = dot(nf,nf);
            sum2 += nf2;
            sum3 += -a*a/(rkj*rkj*(rkj-a)*(rkj-a));
        }*/
    }

    mat nftest = Int->nablaf(R,distRkj);
    mat lphitest = Int->lapphi(R,alpha);
    mat nphitest = Int->nablaphi(R,alpha);
//    mat nfnftest = Int->doublesum(R,distRkj);
    mat sum3mat = Int->suma2(distRkj);
    double sum0test = 0;
    double sum1test = 0;
    double sum2test = 0;
    double sum3test = 0;
    for(int l = 0; l < numpart; l++){
        sum0test += lphitest(l);
        sum1test += 2*dot(nphitest.row(l),nftest.row(l));
        sum2test += dot(nftest.row(l),nftest.row(l));
        sum3test += sum3mat(l);
    }
    cout << "sum0 " << sum0 << " | test " << sum0test << endl;
    cout << "sum1 " << sum1 << " | test " << sum1test << endl;
    cout << "sum2 " << sum2 << " | test " << sum2test << endl;
    cout << "sum3 " << sum3 << " |test " << sum3test << endl;
    cout << "ELtest = " << sum0test + sum1test + sum2test + sum3test << endl;

    cout << "EL = " << sum0 + sum1 + sum2 + sum3 << endl;
    delete Int;
}
