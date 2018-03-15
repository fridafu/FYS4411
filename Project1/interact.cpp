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
{}
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
    double a = 0.043;
    int counter = 0;
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
                        Rtull(o,dude) = (gaussianRNG(genMT64) - 0.5)*rho;
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


