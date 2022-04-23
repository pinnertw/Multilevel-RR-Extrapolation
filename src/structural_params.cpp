// Define structural_parmas (q, h, n, R, N)
// and multilevel_params (alpha, beta, V1, varY0, theta, hmax, c_tilde, c1, weights)
//
#include "define.h"
#include "structural_params.h"
using namespace std;

// STRUCTURAL PARAMETERS
structural_params::structural_params(method_type method_): method(method_){
}

void structural_params::init_n(int M){
    n = vi(R);
    n[0] = 1;
    for (auto i=1; i<R; i++){
        n[i] = n[i-1] * M;
    }
}

// MULTILEVEL PARAMETERS
multilevel_params::multilevel_params(double alpha_, double beta_, double V1_, double varY0_, simulation_type sim_type_): alpha(alpha_), beta(beta_), V1(V1_), varY0(varY0_), sim_type(sim_type_), hmax(1.){
    theta = sqrt(V1 / varY0);
};

double multilevel_params::dp_exponant(int i, int n_i){
    if (dp_exist[i]) return dp_exponant_array[i];
    else{
        dp_exponant_array[i] =  pow(n_i, alpha);
        dp_exist[i] = true;
        return dp_exponant_array[i];
    }
}

void multilevel_params::init_weights(int R, vi& n){
    vd weights = vd(R, 0.);
    vb dp_exist = vb(R, false);
    vd dp_exponant_array = vd(R, 0.);
    double w_total = 0.;
    for (auto i=1; i<R; i++){
        if ((R - (i + 1)) % 2 != 0){
            weights[i] = -1;
        }
        double numerator = pow(n[i], alpha * (R - 1));
        double denominator = 1;
        for (auto j=0; j < i; j++){
            denominator *= dp_exponant(i, n[i]) - dp_exponant(j, n[j]);
        }
        for (auto j=i+1; j<R; j++){
            denominator *= dp_exponant(j, n[j]) - dp_exponant(i, n[i]);
        }
        weights[i] *= numerator / denominator;
        w_total += weights[i];
    }
}
