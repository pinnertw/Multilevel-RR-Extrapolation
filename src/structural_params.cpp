#include "define.h"
using namespace std;

double v(int M, double alpha, double beta, double epsilon){
    if (beta > 1) return pow(epsilon, -2);
    else if (beta == 1) return pow(epsilon, -2) * log(1 / epsilon);
    else return pow(epsilon, -2) * exp((1 - beta) / sqrt(alpha) * sqrt(2 * log(1 / epsilon) * log(M)));// o(epsilon **-eta) for eta > 0
}

enum method_type {Multistep_RR, Multilevel_MC, Multilevel_RR};

class structural_params{
    public:
        // q
        vd q;
        // pi0
        double h;
        vi n;
        int R;
        int N;

        // Method type
        method_type method; 

        // Usage
        double fact_n;
        vvd T;

    void init_n(int M){
        n = vi(R);
        n[0] = 1;
        for (auto i=1; i<R; i++){
            n[i] = n[i-1] * M;
        }
    }
};


enum simulation_type {diffusion, nested};

class multilevel_params{
    public:
        // Given: alpha/beta
        double alpha; // weak error rate, exponant in polynome.
        double beta; // Strong approximation error, L2 error <= V1 h^beta

        // Given or to calculate: V1, var(Y0), theta), TODO
        double V1;
        double varY0;
        double theta = sqrt(V1 / varY0);
        simulation_type sim_type;

        double hmax;
        double c_tilde=1;
        double c1=1;

        // weights
        vd weights;
        double w_total;

        vb dp_exist;
        vd dp_exponant_array;

    double dp_exponant(int i, int n_i){
        if (dp_exist[i]) return dp_exponant_array[i];
        else{
            dp_exponant_array[i] =  pow(n_i, alpha);
            dp_exist[i] = true;
            return dp_exponant_array[i];
        }
    }

    void init_weights(int R, vi& n){
        vd weights = vd(R, 0.);
        vb dp_exist = vb(R, false);
        vd dp_exponant_array = vd(R, 0.);
        //double w_tilde = (R%2 == 1)? 1: -1;
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
            //w_tilde /= dp_exponant(i, n[i]);
        }
    }
};

