#include <bits/stdc++.h>
#include "src/define.h"
#include "src/structural_params.h"
#include "src/models.h"
#include "src/estimator.h"

void test(string str="000"){
    cout << str << endl;
}

class instance{
    public:
        double s0;
        double r;
        double sigma;
        double T;
        double lambda;
    instance(double s0_, double r_, double sigma_, double T_, double lambda_):
        s0(s0_), r(r_), sigma(sigma_), T(T_), lambda(lambda_){
        }
    double payoff_array(vd & results, int M=1){
        double x_min = DBL_MAX;
        for (auto i=0; i< results.size(); i+=M){
            if (results[i] < x_min) x_min = results[i];
        }
        return exp(-r * T) * max(results[results.size()-1] - lambda * x_min, 0.);
    }
};

int main(){
    /* We test on a Lookback option where s0=100, r=0.15, sigma=0.1, T=1 and lambda=1.1.
     * For params, we have
     * alpha = 0.5
     * beta = 1
     * V1 = 3.58
     * var(Y0) = 41
     * Real value of price: I0=29.4987
     */
    double s0=100, r=0.15, sigma=0.1, T=1., lambda=1.1;
    double alpha=0.5, beta=1., V1=3.58, varY0=41.;
    double real_value = 8.89343;

    instance eval(s0, r, sigma, T, lambda);

    /*structural_params
     */
    structural_params sp;

    /* multilevel_params:
     * alpha, beta, V1, varY0, simulation_type
     * simulation_type = diffusion / nested
     */
    test("Init multilevel parameters");
    multilevel_params mlp(alpha, beta, V1, varY0, diffusion);

    /* euler_scheme
     * step_size, X0, b, sigma, T
     */
    test("Init euler scheme");
    euler_scheme model(0.01, s0, r, sigma, T);

    /* estimator
     * method_type, structural_params, multilevel_params;
     * method_type = Multilevel_RR, Multilevel_MC
     */
    test("Init estimator");
    estimator est(Multilevel_RR, sp, mlp);

    vvd results = model.simulations(5);
    for (auto i: results){
        for (auto j: i){
            cout << j << " ";
        }
        cout << endl;
    }

    return 0;
}