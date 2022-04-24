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
        double T1;
        double T2;
        double K1;
        double K2;
    instance(double s0_, double r_, double sigma_, double T1_, double T2_, double K1_, double K2_):
        s0(s0_), r(r_), sigma(sigma_), T1(T1_), T2(T2_), K1(K1_), K2(K2_){
        }

    // multilevel
    double payoff_array(vd & results, int n){
        double sum_ = 0.;
        double S_T1 = results[0];
        for (auto i=1; i<n+1; i++){
            if (results[i] > K2) sum_ += results[i];
        }
        sum_ /= (double) n;
        if (K1 > sum_) return K1 - sum_;
        else return 0.;
    }
};

int main(){
    /* We test on a compound option where s0=100, r=0.03, sigma=0.3, T1=1/12., T2=1/2. and K1=6.5, K2=100.
     * For params, we have
     * alpha = beta = 1
     * V1 = 7.2
     * var(Y0) = 9.09
     * Real value of price: I0=1.36857
     */
    double s0=100, r=0.03, sigma=0.3, T1=1/12., T2=1/2., K1=6.5, K2=100;
    double alpha=1., beta=1., V1=7.2, varY0=9.09;
    double real_value = 1.36857;

    instance eval(s0, r, sigma, T1, T2, K1, K2);

    /*structural_params
     */
    structural_params sp;

    /* multilevel_params:
     * alpha, beta, V1, varY0, simulation_type
     * simulation_type = diffusion / nested
     */
    test("Init multilevel parameters");
    multilevel_params mlp(alpha, beta, V1, varY0, nested);

    /* nested_monte_carlo
     * X0, r, sigma, K1, K2, T1, T2
     */
    test("Init nested monte carlo model");
    nested_monte_carlo model(s0, r, sigma, K1, K2, T1, T2);

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
