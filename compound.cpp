#include <bits/stdc++.h>
#include "src/define.h"
#include "src/structural_params.h"
#include "src/models.h"
#include "src/estimator.h"

void test(string str="000"){
    cout << str << endl;
}

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
