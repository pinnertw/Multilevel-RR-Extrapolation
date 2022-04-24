#include <bits/stdc++.h>
#include "src/define.h"
#include "src/structural_params.h"
#include "src/models.h"
#include "src/estimator.h"

void test(string str="000"){
    cout << str << endl;
}

int main(){
    /* We test on a Vanilla Call option where s0=100, r=0.06, sigma=0.4, T=1 and K=80.
     * For params, we have
     * alpha = beta = 1
     * V1 = 56
     * var(Y0) = 876
     * Real value of price: I0=29.4987
     */
    double s0=100, r=0.06, sigma=0.4, T=1., K=80.;
    double alpha=1., beta=1., V1=56., varY0=876.;
    double real_value = 29.4987;

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
