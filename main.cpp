#include <bits/stdc++.h>
#include "src/define.h"
#include "src/structural_params.h"
#include "src/models.h"
#include "src/estimator.h"

void test(string str="000"){
    cout << str << endl;
}

int main(){
    /*structural_params
     */
    test("Init structural parameters");
    structural_params sp;

    /* multilevel_params:
     * alpha, beta, V1, varY0, simulation_type
     * simulation_type = diffusion / nested
     */
    test("Init multilevel parameters");
    multilevel_params mlp(1., 1., 1., 1., diffusion);

    /* euler_scheme
     * step_size, X0, b, sigma, T
     */
    test("Init euler scheme");
    euler_scheme model(0.01, 80, 0.01, 0.3, 1);

    /* nested_monte_carlo
     * X0, r, sigma, K1, K2, T1, T2
     */
    test("Init nested monte carlo");
    nested_monte_carlo model2(80, 0.01, 0.3, 50, 80, 1/12., 1/2.);

    /* estimator
     * method_type, structural_params, multilevel_params;
     * method_type = Multilevel_RR, Multilevel_MC
     */
    test("Init estimator");
    estimator est(Multilevel_RR, sp, mlp);

    return 0;
}
