#include <bits/stdc++.h>
#include "src/structural_params.h"
#include "src/models.h"

void test(string str="000"){
    cout << str << endl;
}

int main(){
    /*structural_params
     * method_type
     * method_type = Multistep_RR, Multilevel_MC
     */
    test("Init structural parameters");
    structural_params sp(Multistep_RR);

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

    return 0;
}
