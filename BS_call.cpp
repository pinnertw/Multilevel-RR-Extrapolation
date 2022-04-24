#include <bits/stdc++.h>
#include <sys/time.h>
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
        double K;
    instance(double s0_, double r_, double sigma_, double T_, double K_):
        s0(s0_), r(r_), sigma(sigma_), T(T_), K(K_){
        }
    double payoff(double value){
        if (value > K) return exp(-r * T) * (value - K);
        else return 0.;
    }

    double payoff_array(vd & results){
        return payoff(results[results.size()-1]);
    }

    double payoff_array(vvd & results, bool R_is_1=false){
        if (R_is_1) return payoff_array(results[0]);
        return payoff_array(results[0]) - payoff_array(results[1]);
    }

    double payoff_array(vvvd &results, bool R_is_1=false){
        double sum=0.;
        for (vvd result_loc: results){
            sum += payoff_array(result_loc, R_is_1);
        }
        sum /= (double) results.size();
        return sum;
    }

};

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

    instance eval(s0, r, sigma, T, K);

    /*structural_params
     * q
     * h, n, R, N
     * T
     */
    structural_params sp;

    /* multilevel_params:
     * alpha, beta, V1, varY0, simulation_type
     * simulation_type = diffusion / nested
     */
    multilevel_params mlp(alpha, beta, V1, varY0, diffusion);

    /* euler_scheme
     * total_step, X0, b, sigma, T
     */
    euler_scheme model(10, s0, r, sigma, T);

    /* estimator
     * method_type, structural_params, multilevel_params;
     * method_type = Multilevel_RR, Multilevel_MC
     */
    estimator est(Multilevel_RR, sp, mlp);

    // Autotune
    for (int k=1; k<3; k++){
        struct timeval t1, t2;
        double duration1, duration2;
        gettimeofday(&t1, NULL);

        // Autotuning
        double epsilon = pow(2, -k);
        est.auto_tune(epsilon);
        est.sp.init_n(est.M);

        gettimeofday(&t2, NULL);
        duration1 = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec)/1e6);
        
        // Get total number for each simulation. Then Run!
        int base_num = T * est.h_inverse;
        vd payoff_results(256);
        for (int i=0; i<256; i++){
            // level 1
            model.reset_step(base_num);

            vvvd simulations = model.simulations(ceil(est.sp.N * est.sp.q[0]), 1);

            double result = eval.payoff_array(simulations, true);
            
            // level 2 - R
            for (int j=1; j<est.sp.R; j++){
                model.reset_step(base_num * est.sp.n[j]);
                vvvd simulations = model.simulations(ceil(est.sp.N * est.sp.q[j]), est.M);
                result += est.sp.T[j][1] * eval.payoff_array(simulations);
            }
            payoff_results[i] = result;
        }
        double epsilon_L = RMSE(payoff_results, real_value);
        double bias_ = bias(payoff_results, real_value);
        double var_ = var(payoff_results);
        // End time.
        gettimeofday(&t1, NULL);
        duration2 = (t1.tv_sec - t2.tv_sec) + ((t1.tv_usec - t2.tv_usec)/1e6);
        cout << k << "," << duration1 << "," << duration2 << "," << epsilon_L << "," << bias_ << "," << var_ << "," <<
            est.sp.R << "," << est.M << "," << est.h_inverse << "," << est.sp.N << "," << est.cost() << endl;
    }

    return 0;
}
