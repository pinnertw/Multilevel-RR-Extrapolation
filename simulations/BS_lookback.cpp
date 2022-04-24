#include <bits/stdc++.h>
#include <sys/time.h>
#include "../src/define.h"
#include "../src/structural_params.h"
#include "../src/models.h"
#include "../src/estimator.h"

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
    double payoff_array(vd & results, int M){
        double x_min = DBL_MAX;
        for (auto i=0; i< results.size(); i+=M){
            if (results[i] < x_min) x_min = results[i];
        }
        return exp(-r * T) * max(results[results.size()-1] - lambda * x_min, 0.);
    }
    double payoff_array(vvd & results, int M, bool R_is_1=false){
        if (R_is_1) return payoff_array(results[0], 1);
        return payoff_array(results[0], 1) - payoff_array(results[1], M);
    }
    double payoff_array(vvvd & results, int M, bool R_is_1=false){
        double sum=0.;
        for (vvd result_loc: results){
            sum += payoff_array(result_loc, M, R_is_1);
        }
        sum /= (double) results.size();
        return sum;
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
    multilevel_params mlp(alpha, beta, V1, varY0, diffusion);

    /* euler_scheme
     * step_size, X0, b, sigma, T
     */
    euler_scheme model(10, s0, r, sigma, T);

    /* estimator
     * method_type, structural_params, multilevel_params;
     * method_type = Multilevel_RR, Multilevel_MC
     */
    estimator est(Multilevel_RR, sp, mlp);

    for (int k=1; k<10; k++){
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

            double result = eval.payoff_array(simulations, est.M, true);
            
            //cout << result << endl;
            // level 2 - R
            for (int j=1; j<est.sp.R; j++){
                model.reset_step(base_num * est.sp.n[j]);
                vvvd simulations = model.simulations(ceil(est.sp.N * est.sp.q[j]), est.M);
                result += est.sp.T[j][1] * eval.payoff_array(simulations, est.M);
                //cout << est.sp.T[j][1] * eval.payoff_array(simulations) << endl;
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
