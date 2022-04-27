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
        for (auto i=1; i<n+1; i++){
            if (results[i] > K2) sum_ += results[i];
        }
        sum_ /= (double) n;
        if (K1 > sum_) return K1 - sum_;
        else return 0.;
    }
    double payoff_array(vd &results, int n1, int n2){
        return payoff_array(results, n2) - payoff_array(results, n1);
    }
    double payoff_array(vvd & results, int n1){
        double sum_ = 0.;
        for (auto i: results){
            sum_ += payoff_array(i, n1);
        }
        sum_ /= (double) results.size();
        return sum_;
    }
    double payoff_array(vvd &results, int n1, int n2){
        double sum_ = 0.;
        for (auto i: results){
            sum_ += payoff_array(i, n1, n2);
        }
        sum_ /= (double) results.size();
        return sum_;
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
    multilevel_params mlp(alpha, beta, V1, varY0, nested);

    /* nested_monte_carlo
     * X0, r, sigma, K1, K2, T1, T2
     */
    nested_monte_carlo model(s0, r, sigma, K1, K2, T1, T2);

    /* estimator
     * method_type, structural_params, multilevel_params;
     * method_type = Multilevel_RR, Multilevel_MC
     */
#if MLMC
    estimator est(Multilevel_MC, sp, mlp);
#else
    estimator est(Multilevel_RR, sp, mlp);
#endif

    cout << "k,t1,t2,epsilon_L,bias,variance,R,M,h_inverse,N,cost" << endl;
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
        vd payoff_results(256);
        for (int i=0; i<256; i++){
            // level 1
            model.set_n(1, est.sp.n[0]);

            vvd simulations = model.simulations(ceil(est.sp.N * est.sp.q[0]));

            double result = eval.payoff_array(simulations, est.sp.n[0]);
            
            //cout << result << endl;
            // level 2 - R
            for (int j=1; j<est.sp.R; j++){
                model.set_n(est.sp.n[j-1], est.sp.n[j]);
                vvd simulations = model.simulations(ceil(est.sp.N * est.sp.q[j]));
                result += est.sp.T[j][1] * eval.payoff_array(simulations, est.sp.n[j-1], est.sp.n[j]);
            }
            payoff_results[i] = result;
            //cout << result << endl;
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
