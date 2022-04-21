#include <bits/stdc++.h>
#include "params.h"
using namespace std;

#define vvi vector<vector<int>>
#define vi vector<int>
#define vvd vector<vector<double>>
#define vd vector<double>
#define vb vector<bool>
#define d double

enum model{call, lookback, barrier, compound};
enum method{MLMC, MLRR};

class estimator
{
    private:
        model type;
        double x0;
        double r;
        double sigma;
        double K;
        double T;
        function<double(vd &)> payoff_func;
        double alpha;
        double beta;
        // Other params for different options
        double lambda;
        double B;
        double T2;
        double K2;
        // TODO
        // Real value, Variance, Complexity, Cost, auto-tuning
    public:
        estimator(model type_, double x0_, double r_, double sigma_, double K_, double T_, 
                double extra_param_ = 1.1, double extra_param2_ = 2.2){
            type = type_;
            x0 = x0_;
            r = r_;
            sigma = sigma_;
            K = K_;
            T = T_;
            if (type_ == call){
                alpha = 1.;
                beta = 1.;
                payoff_func = [&](vd& results) {
                    int total_steps = results.size();
                    if (results[total_steps-1] > K){
                        return exp(-r * T) * (results[total_steps-1] - K);
                    }
                    else{
                        return 0.;
                    }
                };
            }
            else if (type_ == lookback){
                alpha = 0.5;
                beta = 1.;
                lambda = extra_param_;
                payoff_func = [&](vd& results){
                    double x_min = 1e9;
                    for (double x : results){
                        if(x < x_min) x_min = x;
                    }
                    return exp(-r * T) * max(results[results.size()-1] - lambda * x_min, 0.);
                };
            }
            else if (type_ == barrier){
                alpha = 0.5;
                beta = 0.5;
                B = extra_param_;
                payoff_func = [&](vd& results){
                    bool passed = false;
                    for (double x : results){
                        if (x > B) {
                            passed = true;
                            break;
                        }
                    };
                    if (!passed) return exp(-r * T) * max(results[results.size()-1] - K, 0.);
                    else return 0.;
                };
            }
            else if (type_ == compound){
                alpha = 1.;
                beta = 1.;
                K2 = extra_param_;
                T2 = extra_param2_;
                payoff_func = [&](vd& results){
                    double sum_ = 0.;
                    for (auto x: results){
                        if (x > K2){
                            sum_ += x - K2;
                        }
                    }
                    sum_ /= (double) results.size();
                    if(K > sum_) return K - sum_;
                    else return 0.;
                };
            }
        }
        void simulations(){
            // TODO step size autotune, simulations with different allocations;
            if(type == call){
                double step_size = 0.01;
                int total_step = T / step_size;
                function<double(double &)> b_func = [&](double t){return r * t;};
                function<double(double &)> sigma_func = [&](double t){return sigma * t;};
                euler_scheme generator = euler_scheme(step_size, total_step, x0, b_func, sigma_func);
            }
            else if (type == lookback){
                double step_size = 0.01;
                int total_step = T / step_size;
                function<double(double &)> b_func = [&](double t){return r * t;};
                function<double(double &)> sigma_func = [&](double t){return sigma * t;};
                euler_scheme generator = euler_scheme(step_size, total_step, x0, b_func, sigma_func);
            }
            else if (type == barrier){
                double step_size = 0.01;
                int total_step = T / step_size;
                function<double(double &)> b_func = [&](double t){return r * t;};
                function<double(double &)> sigma_func = [&](double t){return sigma * t;};
                euler_scheme generator = euler_scheme(step_size, total_step, x0, b_func, sigma_func);
            }
            else if (type == compound){
                int n1 = 1;
                int n2 = 5;
                nested_monte_carlo generator = nested_monte_carlo(x0, r, sigma, K, T, K2, T2, n1, n2);
            }
        }

        // Autotune
        // double 
        params find_best_params(double epsilon){
            int optimal_M = 2;
            double min_complexity = DBL_MAX;
            params param = params(optimal_M, 2, alpha, beta);
            for (auto M=2; M <= 10; M++){
                auto_tune(param, M, epsilon);
                if (complexity(param) < min_complexity){
                    min_complexity = complexity(param);
                    optimal_M = M;
                }
            }
            auto_tune(param, optimal_M, epsilon);
            return param;
        }
        // TODO
        void auto_tune(params param, int M, double epsilon); // change other params
        double complexity(params param); // Return complexity
};

int main(int argc, char ** argv)
{
    params param(5, 6, 0.5, 1.);
    param.init_weights();
    function<double(double&)> b = [](double a) {return a * 5.;};
    function<double(double&)> sigma = [](double a) {return 5.;};
    euler_scheme schema = euler_scheme(0.01, 10, 80, b, sigma);
    schema.simulations(10);
    estimator(call, 100, 0.06, 0.4, 80, 1);
    estimator(lookback, 100, 0.15, 0.1, 0, 1, 1.1);
    estimator(barrier, 100, 0., 0.15, 100, 1, 120);
    estimator(compound, 100, 0.03, 0.4, 6.5, 1./12., 100, 1/2.);
    test();
    return 0;
}
