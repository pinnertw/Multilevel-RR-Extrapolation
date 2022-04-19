#include <bits/stdc++.h>
#include "params.h"
using namespace std;

#define vvi vector<vector<int>>
#define vi vector<int>
#define vvd vector<vector<double>>
#define vd vector<double>
#define vb vector<bool>
#define d double

class estimator
{
    private:
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
        estimator(int type_, double x0_, double r_, double sigma_, double K_, double T_, 
                double extra_param_ = 1.1, double extra_param2_ = 2.2){
            x0 = x0_;
            r = r_;
            sigma = sigma_;
            K = K_;
            T = T_;
            /* type_:
             * 1 -> BS vanilla call
             * 2 -> Lookback option
             * 3 -> Barrier option
             * 4 -> Compound option
             */
            if (type_ == 1){
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
            else if (type_ == 2){
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
            else if (type_ == 3){
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
            else if (type_ == 4){
                alpha = 1.;
                beta = 1.;
                T2 = extra_param_;
                K2 = extra_param2_;
                payoff_func = [&](vd& results){
                    // TODO
                    return 0.;
                };
            }
        }
};

int main(int argc, char ** argv)
{
    params param(5, 6, 0.5, 1.);
    param.init_weights();
    function<double(double&)> b = [](double a) {return a * 5.;};
    function<double(double&)> sigma = [](double a) {return 5.;};
    euler_schema schema = euler_schema(0.01, 10, 80, b, sigma);
    schema.simulations(10);
    estimator(1, 100, 0.06, 0.4, 80, 1);
    estimator(1, 100, 0.06, 0.4, 80, 2);
    estimator(1, 100, 0.06, 0.4, 80, 3);
    estimator(1, 100, 0.06, 0.4, 80, 4);
    nested_monte_carlo(100., 0.04, 0.3, 10., 20., 1/6, 1/2, 5, 10).simulations(50);
    test();
    return 0;
}
