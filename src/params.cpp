#include <bits/stdc++.h>
#include "params.h"
using namespace std;


#define vvi vector<vector<int>>
#define vi vector<int>
#define vvd vector<vector<double>>
#define vd vector<double>
#define vb vector<bool>
#define d double

void test(){
    cout << 0 << endl;
}

params::params(int M_, int R_, double alpha_, double beta_){
        // Init params with M as root, and level R
        // Then get the refiners with n_j = M^j;
        M = M_;
        R = R_;
        alpha = alpha_;
        beta = beta_;
        dp_exist = vb(R, false);
        dp_exponant_array = vd(R, 0);
        // INIT
        init_refiners();
        init_weights();
        print_params();
    }

void params::print_params(){
#if PRINT1
    cout << "M, R, alpha, beta"<<endl;
    cout << M << ", " << R << ", " << alpha << ", " << beta << endl;
#endif
}

double params::dp_exponant(int i){
    if (dp_exist[i]) return dp_exponant_array[i];
    else{
        dp_exponant_array[i] = pow(refiners[i], alpha);
        dp_exist[i] = true;
        return dp_exponant_array[i];
    }
}

void params::init_refiners(){
    refiners = vi(R);
    refiners[0] = 1;
    for (auto i = 1; i < R; i++){
        refiners[i] = refiners[i-1] * M;
    }
#if PRINT1
    cout << "Print refiners coefficients" << endl;
    for (auto i = 0; i < R; i++){
        cout << refiners[i] << endl;
    }
#endif
}

void params::init_weights(){
    weights = vd(R, 1.);
    for (auto i = 0; i < R; i++){
        if ( (R - (i + 1)) % 2 != 0){
            weights[i] = -1;
            w_tilde = 1;
        }
        else{
            w_tilde = -1;
        }
        double numerator = pow(refiners[i], alpha * (R - 1));
        double denominator = 1;
        for (auto j = 0; j < i; j++){
            denominator *= dp_exponant(i) - dp_exponant(j);
        }
        for (auto j = i + 1; j < R; j++){
            denominator *= dp_exponant(j) - dp_exponant(i);
        }
        weights[i] *= numerator / denominator;
        w_tilde /= dp_exponant(i);
    }
#if PRINT1
    cout << "Print weight coefficients" << endl;
    for (auto i = 0; i < R; i++){
        cout << weights[i] << ' ';
    }
    cout << endl;
    cout << "w_tilde = " << w_tilde << endl;
#endif
}

// EULER SCHEME
euler_scheme::euler_scheme(double step_size_, int total_step_, double X0_, function<double(double &)> b_, function<double(double &)> sigma_){
    step_size = step_size_;
    sqrt_step_size = sqrt(step_size);
    total_step = total_step_;
    X0 = X0_;
    b = b_;
    sigma = sigma_;
}

void euler_scheme::reset_step(double step_size_, int total_step_){
    step_size = step_size_;
    sqrt_step_size = sqrt(step_size);
    total_step = total_step_;
}

void euler_scheme::get_normal_distribution(int N){
    mt19937_64 generator;
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    random_normal = vd(N*total_step, 0.);
    normal_distribution<double> dist(0., 1.);
    for (auto i=0; i < N*total_step; i++){
        random_normal[i] = dist(generator);
    }
}

vvd euler_scheme::simulations(int N){
    get_normal_distribution(N);
    vvd results = vvd(N, vd(total_step + 1, 0.));
    for (auto i=0; i < N; i++){
        results[i][0] = X0;
        for (auto j = 1; j <= total_step; j++){
            results[i][j] = results[i][j-1] + step_size * b(results[i][j-1]) + sqrt_step_size * sigma(results[i][j-1]) * random_normal[i * total_step + j];
        }
    }
#if PRINT
    for (auto i=0; i < N; i++){
        for (auto j=0; j < N; j++){
            cout << results[i][j] << ' ';
        }
        cout << endl;
    }
#endif
    return results;
}

// NESTED MONTE CARLO
nested_monte_carlo::nested_monte_carlo(double s0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_, int n1_, int n2_){
    s0 = s0_;
    r = r_;
    sigma = sigma_;
    K1 = K1_;
    K2 = K2_;
    T1 = T1_;
    T2 = T2_;
    n1 = n1_;
    n2 = n2_;
}

void nested_monte_carlo::reset_n(int n1_, int n2_){
    n1 = n1_;
    n2 = n2_;
}

void nested_monte_carlo::get_normal_distribution(int N){
    mt19937_64 generator;
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    random_normal = vd(N*(n2 + 1), 0.);
    normal_distribution<double> dist(0., 1.);
    for (auto i=0; i < N*(n2 + 1); i++){
        random_normal[i] = dist(generator);
    }
}

void nested_monte_carlo::simulations(int N){
    get_normal_distribution(N);
    vvd results = vvd(N, vd(n2, 0.));
    for (auto i=0; i < N; i++){
        double S_T1 = s0 * exp((r - 0.5 * sigma * sigma) * T1 + sigma * sqrt(T1) * random_normal[i * (n2 + 1)]);
        for (auto j=0; j < n2; j++){
            results[i][j] = S_T1 * exp((r - 0.5 * sigma * sigma) * (T2 - T1) + sigma * sqrt(T2 - T1) * random_normal[i * (n2 + 1) + j]);
        }
    }
}
