#include "define.h"
#include "models.h"
using namespace std;

// EULER_SCHEME
euler_scheme::euler_scheme(int total_step_, double X0_, double b_, double sigma_, double T_):
    total_step(total_step_), X0(X0_), b(b_), sigma(sigma_), T(T_){
        step_size = (T / (double) total_step);
        sqrt_step_size = sqrt(step_size);
}
void euler_scheme::reset_step(int total_step_){
    total_step = total_step_;
    step_size = (T / (double) total_step);
    sqrt_step_size = sqrt(step_size);
}
void euler_scheme::get_normal_distribution(int N){
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed);
    random_normal = vd(N * total_step, 0.);
    normal_distribution<double> dist(0., 1.);
    for (auto i=0; i<N*total_step; i++){
        random_normal[i] = dist(generator);
    }
}
vvvd euler_scheme::simulations(int N, int M){
    // Return results[N][0/1][total_step]
    get_normal_distribution(N);
    vvvd results = vvvd(N, vvd(2, vd(total_step + 1, 0.)));
    double sum_random = 0.;
    for (auto i=0; i <N; i++){
        results[i][0][0] = X0;
        results[i][1][0] = X0;
        for (auto j=1; j<=total_step; j++){
            // Simulation with little steps
            results[i][0][j] = results[i][0][j-1] + step_size * b * results[i][0][j-1] + 
                sqrt_step_size * sigma * results[i][0][j-1] * random_normal[i * total_step + j-1];
            sum_random += random_normal[i * total_step + j - 1];
            // Simulation with greater steps
            if (j % M == 0){
                results[i][1][j] = results[i][1][j-M] + M * step_size * b * results[i][1][j-M] +
                    sqrt_step_size * sigma * results[i][1][j-M] * sum_random;
                sum_random = 0.;
            }
        }
    }
    return results;
}

// NESTED_MONTE_CARLO
nested_monte_carlo::nested_monte_carlo(double X0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_)
: X0(X0_), r(r_), sigma(sigma_), K1(K1_), K2(K2_), T1(T1_), T2(T2_){
    n1 = 1;
    n2 = 3;
}

void nested_monte_carlo::set_n(int n1_, int n2_){
    n1 = n1_;
    n2 = n2_;
}

void nested_monte_carlo::get_normal_distribution(int N){
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed);
    random_normal = vd(N * (n2+1), 0.);
    normal_distribution<double> dist(0., 1.);
    for (auto i=0; i < N*(n2+1); i++){
        random_normal[i] = dist(generator);
    }
}
vvd nested_monte_carlo::simulations(int N){
    get_normal_distribution(N);
    vvd results = vvd(N, vd((n2+1), 0.));
    double shift1 = exp((r-0.5*sigma*sigma) * T1);
    double shift2 = exp((r-0.5*sigma*sigma) * (T2 - T1));
    for (auto i=0; i<N; i++){
        // results[i][0] = S_T1, results[i][1 ~ n2] = S_T2
        results[i][0] = X0 * shift1 * exp(sigma*sqrt(T1) * random_normal[i*(n2+1)]);
        for (auto j=1; j<n2 + 1; j++){
            results[i][j] = results[i][0] * shift2 * exp(sigma * sqrt(T2 - T1) * random_normal[i * (n2 + 1) + j]);
        }
    }
    return results;
}
