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

// EULER_SCHEME
euler_scheme_MSRR::euler_scheme_MSRR(int total_step_, double X0_, double b_, double sigma_, double T_):
    total_step(total_step_), X0(X0_), b(b_), sigma(sigma_), T(T_){
        step_size = (T / (double) total_step);
        sqrt_step_size = sqrt(step_size);
        set_alpha_r(total_step);
}
void euler_scheme_MSRR::reset_step(int total_step_){
    total_step = total_step_;
    step_size = (T / (double) total_step);
    sqrt_step_size = sqrt(step_size);
    set_alpha_r(total_step);
}

void euler_scheme_MSRR::set_alpha_r(int R){
    if (R == 2){
        alpha_rs = {-1, 2};
    }
    else if (R == 3){
        alpha_rs = {0.5, -4, 4.5};
    }
    else if (R == 4){
        alpha_rs = {-1/6., 4, -13.5, 32/3.};
    }
}

void euler_scheme_MSRR::get_normal_distribution(int N){
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed);
    random_normal = vvd(N, vd((total_step-1)*2, 0.));
    normal_distribution<double> dist(0., 1.);
    for (auto i=0; i<N; i++){
        for (auto j=0; j < (total_step - 1) * 2; j++){
            random_normal[i][j] = dist(generator);
        }
    }
}
vvvd euler_scheme_MSRR::simulations(int N){
    // Return results[N][0/1][total_step]
    get_normal_distribution(N);
    vvvd results = vvvd(N, vvd(total_step, vd(total_step + 1, 0.)));
    double sum_random = 0.;
    for (auto i=0; i <N; i++){
        if (total_step == 2){
            results[i][0][0] = X0;
            results[i][1][0] = X0;
            // U^2
            results[i][1][1] = results[i][1][0] * (1 + step_size * b + sqrt_step_size * sigma * random_normal[i][0]);
            results[i][1][2] = results[i][1][1] * (1 + step_size * b + sqrt_step_size * sigma * random_normal[i][1]);
            // U^1
            results[i][0][1] = results[i][0][0] * (1 + 2 * step_size * b + sqrt_step_size * sigma * (random_normal[i][0] + random_normal[i][1]));
        }
        else if (total_step == 3){
            results[i][0][0] = X0;
            results[i][1][0] = X0;
            results[i][2][0] = X0;
            // U^3
            results[i][2][1] = results[i][2][0] * (1 + step_size * b + sqrt_step_size * sigma * random_normal[i][0]);
            results[i][2][2] = results[i][2][1] * (1 + step_size * b + sqrt_step_size * sigma * (random_normal[i][1] + random_normal[i][2]) / sqrt(2));
            results[i][2][3] = results[i][2][2] * (1 + step_size * b + sqrt_step_size * sigma * random_normal[i][3]);
            // U^2
            double u2_1 = (sqrt(2) * random_normal[i][0] + random_normal[i][1]) / sqrt(3);
            double u2_2 = (random_normal[i][2] + sqrt(2) * random_normal[i][3]) / sqrt(3);
            results[i][1][1] = results[i][1][0] * (1 + 1.5 * step_size * b + sqrt(1.5) * sqrt_step_size * sigma * u2_1);
            results[i][1][2] = results[i][1][1] * (1 + 1.5 * step_size * b + sqrt(1.5) * sqrt_step_size * sigma * u2_2);
            // U^1
            results[i][0][1] = results[i][0][0] * (1 + 3 * step_size * b + sqrt(3)* sqrt_step_size * sigma * (u2_1 + u2_2));
        }
        else if (total_step == 4){
            results[i][0][0] = X0;
            results[i][1][0] = X0;
            results[i][2][0] = X0;
            results[i][3][0] = X0;
            // U^4
            double u4_1 = random_normal[i][0];
            double u4_2 = (random_normal[i][1] + sqrt(2) * random_normal[i][2])/sqrt(3);
            double u4_3 = (random_normal[i][3] * sqrt(2) + random_normal[i][4])/sqrt(3);
            double u4_4 = random_normal[i][5];
            results[i][3][1] = results[i][3][0] * (1 + step_size * b + sqrt_step_size * sigma * u4_1);
            results[i][3][2] = results[i][3][1] * (1 + step_size * b + sqrt_step_size * sigma * u4_2);
            results[i][3][3] = results[i][3][2] * (1 + step_size * b + sqrt_step_size * sigma * u4_3);
            results[i][3][4] = results[i][3][3] * (1 + step_size * b + sqrt_step_size * sigma * u4_4);
            // U^3
            double u3_1 = (random_normal[i][0] * sqrt(3) + random_normal[i][1])/2.;
            double u3_2 = (random_normal[i][2] + random_normal[i][3]) / sqrt(2);
            double u3_3 = (random_normal[i][4] + sqrt(3) * random_normal[i][5])/2.;
            results[i][2][1] = results[i][2][0] * (1 + 4/3. * step_size * b + sqrt(4/3.) * sqrt_step_size * sigma * (random_normal[i][0] * sqrt(3) + random_normal[i][1])/2.);
            results[i][2][2] = results[i][2][1] * (1 + 4/3. * step_size * b + sqrt(4/3.) * sqrt_step_size * sigma * (random_normal[i][1] + random_normal[i][2]) / sqrt(2));
            results[i][2][3] = results[i][2][2] * (1 + 4/3. * step_size * b + sqrt(4/3.) * sqrt_step_size * sigma * random_normal[i][3]);
            // U^2
            double u2_1 = (sqrt(2) * random_normal[i][0] + random_normal[i][1]) / sqrt(3.);
            double u2_2 = (random_normal[i][2] + sqrt(2) * random_normal[i][3]) / sqrt(3.);
            results[i][1][1] = results[i][1][0] * (1 + 2. * step_size * b + sqrt(2) * sqrt_step_size * sigma * u2_1);
            results[i][1][2] = results[i][1][1] * (1 + 2. * step_size * b + sqrt(2) * sqrt_step_size * sigma * u2_2);
            // U^1
            results[i][0][1] = results[i][0][0] * (1 + 4. * step_size * b + 2 * sqrt_step_size * sigma * (u2_1 + u2_2));
        }
    }
    return results;
}
