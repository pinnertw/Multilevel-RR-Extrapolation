#include "define.h"
using namespace std;

class euler_scheme{
    private:
        double step_size;
        double sqrt_step_size;
        double T;
        double X0;
        double b;
        double sigma;
        int total_step;

        vd random_normal;
    public:
        euler_scheme(double step_size_, double X0_, double b_, double sigma_, double T_):
            step_size(step_size_), X0(X0_), b(b_), sigma(sigma_), T(T_){
                sqrt_step_size = sqrt(step_size);
                total_step = (int)(T / step_size);
            }
        void reset_step(double step_size_){
            step_size = step_size_;
            sqrt_step_size = sqrt(step_size);
            total_step = (int)(T / step_size);
        }
        void get_normal_distribution(int N){
            mt19937_64 generator;
            auto seed = chrono::system_clock::now().time_since_epoch().count();
            random_normal = vd(N * total_step, 0.);
            normal_distribution<double> dist(0., 1.);
            for (auto i=0; i<N*total_step; i++){
                random_normal[i] = dist(generator);
            }
        }
        vvd simulations(int N){
            get_normal_distribution(N);
            vvd results = vvd(N, vd(total_step + 1, 0.));
            for (auto i=0; i <N; i++){
                results[i][0] = X0;
                for (auto j=1; j<=total_step; j++){
                    results[i][j] = results[i][j-1] + step_size * b * results[i][j-1] + sqrt_step_size * sigma * results[i][j-1] * random_normal[i * total_step + j];
                }
            }
            return results;
        }
};

class nested_monte_carlo{
    private:
        double X0;
        double r;
        double sigma;
        double K1;
        double K2;
        double T1;
        double T2;
        int n1;
        int n2;

        vd random_normal;
    public:
        nested_monte_carlo(double X0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_)
        : X0(X0_), r(r_), sigma(sigma_), K1(K1_), K2(K2_), T1(T1_), T2(T2_){
        }
        void set_n(int n1_, int n2_){
            n1 = n1_;
            n2 = n2_;
        }

        void get_normal_distribution(int N){
            mt19937_64 generator;
            auto seed = chrono::system_clock::now().time_since_epoch().count();
            random_normal = vd(N * (n2+1), 0.);
            normal_distribution<double> dist(0., 1.);
            for (auto i=0; i < N*(n2+1); i++){
                random_normal[i] = dist(generator);
            }
        }
        void simulations(int N){
            get_normal_distribution(N);
            vvd results = vvd(N, vd((n2+1), 0.));
            for (auto i=0; i<N; i++){
                // results[i][0] = S_T1, results[i][1 ~ n2] = S_T2
                results[i][0] = X0 * exp((r - 0.5 * sigma * sigma) * T1 + sigma * sqrt(T1) * random_normal[i*(n2+1)]);
                for (auto j=1; j<n2 + 1; j++){
                    results[i][j] = results[i][0] * exp((r - 0.5 * sigma * sigma) * (T2 - T1) + sigma * sqrt(T2 - T1) * random_normal[i * (n2 + 1) + j]);
                }
        }
}

