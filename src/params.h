#ifndef PARAMS
#define PARAMS
#include <bits/stdc++.h>
using namespace std;

#define vvi vector<vector<int>>
#define vi vector<int>
#define vvd vector<vector<double>>
#define vd vector<double>
#define vb vector<bool>
#define d double

void test();
class params{
    private:
        int M;
        int R;
        int type;
        vi refiners;
        vd weights;
        double w_tilde;
        double alpha;
        double beta;
        vb dp_exist;
        vd dp_exponant_array;
    public:
        params(int M_, int R_, double alpha_, double beta_);
        void print_params();
        double dp_exponant(int i);
        void init_refiners();
        void init_weights();
};
class euler_schema{
    private:
        double step_size;
        double sqrt_step_size;
        int total_step;
        double X0;
        function<double(double &)> b;
        function<double(double &)> sigma;
        vd random_normal;
    public:
        euler_schema(double step_size_, int total_step_, double X0_, function<double(double &)> b_, function<double(double &)> sigma_);
        void get_normal_distribution(int N);
        vvd simulations(int N);
};
class nested_monte_carlo{
    private:
        double s0;
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
        nested_monte_carlo(double s0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_, int n1_, int n2_);
        void reset_n(int n1_, int n2_);
        void get_normal_distribution(int N);
        void simulations(int N);
};

#endif
