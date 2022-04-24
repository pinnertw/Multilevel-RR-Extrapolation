#pragma once
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
        euler_scheme(double step_size_, double X0_, double b_, double sigma_, double T_);
        void reset_step(double step_size_);
        void get_normal_distribution(int N);
        vvd simulations(int N);
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
        nested_monte_carlo(double X0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_);

        void set_n(int n1_, int n2_);
        void get_normal_distribution(int N);
        vvd simulations(int N);
};
