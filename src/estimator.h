#pragma once
// Estimator constructed with structural_params and multilevel_params.
#include <bits/stdc++.h>
#include "define.h"
#include "structural_params.h"
#include "models.h"
using namespace std;

class estimator{
    public:
        method_type method;
        structural_params sp;
        multilevel_params mlp;
        int M;
        int h_inverse;
    public:
        estimator(method_type method_, structural_params& sp_, multilevel_params& mlp_);
        void init(double epsilon);
        void auto_tune(double epsilon);

        // Theorem 3.6
        void init_q_and_N(double epsilon);
        int unitary_cost(int i);
        double unitary_variance(int i);
        double kappa();
        double cost();
        // Proposition 3.9, bias parameter optimization.
        double h_star(double epsilon);
        // Theorem 3.12
        int R_star(double epsilon);
        void init_T();
};

double RMSE(vd, double);
double bias(vd, double);
double var(vd);
