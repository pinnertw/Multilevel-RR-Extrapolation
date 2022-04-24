// Estimator constructed with structural_params and multilevel_params.
#include <bits/stdc++.h>
#include "define.h"
#include "estimator.h"
using namespace std;

estimator::estimator(method_type method_, structural_params &sp_, multilevel_params &mlp_)
    : method(method_), sp(sp_), mlp(mlp_)
{
    M = 2;
}

void estimator::auto_tune(double epsilon){
    double cost_ = DBL_MAX;
    int M_opt = 1;
    double cost_loc;
    for (auto M_=2; M_<11; M_++){
        M = M_;
        init(epsilon);
        cost_loc = cost();
        if (cost_ > cost_loc)
        {
            cost_ = cost_loc;
            M_opt = M_;
        }
    }
    M = M_opt;
}

void estimator::init(double epsilon){
    sp.R = R_star(epsilon);
    sp.init_n(M);
    mlp.init_weights(sp.R, sp.n);
    init_T();
    sp.h = h_star(epsilon);
    init_q_and_N(epsilon);
}

// Theorem 3.6
void estimator::init_q_and_N(double epsilon){
    // require T, mlp.theta, sp.h, mlp.beta, n
    sp.q = vd (sp.R, 0.);
    double theta_h_beta2 = mlp.theta * pow(sp.h, mlp.beta / 2.);
    sp.q[0] = 1 + theta_h_beta2;
    double sum_q = sp.q[0];
    double sum_N = sp.q[0];

    // qj for j in 2,...,R
    for (auto i=1; i<sp.R; i++){
        sp.q[i] = theta_h_beta2 * unitary_variance(i) * pow(unitary_cost(i), -0.5);

        sum_q += sp.q[i];
        sum_N += theta_h_beta2 * unitary_variance(i) * pow(unitary_cost(i), 0.5);
    }
    // Normalization
    for (auto i=0; i<sp.R; i++){
        sp.q[i] /= sum_q;
    }

    // N
    if(method == Multilevel_RR) sp.N = (1 + 1 / (2. * mlp.alpha * sp.R)) * mlp.varY0 * sum_N * sum_q / epsilon / epsilon;
    if(method == Multilevel_MC) sp.N = (1 + 1 / (2. * mlp.alpha)) * mlp.varY0 * sum_N * sum_q / epsilon / epsilon;
}

int estimator::unitary_cost(int i){
    if (mlp.sim_type == diffusion) return sp.n[i-1] + sp.n[i];
    // Otherwise nested
    return sp.n[i];
}

double estimator::unitary_variance(int i){
    if (mlp.sim_type == diffusion) return (pow(sp.n[i-1], -mlp.beta / 2.) + pow(sp.n[i], -mlp.beta / 2.));
    // Otherwise nested
    return (pow(1 / double(sp.n[i-1]) - 1 / double(sp.n[i]), mlp.beta/2));
}

double estimator::kappa(){
    double sum = sp.q[0];
    for (auto i=1; i < sp.R; i++){
        sum += sp.q[i] * unitary_cost(i);
    }
    return sum;
}

double estimator::cost(){
    return sp.N * kappa() / sp.h;
};

// Proposition 3.9, bias parameter optimization.
double estimator::h_star(double epsilon){
    if (method == Multilevel_RR){
        double n_fact = 0.;
        double alphaR = mlp.alpha * sp.R;
        h_inverse = ceil(pow(1 + 2 * alphaR, 1 / 2. / alphaR) *
            pow(epsilon , -1/alphaR) *
            pow(M, -(sp.R - 1)/2));
        return 1 / (double)h_inverse;
    }
    else if (method == Multilevel_MC){
        h_inverse = ceil(pow(1 + 2 * mlp.alpha, 1 / 2. / mlp.alpha) * 
            pow(epsilon / fabs(mlp.c1), -1/mlp.alpha) * 
            pow(M, -sp.R + 1));
        return 1/ (double) h_inverse;
    }
    else return 1.;
}

// Theorem 3.12
int estimator::R_star(double epsilon){
    if (method == Multilevel_RR){
        double a = pow(mlp.c_tilde, 1/mlp.alpha) * mlp.hmax;
        double b = 0.5 + log(a) / log(M);
        return ceil(b + sqrt(pow(b, 2.) + 2 * log(sqrt(1 + 4 * mlp.alpha)/epsilon)/mlp.alpha/log(M)));
    }
    else if (method == Multilevel_MC){
        return ceil(1 + log(pow(mlp.c_tilde, 1/mlp.alpha) * mlp.hmax)/log(M) + log(sqrt(1+2*mlp.alpha)/epsilon) / mlp.alpha / log(M));
    }
    else return 2;
}
 
void estimator::init_T(){
    // Require sp.R, mlp.weights defined
    sp.T = vvd (sp.R, vd(2, 0.));
    if (method == Multilevel_MC){
        for (auto i=1; i<sp.R; i++){
            sp.T[i][1] = 1;
            sp.T[i][0] = -1;
        }
        sp.T[0][0] = 1;
    }
    else if(method == Multilevel_RR){
        double W_total = mlp.w_total;
        for (auto i=1; i<sp.R; i++){
            W_total -= mlp.weights[i-1];
            sp.T[i][1] = W_total;
            sp.T[i][0] = -W_total;
        }
        sp.T[0][0] = 1;
    }
}

double RMSE(vd results, double real_value){
    double sum=0.;
    for (auto i: results){
        sum += pow(i - real_value, 2.);
    }
    sum /= (double)results.size();
    return sqrt(sum);
}

double bias(vd results, double real_value){
    double sum=0.;
    for (auto i: results){
        sum += i - real_value;
    }
    sum /= (double) results.size();
    return sum;
}

double var(vd results){
    double sum=0.;
    double sum_square=0.;
    for (auto i: results){
        sum += i;
        sum_square += i * i;
    }
    sum /= (double)results.size();
    sum_square /= (double)results.size();
    return sum_square - sum * sum;
}
