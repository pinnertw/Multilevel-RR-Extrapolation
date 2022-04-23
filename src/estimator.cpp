#include <bits/stdc++.h>
#include "define.h"
#include "structural_params.cpp"
using namespace std;

class estimator{
    private:
        method_type method;
        structural_params sp;
        multilevel_params mlp;
        int M;
    public:
        estimator(method_type method_, structural_params sp_, multilevel_params mlp_)
            : method(method_), sp(sp_), mlp(mlp_)
        {
            M = 2;
        }

        void init_T(){
            // Require sp.R, mlp.weights defined
            int R = sp.R;
            sp.T = vvd (sp.R, vd(2, 0.));
            if (method == Multilevel_MC){
                for (auto i=1; i<R; i++){
                    sp.T[i][1] = 1;
                    sp.T[i][0] = -1;
                }
                sp.T[0][0] = 1;
            }
            else if(method == Multilevel_RR){
                double W_total = mlp.w_total;
                for (auto i=1; i<R; i++){
                    W_total -= mlp.weights[i-1];
                    sp.T[i][1] = W_total;
                    sp.T[i][0] = -W_total;
                }
                sp.T[0][0] = 1;
            }
        }

        void init_q(){
            // Thm 3.6
            // require T, mlp.theta, sp.h, mlp.beta, n
            sp.q = vd (sp.R, 0.);
            double sum_q = 0.;

            double theta_h_beta2 = mlp.theta * pow(sp.h, mlp.beta / 2.);

            // q1
            sp.q[0] = 1 + theta_h_beta2;
            sum_q += sp.q[0];
            // qj for j in 2,...,R
            for (auto i=1; i<sp.R; i++){
                sp.q[i] = theta_h_beta2 * 
                    (fabs(sp.T[i][0]) * pow(sp.n[i-1], -mlp.beta / 2.) + fabs(sp.T[i][1]) * pow(sp.n[i], -mlp.beta / 2.))*
                    pow(sp.n[i-1] + sp.n[i], -0.5);
                sum_q += sp.q[i];
            }
            // Normalization
            for (auto i=0; i<sp.R; i++){
                sp.q[i] /= sum_q;
            }
        }

        // Proposition 3.9, bias parameter optimization.
        double h_star(double epsilon){
            if (method == Multilevel_RR){
                double n_fact = 0.;
                double alphaR = mlp.alpha * sp.R;
                return pow(1 + 2 * alphaR, -1 / 2. / alphaR) *
                    pow(epsilon , 1/alphaR) *
                    pow(M, (sp.R - 1)/2);
            }
            else if (method == Multilevel_MC){
                return pow(1 + 2 * mlp.alpha, -1 / 2. / mlp.alpha) * 
                    pow(epsilon / fabs(mlp.c1), 1/mlp.alpha) * 
                    pow(M, sp.R - 1);
            }
            else return 1.;
        }

        // Theorem 3.12
        int R_star(double epsilon){
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

        void tuning(double epsilon){
            sp.h = h_star(epsilon);
            sp.R = R_star(epsilon); 
            sp.init_n(M);
            // Compute N, complexity TODO
            // 
        }
};
