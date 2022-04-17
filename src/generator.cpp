#include <bits/stdc++.h>
using namespace std;


#define vvi vector<vector<int>>
#define vi vector<int>
#define vd vector<double>
#define vb vector<bool>

void test(){
    cout << 0 << endl;
}

class params
{
    private:
        int M;
        int R;
        int type;
        vi refiners;
        vd weights;
        double w_tilde;
        double alpha;
        double beta;

        // DP usage in ram
        vb dp_exist;
        vd dp_exponant_array;


    public:
        params(int M_, int R_, double alpha_, double beta_){
            // Init params with M as root, and level R
            // Then get the refiners with n_j = M^j;
            M = M_;
            R = R_;
            alpha = alpha_;
            beta = beta_;
            // Refiners
            init_refiners();
            // Weights
            dp_exist = vb(R, false);
            dp_exponant_array = vd(R, 0);
            weights = vd(R, 1.);
            init_weights();

            print_params();
        }

        void print_params(){
#if PRINT
            cout << "M, R, alpha, beta"<<endl;
            cout << M << ", " << R << ", " << alpha << ", " << beta << endl;
#endif
        }

        double dp_exponant(int i){
            if (dp_exist[i]) return dp_exponant_array[i];
            else{
                dp_exponant_array[i] = pow(refiners[i], alpha);
                dp_exist[i] = true;
                return dp_exponant_array[i];
            }
        }

        void init_refiners(){
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

        void init_weights(){
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
#if PRINT
            cout << "Print weight coefficients" << endl;
            for (auto i = 0; i < R; i++){
                cout << weights[i] << endl;
            }
            cout << "w_tilde = " << w_tilde << endl;
#endif
        }
};

class estimators
{
    private:
        double h;
    public:
        estimators(double h_){
            h = h_;
        }
};



int main(int argc, char ** argv)
{
    params param(5, 6, 0.5, 1.);
    estimators est(10);
    return 0;
}
