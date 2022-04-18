#include <bits/stdc++.h>
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
#if PRINT1
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
#if PRINT1
            cout << "Print weight coefficients" << endl;
            for (auto i = 0; i < R; i++){
                cout << weights[i] << endl;
            }
            cout << "w_tilde = " << w_tilde << endl;
#endif
        }
};

class euler_schema
{
    private:
        double step_size;
        double sqrt_step_size;
        int total_step;
        double X0;
        function<double(double)> b;
        function<double(double)> sigma;
        vd random_normal;
    public:
        euler_schema(double step_size_, int total_step_, double X0_, function<double(double)> b_, function<double(double)> sigma_){
            step_size = step_size_;
            sqrt_step_size = sqrt(step_size);
            total_step = total_step_;
            X0 = X0_;
            b = b_;
            sigma = sigma_;
        }

        void get_normal_distribution(int N){
            mt19937_64 generator;
            auto seed = chrono::system_clock::now().time_since_epoch().count();
            random_normal = vd(N*total_step, 0.);
            normal_distribution<double> dist(0., 1.);
            for (auto i=0; i < N*total_step; i++){
                random_normal[i] = dist(generator);
            }
        }

        vvd simulations(int N){
            get_normal_distribution(N);
            vvd results = vvd(N, vd(total_step + 1, 0.));
            for (auto i=0; i < N; i++){
                results[i][0] = X0;
                for (auto j = 1; j <= total_step; j++){
                    results[i][j] = results[i][j-1] + step_size * b(results[i][j-1]) + sqrt_step_size * sigma(results[i][j-1]) * random_normal[i * total_step + j];
                }
            }
#if PRINT1
            for (auto i=0; i < N; i++){
                for (auto j=0; j < N; j++){
                    cout << results[i][j] << ' ';
                }
                cout << endl;
            }
#endif
            return results;
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
    function<double(double)> b = [](double a) {return a * 5.;};
    function<double(double)> sigma = [](double a) {return 5.;};
    euler_schema(0.01, 10, 80, b, sigma).simulations(10);
    return 0;
}
