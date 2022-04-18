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
                cout << weights[i] << ' ';
            }
            cout << endl;
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
        function<double(double &)> b;
        function<double(double &)> sigma;
        vd random_normal;
    public:
        euler_schema(double step_size_, int total_step_, double X0_, function<double(double &)> b_, function<double(double &)> sigma_){
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

class nested_monte_carlo
{
    private:
        double s0;
        double r;
        double sigma;
        double K1;
        double K2;
        double T1;
        double T2;
        // Normal dist.
        vd random_normal;
    public:
        nested_monte_carlo(double s0_, double r_, double sigma_, double K1_, double K2_, double T1_, double T2_){
            s0 = s0_;
            r = r_;
            sigma = sigma_;
            K1 = K1_;
            K2 = K2_;
            T1 = T1_;
            T2 = T2_;
        }
        void get_normal_distribution(int N){
            mt19937_64 generator;
            auto seed = chrono::system_clock::now().time_since_epoch().count();
            random_normal = vd(N*2, 0.);
            normal_distribution<double> dist(0., 1.);
            for (auto i=0; i < N*2; i++){
                random_normal[i] = dist(generator);
            }
        }

        void simulations(int N){
            get_normal_distribution(N);
            vvd results = vvd(N, vd(2, 0.));
            for (auto i=0; i < N; i++){
                // TODO: Check how many simulation per one S_T1
                results[i][0] = s0 * exp((r - 0.5 * sigma * sigma) * T1 + sigma * sqrt(T1) * random_normal[2 * i]);
                results[i][1] = results[i][0] * exp((r - 0.5 * sigma * sigma) * (T2 - T1) + sigma * sqrt(T2 - T1) * random_normal[2 * i + 1]);
            }
        }
};

class geometric_brownian_payoff
{
    private:
        double x0;
        double r;
        double sigma;
        double K;
        double T;
        function<double(vd &)> payoff_func;
        double alpha;
        double beta;
        // Other params for different options
        double lambda;
        double B;
        double T2;
        double K2;
        // TODO
        // Real value, Variance, Complexity, Cost, auto-tuning
    public:
        geometric_brownian_payoff(int type_, double x0_, double r_, double sigma_, double K_, double T_, 
                double extra_param_ = 1.1, double extra_param2_ = 2.2){
            x0 = x0_;
            r = r_;
            sigma = sigma_;
            K = K_;
            T = T_;
            /* type_:
             * 1 -> BS vanilla call
             * 2 -> Lookback option
             * 3 -> Barrier option
             * 4 -> Compound option
             */
            if (type_ == 1){
                alpha = 1.;
                beta = 1.;
                payoff_func = [&](vd& results) {
                    int total_steps = results.size();
                    if (results[total_steps-1] > K){
                        return exp(-r * T) * (results[total_steps-1] - K);
                    }
                    else{
                        return 0.;
                    }
                };
            }
            else if (type_ == 2){
                alpha = 0.5;
                beta = 1.;
                lambda = extra_param_;
                payoff_func = [&](vd& results){
                    double x_min = 1e9;
                    for (double x : results){
                        if(x < x_min) x_min = x;
                    }
                    return exp(-r * T) * max(results[results.size()-1] - lambda * x_min, 0.);
                };
            }
            else if (type_ == 3){
                alpha = 0.5;
                beta = 0.5;
                B = extra_param_;
                payoff_func = [&](vd& results){
                    bool passed = false;
                    for (double x : results){
                        if (x > B) {
                            passed = true;
                            break;
                        }
                    };
                    if (!passed) return exp(-r * T) * max(results[results.size()-1] - K, 0.);
                    else return 0.;
                };
            }
            else if (type_ == 4){
                alpha = 1.;
                beta = 1.;
                T2 = extra_param_;
                K2 = extra_param2_;
                payoff_func = [&](vd& results){
                    // TODO
                    return 0.;
                };
            }
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
    function<double(double&)> b = [](double a) {return a * 5.;};
    function<double(double&)> sigma = [](double a) {return 5.;};
    euler_schema schema = euler_schema(0.01, 10, 80, b, sigma);
    schema.simulations(10);
    geometric_brownian_payoff(1, 100, 0.06, 0.4, 80, 1);
    geometric_brownian_payoff(1, 100, 0.06, 0.4, 80, 2);
    geometric_brownian_payoff(1, 100, 0.06, 0.4, 80, 3);
    geometric_brownian_payoff(1, 100, 0.06, 0.4, 80, 4);
    test();
    return 0;
}
