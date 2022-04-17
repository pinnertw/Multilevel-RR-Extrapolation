#include <bits/stdc++.h>
using namespace std;

#define vvi vector<vector<int>>
#define vi vector<int>
#define vd vector<double>

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
        vd coeff_extrapolations;

    public:
        params(int M_, int R_, int type_){
            type = type_;
            if (type == 0){
                // Init params with M as root, and level R
                // Then get the refiners with n_j = M^j;
                M = M_;
                R = R_;
                refiners = vi(R);
                refiners[0] = 1;
                for (auto i = 1; i < R; i++){
                    refiners[i] = refiners[i-1] * M;
                }
            }
            else if (type == 1){
                // Init params with M as root, and level R
                // Then get the refiners with n_j = j;
                M = M_;
                R = R_;
                refiners = vi(R);
                refiners[0] = 1;
                for (auto i = 1; i < R; i++){
                    refiners[i] = i + 1;
                }
            }
        }
        void init_coeff_extrapolations(){

            if(type == 0){
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
    params param(5, 6, 0);
    estimators est(10);
    return 0;
}
