#include "define.h"
using namespace std;

enum method_type {Multistep_RR, Multilevel_MC};
enum simulation_type {diffusion, nested};

class structural_params{
    public:
        // q
        vd q;
        // pi0
        double h;
        vi n;
        int R;
        int N;

        // Method type
        method_type method; 

        // Usage
        vvd T;

        structural_params(method_type);
        void init_n(int M);
};

class multilevel_params{
    public:
        // Given: alpha/beta
        double alpha; // weak error rate, exponant in polynome.
        double beta; // Strong approximation error, L2 error <= V1 h^beta

        // Given or to calculate: V1, var(Y0), theta), TODO
        double V1;
        double varY0;
        double theta;
        simulation_type sim_type;
        double hmax;

        double c_tilde=1;
        double c1=1;

        // weights
        vd weights;
        double w_total;

        vb dp_exist;
        vd dp_exponant_array;

        multilevel_params(double alpha_, double beta_, double V1_, double varY0_, simulation_type sim_type_);
        double dp_exponant(int i, int n_i);
        void init_weights(int R, vi& n);
};
