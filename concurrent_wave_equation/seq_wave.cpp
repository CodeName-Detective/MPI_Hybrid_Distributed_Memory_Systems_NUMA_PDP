#include <iostream>
#include <cmath>
#include <cstring>

#define MAXPOINTS 1000
#define MAXSTEPS 1000
#define MINPOINTS 20
#define PI 3.14159265


void init_param(int& total_steps, int& total_points){
    std::cout<<"Enter total number of points vibrating along the string [20-1000]:";
    std::cin>>total_points;
    std::cout<<"Enter total number of time steps [1-1000]:";
    std::cin>>total_steps;
}

void init_line(const int& total_points, double* const values, double* const old_values){
    double x, fac{2.0*PI}, k{0.0}, tmp{(double)total_points-1};

    for(int i=0; i<total_points; ++i,++k){
        x=k/tmp;
        values[i] = std::sin(fac*x);
    }

    std::memcpy(old_values, values, (total_points)*sizeof(double));
}


void do_math(int idx, double* const new_values, const double* const values, const double* const old_values){
    double dtime{0.3}, c{1.0}, dx{1.0}, tau, sqtau;
    tau = (c * dtime / dx);
    sqtau = tau * tau;
    new_values[idx] = (2.0*values[idx]-old_values[idx]+ (sqtau*(values[idx-1]+ values[idx+1] - (2.0*values[idx]))));
}


void update(const int& total_steps, const int&  total_points, double* const new_values, double* const values, double* const old_values){
    for(int i=0; i<total_steps; ++i){
        for(int j=0; j<total_points;++j){
            if((j==0) || (j==total_points-1)){
                new_values[j] = 0.0;
            }
            else{
                do_math(j, new_values, values, old_values);
            }
        }

        std::memcpy(old_values, values, (total_points)*sizeof(double));
        std::memcpy(values, new_values, (total_points)*sizeof(double));
    }
}


void print_final( const int& total_points,const double* const values){
    for(int i=0; i<total_points; ++i){
        std::cout<<values[i]<<", ";
        if((i+1)%10 == 0){std::cout<<std::endl;}
    }
}

int main(int argc, char* argv[]){
    int total_steps, total_points;
    double *values, *old_values, *new_values;

    init_param(total_steps, total_points);

    values = new double[total_points];
    old_values = new double[total_points];
    new_values = new double [total_points];

    init_line(total_points, values, old_values);
    update(total_points, total_points, new_values, values, old_values);
    print_final(total_points, values);

    delete[] values;
    delete[] old_values;
    delete[] new_values;
    return 0;
}