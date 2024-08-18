#include <iostream>
#include <random>


#define MASTER 0
#define DARTS 32768
#define ROUNDS 1024
#define sqr(x) ((x)*(x))


double dboard(int darts){
    double x_coord, y_coord, pi, random_num, temp;
    int score = 0;

    std::mt19937 random_generator_engine;
    std::uniform_real_distribution<double> float_distribution(0.0, 1.0);

    for(int n=0; n<darts; ++n){
        temp = float_distribution(random_generator_engine);
        x_coord = (2.0 * temp) - 1.0;

        temp = float_distribution(random_generator_engine);
        y_coord = (2.0 * temp) - 1.0;

        /* if dart lands in circle, increment score */
        if ((sqr(x_coord) + sqr(y_coord)) <= 1.0){
            ++score;
        }
    }

    pi = 4.0 * (double)score/(double)darts;
    return pi;
}


int main(int argc, char* argv[]){
    double pi, avg_pi=0;

    for(int i=0; i<ROUNDS; ++i){
        pi = dboard(DARTS);
        avg_pi = ((avg_pi*i)+pi)/(i+1);
    }

    std::cout<<"The calculated pi value is: "<<avg_pi<<std::endl;
  
    return 0;
}