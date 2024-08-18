#include <iostream>
#include <random>

#include "mpi.h"

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
    double pi, avg_pi=0, global_avg_pi=0;
    double* buffer;

    int comm_size, task_rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_size( MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank( MPI_COMM_WORLD, &task_rank);

    int local_rounds = ROUNDS/comm_size;

    for(int i=0; i<local_rounds; ++i){
        pi = dboard(DARTS);
        avg_pi = ((avg_pi*i)+pi)/(i+1);
    }


    
    MPI_Send( &avg_pi, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);

    if(task_rank==MASTER){
        buffer = new double[comm_size];
        for(int i=0; i<comm_size; ++i){
            MPI_Recv(buffer+i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Recv(&buffer[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for(int i=0; i<comm_size; ++i){
            global_avg_pi = ((global_avg_pi*i)+buffer[i])/(i+1);
        }

        std::cout<<"The calculated pi value is: "<<global_avg_pi<<std::endl;

        delete[] buffer;
    }

    MPI_Finalize();
  
    return 0;
}