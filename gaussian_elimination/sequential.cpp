#include <algorithm>
#include <iostream>
#include <random>
#include <vector>


void print_matrix(const float* matrix, int dim){
    for(int i=0; i<dim; ++i){
        for(int j=0; j<dim; ++j){
            std::cout<<matrix[i*dim+j]<<" ";
        }
        std::cout<<std::endl;
    }
}


int main(){
    std::mt19937 mt(123);
    std::uniform_real_distribution dist(1.0f, 2.0f);

    std::vector<float> m;
    const int dim = 1 << 3;
    std::generate_n(std::back_inserter(m), dim*dim, [&] {return dist(mt);});

    // Gaussian Elimination
    for(int row=0; row<dim; ++row){
        float pivot = m[row*dim+row];

        // Divide the row to the right of the pivot by pivot 
        for(int col=row; col<dim; ++col){
            m[row*dim+col] =  m[row*dim+col]/pivot;
        }

        // Eliminate Pivot column from the remaining rows
        for(int eliminate_row = row+1; eliminate_row<dim; ++eliminate_row){
            float scale = m[eliminate_row*dim+row];

            //remove the pivot
            for(int col=row; col<dim; ++col){
                m[eliminate_row*dim+col] = m[eliminate_row*dim+col] - (m[row*dim+col]*scale); 
            }
        }
    }

    print_matrix(m.data(), dim);

    return 0;
}