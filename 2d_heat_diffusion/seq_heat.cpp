#include <iostream>
#include <fstream>
#include <iomanip>

#define NXPROB 100
#define NYPROB 100


struct Params {
    float x;
    float y;
    int nts;
};


void intidat(const int nx, const int ny, float (* const u1)[NYPROB]){
    for(int ix=0; ix<nx; ++ix){
        for(int iy=0; iy<ny; ++iy){
            u1[ix][iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
        }
    }
}


void prtdat(const int nx, const int ny, const float (* const u1)[NYPROB], const char* fnam){
    std::ofstream fout;
    fout << std::fixed << std::setprecision(3);
    fout.open(fnam, std::ios::out);

    for(int iy=ny-1; iy>=0; --iy){
        for(int ix=0; ix<nx; ++ix){
            fout<<u1[ix][iy];
            if(ix != nx-1){
                fout<<" ";
            }
            else{
                fout<<"\n";
            }
        }
    }
    fout.close();
}


void update(const int nx, const int ny, const float (* const u1)[NYPROB], float (* const u2)[NYPROB], const Params& parms){
    for(int ix=1; ix<nx-1; ++ix){
        for(int iy=1; iy<ny-1; ++iy){
            u2[ix][iy] = u1[ix][iy]+ (parms.x*(u1[ix+1][iy]+u1[ix-1][iy]-2*u1[ix][iy])) + (parms.y*(u1[ix][iy+1]+u1[ix][iy-1]-2*u1[ix][iy]));
        }
    }
}

int main(int argc, char* argv[]){

    Params parms = {0.1, 0.1, 100};
    float u[2][NXPROB][NYPROB];

    std::cout<<"Starting serial version of 2D heat example..."<<std::endl;
    std::cout<<"Using "<<NXPROB<<"x"<<NYPROB<<" grid"<<std::endl;

    std::cout<<"Initializing grid and creating input file:"<<std::endl;
    intidat(NXPROB, NYPROB, &u[0][0]); // In this context u[0] will decay to &u[0][0] automatically in function call. Just like int arr[5]; 
    prtdat(NXPROB, NYPROB, &u[0][0], "seq_initial.dat");

    for(int ix=0; ix<NXPROB; ++ix){
        u[1][ix][0] = u[0][ix][0];
        u[1][ix][NYPROB-1] = u[0][ix][NYPROB-1];
    }

    for(int iy=0; iy<NYPROB; ++iy){
        u[1][0][iy] = u[0][0][iy];
        u[1][NXPROB-1][iy] = u[0][NXPROB-1][iy];
    }

    std::cout<<"Iterating over "<<parms.nts<<" time steps..."<<std::endl;
    int iz=0;

    for(int it=1; it<=parms.nts; ++it){
        update(NXPROB, NYPROB, &u[iz][0], &u[1-iz][0], parms);
        iz = 1-iz;
    }

    prtdat(NXPROB, NYPROB, &u[1][0], "seq_final.dat");
    std::cout<<"output file created!";
    return 0;
}