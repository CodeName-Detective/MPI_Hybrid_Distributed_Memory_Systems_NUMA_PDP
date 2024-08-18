#include <iostream>
#include <cmath>

#define LIMIT     2500000+11     /* Increase this to find more primes */
#define PRINT     100000      /* Print a line after this many numbers */

bool isprime(const int& n){
    int sqr_root;
    if(n>10){
        sqr_root = (int) std::sqrt(n);
        for(int i=3; i<=sqr_root; i+=2){
            if((n%i)==0){return false;}
        }
        return true;
    }
    else{
        return false;
    }
}



int main(int argc, char* argv[]){
    int prime_counter, recent_prime;
    std::cout<<"Starting! Numbers to be scanned= "<<LIMIT<<std::endl;

    prime_counter=4; //Assume the primes less than 10 (2,3,5,7) are counted here

    for(int i=11; i<=LIMIT; i+=2){
        if(isprime(i)){
            ++prime_counter;
            recent_prime=i;
        }
    }

    std::cout<<"Done. Largest prime is "<<recent_prime<<". Total number of primes are "<<prime_counter<<"."<<std::endl;
    return 0;
}