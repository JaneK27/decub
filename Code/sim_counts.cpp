//g++ codon.cpp -o DECUB -O3 -std=c++11

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>




// declaring the functions
void   simulator( std::vector<int>& counts );
void   stationary_vector( std::vector<double>& beta, std::vector<double>& phi, std::vector<double>& sd );
double log_likelihood( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts );
void   mcmc(std::string ofile, std::vector<int>& counts);
double rnorm(double mean, double sd);
double runif(double min, double max);
double rtrian(double min, double peak, double max);
void   rmapping( void );

//global variables

std::vector <int> mapping = {1,3,2,3,4,5,6,4,7,12,8,13,16,17,0,17,18,20,19,20,21,22,23,24,9,10,10,11,25,26,27,26,29,30,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,46,44,46,12,12,14,15,45,47,49,48,28,50,26,51};

// number of different phis in mapping
int K; 

int main(int argc, char** argv) {

  K=52;
  //K = atoi(argv[2]);

  // initializing the vector of counts
  std::vector<int> counts(64); 
  
  // simulating some data
  //simulator(counts);

  // input and output files  
  std::string cfile = argv[1];
  std::string ofile = argv[2];
  
  //K = atoi(argv[3]);
  std::cout << "\n" << "  K: " << K << "\n" << std::endl;
  //rmapping();

  // reads counts from a file
  read_counts (cfile, counts);
  // estimating the parameters
  mcmc (ofile, counts);
  
}


// stationary vector
// beta is the nucleotide composition vector
// phi is the codon fitness coefficient vector
void stationary_vector( std::vector<double>& beta, std::vector<double>& phi, std::vector<double>& sd){

  // calculating the stationary vector
  // sum is the normalization constant
  // g counts the number of Cs ang Gs in a certain codon 
  // c counts de codons
  double sum = 0.0;
  int c,g;
  
  // c1, c2 and c3 are the nucleotide bases at the 1st, 2nd and 3rd codon position
  for (int c1=0; c1<4; ++c1){ 
    for (int c2=0; c2<4; ++c2){
      for (int c3=0; c3<4; ++c3){
        c     = c1*16+c2*4+c3;
        g     = 1*((c1 == 1) | (c1 == 2)) + 1*((c2 == 1) | (c2 == 2)) + 1*((c3 == 1) | (c3 == 2));
        sd[c] = beta[c1]*beta[c2]*beta[c3]*phi[mapping[c]];        
        sum += sd[c];
        // std::cout << "sd[c]=" << sd[c] << " " <<  g <<  "\n";
      }
    }
  }
  
  // normalizing the stationary vector
  for (int i=0; i<64; ++i){
    sd[i]/=sum;
    //std::cout << "sd[" << i << "]=" << sd[i] << " k=" << sum << "\n\n";
  }

}

// simulates counts
void simulator ( std::vector<int>& counts ){
  
  // defining the initial parameter values 
  // random mapping 
  rmapping();
  
  // gc-bias
  double sgamma = runif(1,2);

  // fitness coefficients
  std::vector<double> sbeta(4);
  sbeta[0] = 1.0;
  sbeta[1] = rtrian(0.5,1.3,3);
  sbeta[2] = rtrian(0.5,1.3,3);
  sbeta[3] = rtrian(0.5,1.3,3);

  // base composition
  std::vector<double> spi(4);
  spi[0] = sgamma / (sgamma + sbeta[1] + sbeta[2] + (sbeta[3]*sgamma));
  spi[1] = (sbeta[1]*spi[0])/sgamma;
  spi[2] = (sbeta[2]*spi[0])/sgamma;
  spi[3] = sbeta[3]*spi[0];

  std::cout << "\n  Simulated parameters: \n\n";

  std::cout << "  pi[" << 0 << "] = " << spi[0] << "\n";
  std::cout << "  pi[" << 1 << "] = " << spi[1] << "\n";
  std::cout << "  pi[" << 2 << "] = " << spi[2] << "\n";
  std::cout << "  pi[" << 3 << "] = " << spi[3] << "\n\n";

  std::cout << "  gamma = " << sgamma << "\n\n";

  std::cout << "  beta[" << 0 << "] = " << sbeta[0] << "\n";
  std::cout << "  beta[" << 1 << "] = " << sbeta[1] << "\n";
  std::cout << "  beta[" << 2 << "] = " << sbeta[2] << "\n";
  std::cout << "  beta[" << 3 << "] = " << sbeta[3] << "\n\n";

  std::vector<double> sphi(K,1);
  for (int i=1; i<K; ++i){ 
    sphi[i] += rtrian(-1,0,2); 
    std::cout << "  phi[" << i << "] = " << sphi[i] << "\n";
  }
  std::cout << " " << "\n";
  // number of sites
  int S = 10000000;
  
  // calculating the stationary distribution 
  std::vector<double> sd(64);
  stationary_vector(sbeta, sphi, sd);

  // calculating the counts
  for (int i=0; i<64; ++i){
    counts[i] = round(sd[i]*S);
    // std::cout << "  count[" << i << "] = " << counts[i] << "\n";
    //std::cout << counts[i] << ", ";
  }

}
    
// simulates from the normal distribution
double rnorm(double mean, double sd){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::normal_distribution<double> distribution(mean,sd);
  double random = distribution(generator);
  return random;
}

// simulates from the uniform distribution in 0:1
double runif(double min, double max){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(min,max);
  double random = distribution(generator);
  return random;
}

// simulates from triangle
double rtrian(double min, double peak, double max){
  double U=runif(0,1);
  double fpeak=(peak - min)/(max - min);
  if (U < fpeak){
    return min + sqrt(U*(max - min)*(peak - min));
  }
  else{
    return max - sqrt((1 - U)*(max - min)*(max - peak));
  }
}

// simulates a random mapping
void rmapping( void ){

  if (K == 64){
    for (int i=0; i<K; ++i){
      mapping[i] = i;
    }

  } else {
    // first elements are just the indexes
    // this is to guarentee all indexes are used
    for (int i=0; i<K; ++i){
      mapping[i] = i;
    }    

    // now we random chose some indexes from 0:(K-1)
    std::random_device rd; 
    std::mt19937 generator(rd()); 
    std::uniform_int_distribution<int> distribution {0, K-1};

    for (int i=K; i<64; ++i){
      mapping[i] =  distribution(generator);
    }


    //for (int i=0; i<64; ++i){
    //  std::cout << "pos1 " << i  << " = " <<  mapping[i] << "\n";
    //}
    std::cout << "  New mapping: " << "\n";
    // we then shuffle the vector
    std::shuffle(begin(mapping)+1, end(mapping), generator);
    for (int i=0; i<64; ++i){
      std::cout << "  pos_" << i  << " = " <<  mapping[i] << "\n";
    }
    std::cout << "   " << "\n";
  }


}