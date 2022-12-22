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
void   read_counts( std::string cfile, std::vector<int>& counts );
void   simulator( std::vector<int>& counts );
void   stationary_vector( std::vector<double>& beta, std::vector<double>& phi, std::vector<double>& sd );
double log_likelihood( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts );
void   mh_beta( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, double & llk, std::vector<double> & tunning_beta, std::vector<int> & ap_beta );
void   mh_phi( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, double & llk, std::vector<double> & tunning_phi, std::vector<int> & ap_phi );
void   mcmc(std::string ofile, std::vector<int>& counts);
double rnorm(double mean, double sd);
double runif(double min, double max);
double rtrian(double min, double peak, double max);


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

// log log_likelihood
// counts is a vector that counts each codon in a protein coding sequence
double log_likelihood(std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts ){
 
  // empty stationary vector    
  std::vector<double> sd(64);

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
      }
    }
  }
  
  // calculating the log log likelihood
  double llk = 0.0;
  for (int i=0; i<64; ++i){
    llk += counts[i]*log(sd[i]/sum);
  }
  
  // returning the log likelihood
  return llk;
    
}



// metropolis-hastings algorithm on phi
// llk is the log likelihood
// tunning_phi is the variance of the normal proposal
// ap_phi is the number of accepted changes during the mcmc: optimal when ap_beta/Total_iterates~0.234
void mh_phi( std::vector<double>& beta,std::vector<double>& phi, std::vector<int>& counts, double & llk,
             std::vector<double> & tunning_phi, std::vector<int> & ap_phi ){

  // performing a mh step per phi
  // first phi is fixed = 1.0 (that is why i=1 and not 0)
  for (int i=1; i<K; ++i){

    // Normal proposal: proposes a new phi that is strictly positive
    double phip = rnorm(phi[i],tunning_phi[i]);
    if ( phip < 0.0 ) { continue; }
  
    // updating the new phi onto the phi vector
    std::vector<double> phi1=phi;
    phi1[i] = phip;

    // calculating the new likelihood
    double llk1 = log_likelihood( beta, phi1, counts);
  
    // Accept-reject step: proposal ratio = 1 
    double acceptance_prob = exp( llk1 - llk );
    //std::cout << "llk=" << llk << " lik1=" << llk1 << " ratio=" << acceptance_prob  << " phi=" << phip <<  "\n\n";
    if (runif(0,1) < acceptance_prob) {
      phi       = phi1;
      llk       = llk1; 
      ap_phi[i] += 1;    
    } 

  }   
}


// metropolis-hastings algorithm on beta
// llk is the log likelihood
// tunning_beta is the variance of the normal proposal
// ap_beta is the number of accepted changes during the mcmc: optimal when ap_beta/Total_iterates~0.234
void mh_beta( std::vector<double>& beta,std::vector<double>& phi, std::vector<int>& counts, double & llk,
              std::vector<double>  & tunning_beta, std::vector<int> & ap_beta ){
  
  // performing a mh step per beta
  // first beta is fixed = 1.0 (that is why i=1 and not 0)
  for (int i=1; i<4; ++i){

    // Normal proposal: proposes a new phi that is strictly positive
    double betap = rnorm(beta[i],tunning_beta[i]);
    if ( betap < 0.0 ) { continue; }
  
    // updating the new phi onto the phi vector
    std::vector<double> beta1=beta;
    beta1[i] = betap;

    // calculating the new likelihood
    double llk1 = log_likelihood( beta1, phi, counts);
  
    // Accept-reject step: proposal ratio = 1 
    double acceptance_prob = exp( llk1 - llk );
    //std::cout << "llk=" << llk << " lik1=" << llk1 << " ratio=" << acceptance_prob  << " phi=" << phip <<  "\n\n";
    if (runif(0,1) < acceptance_prob) {
      beta      = beta1;
      llk       = llk1; 
      ap_beta[i] += 1;    
    } 

  }   

}

// MCMC sampler
// ofile is the output file
// counts is a vector of observed counts
void mcmc (std::string ofile, std::vector<int>& counts){

  // setting the number of generations
  //int I = 10000;
  int I = 1000000;

  // opening the output file 
  std::ofstream myfile;
  myfile.open(ofile);

  // setting the header
  myfile << "Gen llk ";
  for (int i=0; i<4; ++i){ myfile << "beta" << i << " "; }
  for (int i=0; i<K; ++i){ myfile << "phi" << i << " "; }
  myfile << "\n";

  // defining the initial parameter values 
  std::vector<double> beta(4,1);
  beta[1] = runif(0,2);
  beta[2] = runif(0,2);
  beta[3] = runif(0,2);

  std::vector<double> phi(K,1);
  for (int i=1; i<K; ++i){ 
    phi[i] += rtrian(0.5,1,1.5); 
  }
  
  // calculating the likelihood
  double llk = log_likelihood(beta, phi, counts); 
  std::cout << "\n  Initial likelihood: " <<llk<< "\n";

  // setting the tunning and acceptance probability vectors
  // these will be use to tune the proposal so we have optimal mcmc 
  std::vector<double> tunning_beta(4,1.0);
  std::vector<int> ap_beta(4,0);
  std::vector<double> tunning_phi(K,1.0);
  std::vector<int> ap_phi(K,0);

  // running a single MCMC chain
  for (int s=0; s<I; ++s){

    // metropolis-hasting steps
    mh_beta( beta, phi, counts, llk, tunning_beta, ap_beta );
    mh_phi(  beta, phi, counts, llk, tunning_phi, ap_phi );

    // tunning the proposal 
    if (s%100==0){
      
      // base frequency
      for (int i=1; i<4; ++i){  
        tunning_beta[i] /= (2.0*0.234/(0.234+1.0*ap_beta[i]/(s+1))); 
      }

      // fitness coefficients
      for (int i=1; i<K; ++i){ 
        tunning_phi[i] /= (2.0*0.234/(0.234+1.0*ap_phi[i]/(s+1))); 
      }

      // saving the sample to the output file
      if (s%1000==0){

        myfile << s << " " << llk << " ";
        std::cout <<  s << " " << llk << " \n";

        // base frequency
        for (int i=0; i<4; ++i){  
          myfile << beta[i] << " ";
        }

        // fitness coefficients
        for (int i=0; i<K; ++i){ 
          myfile << phi[i] << " ";
        }
        myfile << "\n";

      }

    }
  }

  // closing the file
  myfile.close();
  std::cout << "  The analyses are finished! Check mcmc.txt file.\n\n";
}

// reads a vector of counts
// cfile is the path to the count file
// counts is a vector of integers that will carry the observed counts 
void read_counts( std::string cfile, std::vector<int>& counts ) {

  //reads the file from the terminal
  std::ifstream myfile (cfile); 

  //sanity check for counts file existence
  if ( myfile ) {    
    std::cout << "\n" << "Counts found in the " << cfile << " file:" << "\n";
    
    // specifies an integer and sets an indexing
    int count;
    int p = 0;

    // reads each cunt into the counts vector  
    while ( myfile >> count ) {
        counts[p] = count;
        p += 1;
    }

    //closes the file
    myfile.close(); 

    // prints the counts
    for (int i=0; i<counts.size(); i++){ 
      std::cout << "[" << i << "] " << counts[i] << "\n";
    }
    std::cout << "\n"; // adds empty line at the end of output

    // some check ups here...
    std::cout << "  \n\n  Some check ups needed here: reads counts from file; length of counts, etc.\n  The analyses will start now.\n\n";

  } else {
    //prints if file not opened
    std::cout << "  Unable to open " << cfile << "."; 
    
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

