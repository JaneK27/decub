/* A Bayesian estimator  that quantifies the mutational and selective biases that drive codon preferences in animals.

    author: Ioanna Kotari & Rui Borges
    date: 24.12.2022
    contact: ioanna.kotari@hotmail.com */

// Compile command:
// g++ decub.cpp -o DECUB -O3 -std=c++11 

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

// declaring the functions
void   ReadControlFile(const std::string control_file,std::string& cfile, std::string& ofile);
void   read_counts( std::string cfile, std::vector<int>& counts );
void   read_mapping( std::string mfile, std::vector<int>& mapping );
void   simulator( std::vector<int>& counts );
double log_likelihood( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, std::vector<int>& mapping );
void   mh_beta( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, double & llk, std::vector<double> & tunning_beta, std::vector<int> & ap_beta, std::vector<int>& mapping );
void   mh_phi( std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, double & llk, std::vector<double> & tunning_phi, std::vector<int> & ap_phi, std::vector<int>& mapping );
void   mcmc(std::string ofile, std::vector<int>& counts, std::vector<int>& mapping );
double rnorm(double mean, double sd);
double runif(double min, double max);
double rtrian(double min, double peak, double max);


//global variables
int mcmc_gen, sample_freq, map;

// number of different phis in mapping
int K; 


int main(int argc, char** argv) {
  if ( argc > 1 ) {
    // gets the control file name  
    std::string control_file = argv[1]; 
    std::string cfile;
    std::string ofile; 
    // initialising the vector of counts
    std::vector<int> counts(64);

    ReadControlFile(control_file,cfile,ofile);
  
    // initialising the mapping vector
    std::vector<int> mapping;

    if (map==0) {
      std::cout << "\n\n  DECUB\n\n  DECUB is running with an input mapping \n";

      // check if mapping file exists
      if ( argc > 2 ) {    
      // reads mapping from a file
      std::string mfile = argv[2];
      read_mapping (mfile, mapping);
      
      // calculate K
      std::vector<int> nums=mapping;
      std::sort(nums.begin(), nums.end());
      K = std::unique(nums.begin(), nums.end()) - nums.begin();

      // reads counts from a file
      read_counts (cfile, counts);

      // prints out rest of parameters
      std::cout << " \n You specified the mcmc to run for "  << mcmc_gen <<  "  generations with a sampling frequency of " << sample_freq << ". \n\n";
    
      mcmc (ofile, counts, mapping);

      std::cout << "\n\n  DECUB has finished!\n\n";
      } else {
      //prints if file not opened
      std::cerr << " \n\n A file containing the mapping hasn't been given. \n Please provide a path to the mapping file or choose one of the existing ones. \n\n To choose one of the existing ones please change the mapping in the control file to 1 for the amino acid mapping, 2 for the Arthropoda mapping, or 3 for the Chordata one. \n\n  "; 
      } 

    // runs the moran model with selection
    } else if (map==1) {
      std::cout << "\n\n  DECUB\n\n  DECUB is running with the amino acid mapping \n";

      mapping = {1,2,1,2,3,3,3,3,4,5,4,5,6,6,0,6,7,8,7,8,9,9,9,9,4,4,4,4,10,10,10,10,11,12,11,12,13,13,13,13,14,14,14,14,15,15,15,15,16,17,18,17,5,5,5,5,19,20,21,20,10,22,10,22}; 
    
      // calculate K
      std::vector<int> nums=mapping;
      std::sort(nums.begin(), nums.end());
      K = std::unique(nums.begin(), nums.end()) - nums.begin();

      // reads counts from a file
      read_counts (cfile, counts);

      // prints out rest of parameters
      std::cout << " \n You specified the mcmc to run for "  << mcmc_gen <<  "  generations with a sampling frequency of " << sample_freq << ". \n\n";
    
      mcmc (ofile, counts, mapping);

      std::cout << "\n\n  DECUB has finished!\n\n";
    } else if (map==2) {
      std::cout << "\n\n  DECUB\n\n  DECUB is running with the Arthropoda mapping \n";

      mapping = {1,2,1,3,4,4,4,4,5,7,5,8,9,9,0,10,11,12,11,13,14,15,15,15,5,5,6,5,16,17,18,18,19,20,19,21,22,22,23,22,24,24,25,24,26,27,28,28,29,32,30,32,7,7,7,7,31,33,34,33,18,35,18,36};
    
      // calculate K
      std::vector<int> nums=mapping;
      std::sort(nums.begin(), nums.end());
      K = std::unique(nums.begin(), nums.end()) - nums.begin();

      // reads counts from a file
      read_counts (cfile, counts);

      // prints out rest of parameters
      std::cout << " \n You specified the mcmc to run for "  << mcmc_gen <<  "  generations with a sampling frequency of " << sample_freq << ". \n\n";
    
      mcmc (ofile, counts, mapping);
      std::cout << "\n\n  DECUB has finished!\n\n";
    } else if (map==3) {
      std::cout << "\n\n  DECUB\n\n  DECUB is running with the Chordata mapping \n";

      mapping = {1,3,2,3,4,5,6,4,7,12,8,13,16,17,0,17,18,20,19,20,21,22,23,24,9,10,10,11,25,26,27,26,29,30,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,46,44,46,12,12,14,15,45,47,49,48,28,50,26,51};
    
      // calculate K
      std::vector<int> nums=mapping;
      std::sort(nums.begin(), nums.end());
      K = std::unique(nums.begin(), nums.end()) - nums.begin();
      
      // reads counts from a file
      read_counts (cfile, counts);

      // prints out rest of parameters
      std::cout << " \n You specified the mcmc to run for "  << mcmc_gen <<  "  generations with a sampling frequency of " << sample_freq << ". \n\n";
    
      mcmc (ofile, counts, mapping);

      std::cout << "\n\n  DECUB has finished!\n\n";
    } else {
      std::cerr << "  Mapping is incorrectly specified. Please, check the control file.\n\n";
    }
  } else {
      std::cerr << "  You didn't provide any file. Please add a control file to run. \n\n";
  }
}

void ReadControlFile(const std::string control_file,std::string& cfile, std::string& ofile) {

  //reads the control file 
  std::string word;

  std::ifstream confile(control_file);
  if (!confile.is_open()) {
    std::cerr << "  Unable to open " << control_file << ".\n\n";
  }

  // updates parameters
  confile >> word >> cfile >> word >> ofile >> word >> map >> word >> mcmc_gen >> word >> sample_freq; 
  confile.close();

}


// log log_likelihood
// counts is a vector that counts each codon in a protein coding sequence
double log_likelihood(std::vector<double>& beta, std::vector<double>& phi, std::vector<int>& counts, std::vector<int>& mapping ){
 
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
             std::vector<double> & tunning_phi, std::vector<int> & ap_phi, std::vector<int>& mapping ){

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
    double llk1 = log_likelihood( beta, phi1, counts, mapping);
  
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
              std::vector<double>  & tunning_beta, std::vector<int> & ap_beta, std::vector<int>& mapping ){
  
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
    double llk1 = log_likelihood( beta1, phi, counts, mapping);
  
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
// mapping is the mapping vector 
void mcmc (std::string ofile, std::vector<int>& counts, std::vector<int>& mapping){

  // setting the number of generations
  //int I = 1000000;

  // opening the output file 
  std::ofstream myfile;
  myfile.open(ofile);

  // setting the header
  myfile << "Gen lnL ";
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

  std::cout << "  DECUB has started! \n\n";
  
  // calculating the likelihood
  double llk = log_likelihood(beta, phi, counts, mapping); 
  std::cout << "\n  Initial likelihood: " <<llk<< "\n";

  // setting the tunning and acceptance probability vectors
  // these will be use to tune the proposal so we have optimal mcmc 
  std::vector<double> tunning_beta(4,1.0);
  std::vector<int> ap_beta(4,0);
  std::vector<double> tunning_phi(K,1.0);
  std::vector<int> ap_phi(K,0);

  // running a single MCMC chain
  for (int s=0; s<mcmc_gen; ++s){

    // metropolis-hasting steps
    mh_beta( beta, phi, counts, llk, tunning_beta, ap_beta, mapping );
    mh_phi(  beta, phi, counts, llk, tunning_phi, ap_phi, mapping );

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
      if (s%sample_freq==0){

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
  std::cout << "  The analyses are finished! Check the " << ofile << " file.\n\n";
}

// reads a vector of counts
// cfile is the path to the count file
// counts is a vector of integers that will carry the observed counts 
void read_counts( std::string cfile, std::vector<int>& counts ) {

  //reads the file from the terminal
  std::ifstream myfile (cfile); 

  //sanity check for counts file existence
  if ( myfile ) {    
    std::cout << "  \n\n  Reading counts from a file.... \n\n";
    std::cout << "\n" << "Counts found in the " << cfile << " file:" << "\n";
    
    // specifies an integer and sets an indexing
    int count;
    int p = 0;

    // reads each count into the counts vector  
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
    std::cout << "  \n\n  Read counts successfully! \n\n\n  The analysis will start now.\n\n";

  } else {
    //prints if file not opened
    std::cout << "  Unable to open " << cfile << "."; 
    
  }
}

// reads the mapping vector
// mfile is the path to the mapping file
// mapping is a vector of integers that will carry the phi categories 
void read_mapping( std::string mfile, std::vector<int>& mapping ) {

  //reads the file from the terminal
  std::ifstream myfile (mfile); 

  //sanity check for mapping file existence
  if ( myfile ) {    
    std::cout << "  \n\n  Reading the mapping from file.... \n\n";
    std::cout << "\n" << "The mapping found in the " << mfile << " file:" << "\n";
    
    
    // specifies an integer and sets an indexing
    int numphi;

    while ( myfile >> numphi ) {
      mapping.push_back(numphi);
    }

    //closes the file
    myfile.close(); 

    if ( mapping.size() == 64) {
      // prints the counts
      for (int i=0; i<mapping.size(); i++){ 
        std::cout << "[" << i << "] " << mapping[i] << "\n";
      }
      std::cout << "\n"; // adds empty line at the end of output

      std::cout << "  \n\n  Read the mapping successfully! \n\n";
    } else {
    std::cerr << " Mapping has less than 64 values or some of the values are not integers (eg. NA). Please, specify the phi category for all codons. \n\n"; 
    exit(-1);
    }
  } else {
    //prints if file not opened
    std::cout << "  Unable to open " << mfile << "."; 
    
  }
}


////////////////////////
// statistical functions
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

