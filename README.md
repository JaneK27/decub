# DECUB

DECUB (disentangling the effects of codon usage bias) is a bayesian estimator that quantifies the mutational and selective biases that drive codon preferences in animals.

## Version

0.1 (December 2022)

## Downloading and compiling BESS

This GitHub repository consists of all the necessary files to compile and run DECUB. You can download these files manually or directly from the terminal by running this command: 

```

git clone git@github.com:JaneK27/decub.git decub

```

This will install a folder in your computer named decub.  DECUB is implemented in C++ and needs to be compiled. Compiling will produce an executable that can then be used to run DECUB. Here is an example of how to compile it using the `g++` compiler (here we used c++ version `v.11`)

```

g++ decub.cpp -o DECUB -O3 -std=c++11

```

This step should have created an executable called `DECUB` in your working directory. Otherwise you can download and use the already precompiled executable called `DECUB` found in this repository.

## Control file

The control file (ending in `.ctl` in this repository) contains all the necessary parameters to run `DECUB`: this includes the location of the data, the parameters for the software and the MCMC.

* ``` countFile ```: path to the file holding the count number of each codon found in the species of interest (more details below).
* ``` outputFile ```: the desired name of the output file (please include the path where you would like to save the file).
* ``` mapping ```: a number specifying the type of mapping you would like to use. Choose `1` for  an amino acid mapping (not accounting for variation in synonymous sites), `2` for an arthropods-specific mapping that accounts for synonymous sites, or `3` for a chordates-specific one. You can also choose `0` and input a file with a custom mapping. More information on what the mapping is and how it works, please see below.
* ``` MCMC_generations ```: the number of generations for which the MCMC is intended to run.
* ``` sample_frequency ```: the sample frequency determines how often the chain is sampled.

## How it works

DECUB (disentangling the effects of codon usage bias) is a software that uses a Bayesian estimator that can explain the observed variability of codon composition in light of mutational and selection biases. It allows the disentangling of mutation, selection on amino acids and synonymous codons, as well as GC-bias gene conversion (gBGC). 

### count file
DECUB needs as input a count file and a mapping. A count file is a .txt file that contains all the counts of each codon for a specific species. It should be space-delimited and follow the alphabetical sequence of codons (see below), including the stop codons. In the folder `example_files` you can also find an example of this file for _Drosophila melanogaster_ called `D.mel_count.txt`. 

| AAA 	| AAC 	| AAG 	| AAT 	| ACA 	| ACC 	| ACG   | ... 	| TGT 	| TTA 	| TTC 	| TTG 	| TTT 	|
|-----	|-----	|-----	|-----	|-----	|-----	|-----	|-----	|-----	|-----	|-----	|----- |-----  |

### mapping
The mapping refers to the process of assigning a codon category to a codon. It can be assume no preferences between codons (a single category for every codon) or allow for more complex scenarios. The maximum number of possible categories on DECUB is currently 52. Already in the software you can find a premade mapping that only accounts for non-synonymous variation (one category per amino acid + the three stop codons), and mappings made for arthropods and chordates (for more information on how these were calculated, please refer to the citation (_in preparation_)). Finally, you can input your own mapping following the same rules for the count file (following the alphabetical order of the codons, put a number that signifies the category of each codon). An example of a random mapping of 10 categories can be found in `example_files`.

After providing a count file and the appropriate mapping, DECUB uses the MCMC algorithm to produce the coefficients of mutational biases (mutations and gBGC combined) and of each codon category specified in the mapping. These describe the different frequences of codons found in each species based on the input count files. 

## Running DECUB

To run `DECUB`, specify the parameters in the control file and place it in the same location as the executable. Then, open the terminal and run the executable `DECUB` followed by the name of the control file:


```

./DECUB decub.ctl

```

DECUB will immediately print out a short description of the count file. If you have chosen to input a custom mapping, please include the path to the mapping file in the terminal too, like this:

```

./DECUB decub.ctl path_to_mapping/mapping.txt

```
In that case, DECUB will also print out the description of the mapping file. Please, confirm that this information is compliant to your data. To ensure DECUB is running, you should get the message `DECUB has started!`. When the analysis finishes, you should get the message `DECUB has finished!`. DECUB periodically writes to the `output_file`, based on the sample frequency you have specified in the control file.

## Output file

The output file includes information on the estimated parameters. Each line is a sample from the posterior. Here are the first lines of the output from running DECUB using the example count file with the example mapping provided in the `example_files` folder.

```
Gen lnL             beta0   beta1       beta2       beta3    phi0   phi1    phi2    phi3    phi4    phi5    phi6    phi7    phi8        phi9 
0   -3.68745e+07    1       0.135596    1.1065      0.633355 1      1.78291 2.06733 2.26921 2.83081 1.97883 2.15199 1.81068 2.26096     2.06783 
100 -3.08055e+07    1       0.944272    0.983288    0.859995 1      2.62513 1.83416 2.45937 1.71723 2.47216 1.85527 1.39234 1.16002     1.49529 
200 -3.07875e+07    1       1.00467     1.07403     0.859995 1      2.71318 1.82721 2.32623 1.62061 2.42503 1.6848  1.38745 1.18789     1.49959 
300 -3.0782e+07     1       1.02408     1.08828     0.860523 1      2.53877 1.70705 2.22708 1.46034 2.23427 1.57305 1.26659 1.11255     1.43925 
400 -3.07745e+07    1       1.03566     1.10036     0.895608 1      2.38643 1.56641 1.9772  1.33589 2.04071 1.42997 1.13632 0.996533    1.27376 
500 -3.07678e+07    1       1.04762     1.11367     0.895608 1      2.12041 1.42427 1.77603 1.13461 1.80089 1.25374 1.01949 0.904104    1.11836 


```

## Citation
_The patterns of codon usage between chordates and arthropods are different but co-evolving with mutational biases._ Kotari, Kosiol, and Borges (In preparation)

## License
This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.