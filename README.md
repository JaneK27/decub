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

This step should have created an executable called `DECUB` in your working directory.

## How it works

## Control file

The control file (ending in `.ctl` in this repository) includen all the necessary parameters to run `DECUB`: this includes the location of the data, the parameters for the software and the MCMC.

* ``` countFile ```:
* ``` outputFile ```:
* ``` mapping ```:
* ``` MCMC_generations ```:
* ``` sample_frequency ```:

## Citation
_The patterns of codon usage between chordates and arthropods are different but co-evolving with mutational biases._ Kotari, Kosiol, and Borges (In preparation)

## License
This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.