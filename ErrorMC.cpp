/*
================================================================================
   ErrorMC
   A Monte Carlo Error Propagation
================================================================================

   *** File with functions implementations ***

   Important notes and descriptions:
   ....


   Dependencies:

   * Gnu Scientific Library (GSL) for
     - Random Number Generation
     - Random Distributions


  Current version: v1
  Created by PVRG Silva, Pelotas, 11/03/2020
  Last update by PVRG Silva, Pelotas, 11/03/2020

--------------------------------------------------------------------------------
  Code by Paulo V.R.G. Silva
  pvrgsilva@ufpel.edu.br | pvrecchia@gmail.com
  High and Medium Energy Group
  Grupo de Altas e Medias Energias (GAME)
  Universidade Federal de Pelotas, Pelotas - RS (Brasil)
--------------------------------------------------------------------------------

*/

#include<iostream>
#include<cstdio>
#include<vector>
#include<cmath>

//GSL Libraries
#include<gsl/gsl_rng.h> // GSL Random number generation
#include<gsl/gsl_randist.h>  // GSL Random distributions

#include<ErrorMC.h> //header file

//==============================================================================
// Function Implementations
//==============================================================================

//Parameter Generator
//Gaussian distribution
//Inputs: vectors with parameters central values and uncertainties,
// pointer to random number generator (GSL) and vector to store result
int GenParMC
(vector<double> par_central, vector<double> par_sigma,
gsl_rng *r ,vector<double> &result)
}
  int npar = par_central.size();
  int npar_sigma = par_sigma.size();

  if(npar != npar_sigma){
    cout << "Parameters and Sigma must have same number of entries.\n";
    return -1;
  }


  for(int i=0; i<npar;i++){
    result.push_back(gsl_ran_gaussian(rand,par_sigma.at(i)) + par_central.at(i));
  }

  return 0;

}
//end of function GenParMC

//------------------------------------------------------------------------------
//Error Propagation
