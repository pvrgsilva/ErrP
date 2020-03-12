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

// //GSL Libraries
#include<gsl/gsl_rng.h> // GSL Random number generation
#include<gsl/gsl_randist.h>  // GSL Random distributions

// #include "ErrorMCdef.h" //header file
#include "ErrorMC.h" //header file

//==============================================================================
// Function Implementations
//==============================================================================

//constructor
ErrorMC::ErrorMC(){

  //creating random number generator
  //(for now, using default parameters)

  unsigned long int seed = 0; //seed

  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc(T);
  gsl_rng_set(r,seed);

}

void ErrorMC:SetParCentral(std::vector<double> user_info){
  par_central = user_par;
}

void ErrorMC:SetParSigma(std::vector<double> user_info){
  par_sigma = user_info;
}

void ErrorMC:SetExtraPar(void *user_info){
  extra_par = user_info;
}

//Parameter Generator
//Gaussian distribution
//Inputs: vectors with parameters central values and uncertainties,
// pointer to random number generator (GSL) and vector to store result
int ErrorMC::GenParMC
(std::vector<double> par_central, std::vector<double> par_sigma,
 gsl_rng *r ,std::vector<double> &result)
{
  int npar = par_central.size();
  int npar_sigma = par_sigma.size();

  if(npar != npar_sigma){
    std::cout << "Parameters and Sigma must have same number of entries.\n";
    return -1;
  }


  for(int i=0; i<npar;i++){
    result.push_back(gsl_ran_gaussian(r,par_sigma.at(i)) + par_central.at(i));
  }

  return 0;

}
//end of function GenParMC

//------------------------------------------------------------------------------
//Error Propagation Calculation
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
double ErrorMC::ErrorCalc
(double x, std::vector<double> par_central, std::vector<double> par_sigma,
 double(*model)(double,std::vector<double>,void*),
 // funcmodel_t model,
 void *extra_par, int Nmax, double &average)
{

  // //creating random number generator
  // //(for now, using default parameters)
  //
  // unsigned long int seed = 0; //seed
  //
  // const gsl_rng_type *T = gsl_rng_default;
  // gsl_rng *r = gsl_rng_alloc(T);
  // gsl_rng_set(r,seed);

  std::vector<double> Ymc, parmc;
  double y; int status;

// loop for calculations
  for(int i=0;i<Nmax;++i){
    status = GenParMC(par_central,par_sigma,r,parmc);
    if(status==-1) continue; //-1 means error, and should skip lines below
    y = model(x,parmc,extra_par);//
    Ymc.push_back(y);
    parmc.clear();
  }

  //mean value and standard deviation

  double sdev = 0;
  double sum1 = 0;
  double sum2 = 0;

  int Ndata = Ymc.size();

//mean

  for(int i=0;i<Ndata;++i){
    y = Ymc.at(i);
    sum1 = sum1+y;
    // sum2 = sum2 + y*y;
  }

  average = sum1/Ndata;

  // const double *aux = &Ymc[0];
  // average = gsl_stats_mean(aux,1,Ndata);

// standard deviation


  double Delta = 0;//abs(sum2 - sum1*sum1);

  for(int i=0;i<Ndata;++i){
    Delta = Ymc.at(i)-average;
    sum2 = sum2 + Delta*Delta;
  }

  sdev = sqrt(sum2/(Ndata-1.0)); // problem?

  // sdev = gsl_stats_sd_m(aux,1,Ndata,average);

  gsl_rng_free(r); // free memory associated to random number generator
  Ymc.clear();
  parmc.clear();

  return sdev;
}
//end of function ErrorCalc
