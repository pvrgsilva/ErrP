/*
================================================================================
   ErrorMC
   A Monte Carlo Error Propagation
================================================================================

   *** Header file ***

   Class definition

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


#ifndef ERRORMC_H
#define ERRORMC_H

#include<iostream>
#include<cstdio>
#include<vector>
#include<cmath>
#include<string>

//GSL Libraries
#include<gsl/gsl_rng.h> // GSL Random number generation
#include<gsl/gsl_randist.h>  // GSL Random distributions
#include<gsl/gsl_statistics_double.h> // for tests

//#include "ErrorMCdef.h" //file with functions implementations


typedef double (*funcmodel_t) (double,std::vector<double>,void*);

//#include "ErrorMC.cpp" //file with functions implementations

//=============================================================================
// Class declaration
//=============================================================================

class ErrorMC {
// private:
  std::vector<double> par_central;
  std::vector<double> par_sigma;
  // std::vector<double> result;
  void *extra_par;
  funcmodel_t model;
  int Nmax;
  unsigned long int seed;
  gsl_rng *r;
  const gsl_rng_type *rng_type;
  // string rand_method;
public:
  double average;

  ErrorMC(); //constructor (create others?)
  ~ErrorMC(); //destructor
  void SetModel();
  void SetN(int);
  void SetSeed(unsigned long int);
  // void SetDist(??);
  void SetRandMethod(std::string);
  void SetParCentral(std::vector<double>);
  void SetParSigma(std::vector<double>);
  // void SetParCov(); //to be implemented
  void SetExtraPar(void *);
  void SetModel(funcmodel_t);

//Parameter Generator
  // int GenParMC(std::vector<double>,std::vector<double>,gsl_rng*,std::vector<double>&);
  int GenParMC(std::vector<double>&);

  //Error Propagation Calculation
  //Inputs: independent variable of the model (x)
  //vectors with parameters central values and uncertainties,
  //function with model and extra parameters (both defined by user),
  //maximum number of points to be calculated (user's choice).
  //It creats random number generator and performs MC to calculate
  //Returns: error, and mean value (by reference)
    double ErrorCalc
    // (double,std::vector<double>,std::vector<double>,
    // double(*)(double,std::vector<double>,void*),
    // //funcmodel_t,
    // void*,int,double&);
    (double,
     // double(*)(double,std::vector<double>,void*),
     // funcmodel_t,
     // int,
     double&);
};

//=============================================================================
// Function Headers
//=============================================================================



//------------------------------------------------------------------------------


#endif
