/*
================================================================================
   ErrP
   A Class for (Numerical) Error Propagation
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


#ifndef ERRP_H
#define ERRP_H

#include<iostream>
#include<cstdio>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>

//GSL Libraries
#include<gsl/gsl_vector.h> // GSL Vectors
#include<gsl/gsl_matrix.h> // GSL Matrices
#include<gsl/gsl_linalg.h> // GSL Linear Algebra
#include<gsl/gsl_rng.h> // GSL Random number generation
#include<gsl/gsl_randist.h>  // GSL Random distributions
#include<gsl/gsl_statistics_double.h> // for tests

//#include "ErrPdef.h" //file with functions implementations


typedef double (*funcmodel_t) (double,std::vector<double>,void*);

//#include "ErrP.cpp" //file with functions implementations

//=============================================================================
// Class declaration
//=============================================================================

class ErrP {
// private:
  std::vector<double> par_central;
  std::vector<double> par_sigma;
  // std::vector<double> result;
  void *extra_par;
  funcmodel_t model;
  int Nmax;
  int Npar;
  unsigned long int seed;
  gsl_rng *r;
  gsl_matrix *covmatrix; //new
  gsl_matrix *covmatrix_fc; //after Cholesky factorization (used only for MC)
  const gsl_rng_type *rng_type;

  // string rand_method;
public:
  double average;

  ErrP(); //constructor (create others?)
  ErrP(std::string, unsigned long int); //constructor 2
  ErrP(std::string, unsigned long int, int); //constructor 3
  ~ErrP(); //destructor
  void SetModel();
  void SetN(int);
  void SetSeed(unsigned long int);
  // void SetDist(??);
  void SetRandMethod(std::string);
  void SetParCentral(std::vector<double>);
  void SetParSigma(std::vector<double>);
  void SetParCov(gsl_matrix *); //to be implemented
  void SetParCov(const char *);
  void SetExtraPar(void *);
  void SetModel(funcmodel_t);

//------------------------------- DERIVATIVES --------------------------------//
//Gradient of function model
//Inputs: independent variable of the model
//        and vector to store gradient components
  int CalcGrad(double,std::vector<double>&);

//------------------------------- MONTE CARLO --------------------------------//
//Parameter Generator (Monte Carlo)
  int GenParMC(std::vector<double>&);

//Parameter Generator (Monte Carlo) with Covariance
  int GenParMCCov(std::vector<double>&);

//Error Propagation Calculation Monte Carlo
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
  double ErrorCalcMC(double,double&);

//Error Propagation Calculation with Monte Carlo with Covariance
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
  double ErrorCalcCovMC(double,double&);

};
//------------------------------------------------------------------------------


#endif
