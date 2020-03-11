/*
================================================================================
   ErrorMC
   A Monte Carlo Error Propagation
================================================================================

   *** Header file ***

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

//GSL Libraries
#include<gsl/gsl_rng.h> // GSL Random number generation
#include<gsl/gsl_randist.h>  // GSL Random distributions
#include<gsl/gsl_statistics_double.h> // for tests

#include "ErrorMC.cpp" //file with functions implementations

//=============================================================================
// Function Headers
//=============================================================================

//Parameter Generator
int GenParMC(std::vector<double>,std::vector<double>,
             gsl_rng*,std::vector<double>&);

//------------------------------------------------------------------------------
//Error Propagation Calculation
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
double ErrorCalc
(double,std::vector<double>,std::vector<double>,
 double(*)(double,std::vector<double>,void*),void*,int,double&);

#endif
