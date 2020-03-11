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

//=============================================================================
// Function Headers
//=============================================================================

//Parameter Generator
int GenParMC(vector<double>, vector<double>,gsl_rng*,vector<double>&);


#endif
