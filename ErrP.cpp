/*
================================================================================
   ErrP
   A Class for (Numerical) Error Propagation
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

#include "ErrP.hpp" //header file

//==============================================================================
// Function Implementations
//==============================================================================

//default constructor (default parameters)
ErrP::ErrP(){

  statepar=false;
  statesig=false;
  statecov=false;
  statefunc=false;

  //creating random number generator

  seed = 0; //seed

  rng_type = gsl_rng_default;

 //Tried move to GenParMC, but not worked
  r = gsl_rng_alloc(rng_type);
  gsl_rng_set(r,seed);

}


//Constructor (w/ parameters passed by the user)
ErrP::ErrP(std::string user_method, unsigned long int user_seed){

  statepar=false;
  statesig=false;
  statecov=false;
  statefunc=false;

  //creating random number generator

  seed = user_seed; //seed

  SetRandMethod(user_method);

  // rng_type = gsl_rng_default;

 //Tried move to GenParMC, but not worked
  r = gsl_rng_alloc(rng_type);
  gsl_rng_set(r,seed);

}

//Constructor (w/ parameters passed by the user)
ErrP::ErrP(std::string user_method, unsigned long int user_seed, int Nuser){

  statepar=false;
  statesig=false;
  statecov=false;
  statefunc=false;

  //creating random number generator

  seed = user_seed; //seed
  Nmax = Nuser;

  SetRandMethod(user_method);

  // rng_type = gsl_rng_default;

 //Tried move to GenParMC, but not worked
  r = gsl_rng_alloc(rng_type);
  gsl_rng_set(r,seed);

}

//destructor
ErrP::~ErrP(){
  // gsl_rng_free(r); // free memory associated to random number generator
  std::cout << "Releasing memory...\n";
}

void ErrP::SetParCentral(std::vector<double> user_info){
  par_central = user_info;
  Npar = par_central.size();

  statepar=true;

  // int Npar = par_central.size();
  // std::cout << "Testing passing parameters\n";
  // for(int i = 0;i<Npar;i++){
  //   std::cout << par_central[i] << "\n";
  // } //working
}

void ErrP::SetParSigma(std::vector<double> user_info){
  par_sigma = user_info;

  statesig=true;
  // int Npar = par_sigma.size();
  // std::cout << "Testing passing parameter errors\n";
  // for(int i = 0;i<Npar;i++){
  //   std::cout << par_sigma[i] << "\n";
  // } //working

}

int ErrP::SetParCov(gsl_matrix *user_matrix){

  // if(Npar == 0){
  if(statepar==false){
    std::cout << "\nPlease, set parameter values first\n";
    return -1;
  }else{
    covmatrix = gsl_matrix_alloc(Npar,Npar);
    gsl_matrix_memcpy(covmatrix, user_matrix);

    //perform Cholesky decomposition (needed for error propagation)
    covmatrix_fc = gsl_matrix_alloc(Npar,Npar);
    gsl_matrix_memcpy(covmatrix_fc, covmatrix);
    gsl_linalg_cholesky_decomp1(covmatrix_fc);
    statecov==true;
    return 0;
  }

}


int ErrP::SetParCov(const char *file_path){

  // if(Npar == 0){
  if(statepar==false){
    std::cout << "\nPlease, set parameter values first\n";
    return -1;
  }else{
    covmatrix = gsl_matrix_alloc(Npar,Npar);
    FILE *mfile = fopen(file_path,"r");
    int status = gsl_matrix_fscanf(mfile,covmatrix);

    if(status==0){
      std::cout << "\nMatrix elements read with success from file.\n";
    }
    fclose(mfile);

    //perform Cholesky decomposition (needed for error propagation)
    covmatrix_fc = gsl_matrix_alloc(Npar,Npar);
    gsl_matrix_memcpy(covmatrix_fc, covmatrix);
    gsl_linalg_cholesky_decomp1(covmatrix_fc);
    statecov=true;
    return 0;
  }



}


void ErrP::SetExtraPar(void *user_info){

  // std::cout << "Testing Extra parameter copy\n";

  extra_par = user_info;

  // std::cout << "Copied?\n";
  // std::cout << "Size: " << sizeof(extra_par) << "\n";
  // std::cout << "Info " << *(float*)extra_par << "\n";
  // std::cout << "Yes!!\n";
  //working!!!
}

void ErrP::SetN(int Nuser){
  Nmax = Nuser;
}

void ErrP::SetSeed(unsigned long int seed_user){
  seed = seed_user;
}

void ErrP::SetRandMethod(std::string user_method){

  if(user_method.compare("default")==0){
    rng_type = gsl_rng_default;
  }
  if(user_method.compare("mt19937")==0){
    rng_type = gsl_rng_mt19937;
  }
  if(user_method.compare("taus")==0){
    rng_type = gsl_rng_taus;
  }
  if(user_method.compare("taus2")==0){
    rng_type = gsl_rng_taus2;
  }
  if(user_method.compare("ranlux")==0){
    rng_type = gsl_rng_ranlxs0;
  }
  if(user_method.compare("ranlux1")==0){
    rng_type = gsl_rng_ranlxs1;
  }
  if(user_method.compare("ranlux2")==0){
    rng_type = gsl_rng_ranlxs2;
  }
  //add more ...

}

void ErrP::SetModel(funcmodel_t user_func){
  model = user_func;
  statefunc = true;
}


void ErrP::UseCovariance(bool info_user){
  statecov = info_user;

  std::cout << "Use cov matrix is set to " << statecov << "\n";
}

//------------------------------- DERIVATIVES --------------------------------//
//Gradient of function model
//Inputs: independent variable of the model
//        and vector to store gradient components
int ErrP::CalcGrad(double x,std::vector<double> &grad){


  // if(statefunc==false){
    // std::cout << "Function model not informed\n";
    // return -1;
  // }else{
  double eps = 10e-8;

  // std::vector<double> par_aux = par_central;
  std::vector<double> par_deriv;
  par_deriv = par_central;
  double deriv,delta;

  for(int i = 0; i<Npar; i++){
    par_deriv.at(i) = par_central.at(i) + eps;
    delta = model(x,par_deriv,extra_par) - model(x,par_central,extra_par);
    deriv = delta/eps;
    grad.push_back(deriv);
    par_deriv.at(i) = par_central.at(i);
  }

  par_deriv.clear();
// }

  return 0;

}

//Needs new tests!!
double ErrP::ErrorCalcGrad(double x){

  gsl_matrix *cov_aux = gsl_matrix_alloc(Npar,Npar);

  if(statefunc==false){
    std::cout << "Function model not informed\n";
    return -1;
  }else if(statepar==false){
    std::cout << "Parameters not informed\n";
    return -1;
  }else if(statesig==false && statecov==false){
    std::cout << "Uncertainties and covariance matrix not informed\n";
    return -1;
  }else{

    if(statesig==true && statecov==false){
    std::cout << "Error propagation without covariance\n";
    // gsl_matrix_set_zero(cov_aux);
    double error_par=0;
    for(int i=0;i<Npar;i++){
      error_par = par_sigma.at(i);
      for(int j=0;j<Npar;j++){
        if(i==j) {gsl_matrix_set(cov_aux,i,j,error_par*error_par);}
        else {gsl_matrix_set(cov_aux,i,j,0);}
      }
    }
  }else// if(statesig==false && statecov==true){
    {

    std::cout << "Error propagation with covariance\n";
    gsl_matrix_memcpy(cov_aux,covmatrix);
  }
    std::vector<double> grad_aux;
    int state = CalcGrad(x,grad_aux);

    if(state!=0){
      std::cout << "Gradient not calculated\n";
      return 0;
    }else{
      gsl_vector *grad = gsl_vector_alloc(Npar);
      for(int i=0;i<Npar;i++){
        gsl_vector_set(grad,i,grad_aux.at(i));
      }
      gsl_vector *prodSigmaGrad = gsl_vector_alloc(Npar);

//do product CovMatrix*Grad(model) = vec(v)
      // gsl_blas_dsymv(CblasUpper,1,cov_aux,grad,0,prodSigmaGrad);
      gsl_blas_dgemv(CblasNoTrans,1,cov_aux,grad,0,prodSigmaGrad);

      double result=0;

//do product Transpose(Grad(model))*vec(v)
      gsl_blas_ddot(grad,prodSigmaGrad,&result);

      double error = sqrt(result);

      grad_aux.clear();
      gsl_vector_free(grad);
      gsl_vector_free(prodSigmaGrad);
      gsl_matrix_free(cov_aux);
      return error;
    }
}


}



//------------------------------- MONTE CARLO --------------------------------//

//Parameter Generator (Monte Carlo)
//Gaussian distribution
//Returns parameters in vector passed as argument
int ErrP::GenParMC(std::vector<double> &result)
{
  // int npar = par_central.size();
  int npar_sigma = par_sigma.size();

  if(Npar != npar_sigma){
    std::cout << "Parameters and Sigma must have same number of entries\n";
    return -1;
  }else if(statesig==false){
    std::cout << "Sigma not informed\n";
    return -1;
  }//else if(statepar==false){
  //   std::cout << "Parameters not informed\n";
  //   return -1;
  // }else if(statefunc==false){
  //   std::cout << "Function model not informed\n";
  //   return -1;
  // }
  else{

  for(int i=0; i<Npar;i++){
    if(par_sigma.at(i)!=0){
      result.push_back(gsl_ran_gaussian(r,par_sigma.at(i)) + par_central.at(i));
    }else{
      result.push_back(par_central.at(i)); // if error = 0 (fixed parameter), do not change
    }

  }

  return 0;
}

}
//end of function GenParMC



//Parameter Generator (Monte Carlo) with Covariance
//Gaussian distribution
//Returns parameters in vector passed as argument
int ErrP::GenParMCCov(std::vector<double> &result)
{
  // int npar = par_central.size();
  int npar_sigma = par_sigma.size();

  if(Npar != npar_sigma){
    std::cout << "Parameters and Sigma must have same number of entries.\n";
    return -1;
  }else if(statecov==false){
    std::cout << "Covariance Matrix not informed\n";
    return -1;
  }else{


  gsl_vector *res_aux = gsl_vector_alloc(Npar);
  gsl_vector *mu = gsl_vector_alloc(Npar);

  for(int i=0; i<Npar;i++){
    gsl_vector_set(mu,i,par_central.at(i));
  }

  gsl_ran_multivariate_gaussian(r,mu,covmatrix_fc,res_aux);

  for(int i=0; i<Npar;i++){
    result.push_back(gsl_vector_get(res_aux,i));
  }

  gsl_vector_free(res_aux);
  gsl_vector_free(mu);

  return 0;
}

}
//end of function GenParMCCov


//------------------------------------------------------------------------------
//Error Propagation Calculation with Monte Carlo
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
// double ErrP::ErrorCalc
// (double x, std::vector<double> par_central, std::vector<double> par_sigma,
//  double(*model)(double,std::vector<double>,void*),
//  // funcmodel_t model,
//  void *extra_par, int Nmax, double &average)
double ErrP::ErrorCalcMC(double x,double &average)
{
  std::vector<double> Ymc, parmc;
  double y; int status;

  if(statefunc==false){
    std::cout << "Function model not informed\n";
    return -1;
  }else if(statepar==false){
    std::cout << "Parameters not informed\n";
    return -1;
  }else{
    if(statecov==false){
      std::cout << "Calculating without covariance matrix\n";
// loop for calculations
      for(int i=0;i<Nmax;++i){
        status = GenParMC(parmc);//(par_central,par_sigma,r,parmc);
        if(status!=0) continue; //-1 means error, and should skip lines below
        y = model(x,parmc,extra_par);//
        Ymc.push_back(y);
        parmc.clear();
      }
    }
    if(statecov==true){
      std::cout << "Calculating with covariance matrix\n";
  // loop for calculations
      for(int i=0;i<Nmax;++i){
        status = GenParMCCov(parmc);//(par_central,par_sigma,r,parmc);
        if(status!=0) continue; //-1 means error, and should skip lines below
        y = model(x,parmc,extra_par);//
        Ymc.push_back(y);
        parmc.clear();
      }
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

    sdev = sqrt(sum2/(Ndata-1.0));

  // sdev = gsl_stats_sd_m(aux,1,Ndata,average);

  // gsl_rng_free(r); // free memory associated to random number generator
    Ymc.clear();
    parmc.clear();

    return sdev;
}
}
//end of function ErrorCalcMC

//------------------------------------------------------------------------------
//ATTENTION: With boolean flags, this function is not necessary.
//Error Propagation Calculation with Monte Carlo with Covariance
//****Working. Error still larger than from analytical calculations,
// but still smaller than without covariance. Check method!!
//
//Inputs: independent variable of the model (x)
//vectors with parameters central values and uncertainties,
//function with model and extra parameters (both defined by user),
//maximum number of points to be calculated (user's choice).
//It creats random number generator and performs MC to calculate
//Returns: error, and mean value (by reference)
// double ErrP::ErrorCalc
// (double x, std::vector<double> par_central, std::vector<double> par_sigma,
//  double(*model)(double,std::vector<double>,void*),
//  // funcmodel_t model,
//  void *extra_par, int Nmax, double &average)
double ErrP::ErrorCalcCovMC(double x,double &average)
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
    status = GenParMCCov(parmc);//(par_central,par_sigma,r,parmc);
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

  sdev = sqrt(sum2/(Ndata-1.0));

  // sdev = gsl_stats_sd_m(aux,1,Ndata,average);

  // gsl_rng_free(r); // free memory associated to random number generator
  Ymc.clear();
  parmc.clear();

  return sdev;
}
//end of function ErrorCalcCovMC
