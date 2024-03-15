// #ifndef AO_H
// #define AO_H
#pragma once
#include <iostream> 
#include <vector> 
#include <stdexcept>
#include <fstream>
#include <functional> 
#include <algorithm> 
#include <map> 
#include <armadillo> 
// #include "molecule.h"
// class AO;

class AO {
    private: 
        std::string shell_type; // s, p, d, etc.
        arma::vec R_; // center 
        arma::vec alphas_; // exponens 
        arma::vec ds_; // contraction coefficients 
        // aram::vec Ns_; // don't need, normalize upon init in constructor 
        arma::uvec lmn_;

    
    public: 
        // AO(); 
        // AO(const std::string& shell, arma::vec& R, arma::vec& alphas, arma::vec& d_coeff, arma::uvec& lmn); 
        AO(const std::string& shell, const arma::vec& R, const arma::vec& alphas, const arma::vec& d_coeff, const arma::uvec& lmn);

        
      
        // double overlap3d(const & other) const;
        
        void print_AO() const;

        std::string shell() const; 
        arma::vec R() const; 
        arma::vec alphas() const;
        arma::vec ds() const;
        arma::uvec lmn() const; 



}; 
// #endif
double overlap1d(double xa, double xb, double alpha, double beta, int la, int lb);
double overlap3d(const arma::vec& Ra, const arma::vec& Rb, double alphas, double betas, const arma::uvec& la, const arma::uvec& lb); 

double evaluate_contracted_overlap(const AO& ao1, const AO& ao2); 


double choose(int n, int k); 
int double_factorial(int n); 
