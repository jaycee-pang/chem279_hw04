#pragma once 
#include <iostream> 
#include <vector> 
#include <stdexcept>
#include <fstream>
#include <functional> 
#include <algorithm> 
#include <map> 
#include <armadillo> 
using namespace arma;

double choose(int n, int k); 
int factorial(int n); 
int double_factorial(int n); 

class PrimitiveGaussian {
    private:

        std::vector<double> R_ ;  // center
        double alpha_; // (exponent)
        int l_; 
        double N_; // normalization coeff 
        std::vector<int> lmn_; 
        std::vector<double> coefficients_; 
        


    public: 
        PrimitiveGaussian();
        PrimitiveGaussian(double x, double y, double z, double alpha, int l, std::vector<int> lmn);
        PrimitiveGaussian(double x, double y, double z, int l, double alpha);
        double overlap1d(double xa, double xb, double alpha, double beta, int la, int lb) const; 
        int dim_func() const;
        // std::vector<std::vector<double>> overlap_Sab(const PrimitiveGaussian & other) const;
        double normalize_primitive();
        double overlap3d(const PrimitiveGaussian& other) const ; 
        std::vector<double> R() const;
        double alpha() const;
        double N() const;
        int l() const;
        std::vector<int> lmn() const; 
        void set_lmn(std::vector<int> lmn); 
        void set_alpha(double new_alpha);
        void set_l(int l);
        void set_R(std::vector<double> R);
        void set_N(double N); 

};
std::ostream& operator<<(std::ostream& os, const PrimitiveGaussian& p);
class ContractedGaussian {
    private: 
        // shell type? use enum for this too maybe 
        std::vector<double> R_; // xyz coords center
        int Z; // atomic number 
        std::vector<PrimitiveGaussian> primitives_; // vector of primitive Gaussians 
        std::vector<int> Ls_; 
        std::vector<int> lmn_;

        std::vector<double> ds_;
        std::vector<double> alphas_;  
        std::vector<double> Ns_; // norm constants for primitives 
        std::string shell_; 
    public: 
        ContractedGaussian();
        ContractedGaussian(std::vector<int> lmn, std::string shell);
        void addPrimitives(const std::vector<PrimitiveGaussian>& p, const std::vector<double>& contraction_coeffs);

        std::vector<double> getN();
        std::vector<std::vector<double>> normalized_overlap() const;

        std::vector<double> alphas() const;
        std::vector<int> Ls() const; // l+m+n = L;
        std::vector<double> R() const;
        std::vector<double> ds() const;
        std::vector<int> lmn() const; 
        void set_R(const std::vector<double>& center);
        std::string shell() const; 
        const std::vector<PrimitiveGaussian>& primitives() const;
        // std::vector<std::vector<double>> contracted_overlap(const ContractedGaussian & other) const;
        double contracted_overlap(const ContractedGaussian & other) const; 
        void print_overlap_matrix(const std::vector<std::vector<double>>& overlap_matrix) const;
}; 

std::ostream& operator<<(std::ostream& os, const ContractedGaussian& c); 


