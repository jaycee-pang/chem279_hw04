#pragma once 
#include <iostream> 
#include <vector> 
#include <stdexcept>
#include <fstream>
#include <functional> 
#include <algorithm> 
#include <map> 
#include <armadillo> 
#include "molecule.h"
arma::vec C_alphas; arma::vec H_alphas; arma::vec O_alphas;
arma::vec N_alphas; arma::vec F_alphas;
arma::mat H1s_coeffs; 
arma::mat C2s_coeffs; arma::mat C2p_coeffs; 
arma::mat O2s_coeffs; arma::mat O2p_coeffs;
arma::mat N2s_coeffs; arma::mat N2p_coeffs;
arma::mat F2s_coeffs; arma::mat F2p_coeffs; 
class AO {
    private: 
        std::string shell_type; // s, p, d, etc.
        arma::vec R_; // center 
        arma::vec alphas_; // exponens 
        arma::vec ds_; // contraction coefficients 
        // aram::vec Ns_; // don't need, normalize upon init in constructor 
        arma::vec lmn_;

    
    public: 
        AO(); 
        AO(std::string& shell, arma::vec& R, arma::vec& alphas, arma::vec& d_coeff, arma::uvec& lmn); 

        double overlap1d(double xa, double xb, double alpha, double beta, int la, int lb) const;
      
        // double overlap3d(const & other) const;
        
        void print_AO() const;

        std::string shell() const; 
        arma::vec R() const; 
        arma::vec alphas() const;
        arma::vec ds() const;
        arma::uvec lmn() const; 



}; 
double overlap3d(const arma::vec& Ra, const arma::vec& Rb, double alphas, double betas, const arma::uvec& la, const arma::uvec& lb); 

double evaluate_contracted_overlap(const AO& ao1, const AO& ao2); 


double choose(int n, int k); 
int double_factorial(int n); 

void getBasis_data(const std::string&atom_type); 
// std::ostream& operator<<(std::ostream& os, const PrimitiveGaussian& p) {
//     std::vector<double> center = p.R(); 
//     std::vector<int> lmn = p.lmn(); 
    
//     os << "Center: " << "(" << center[0] << ", " << center[1] << ", " << center[2] << ")\t" 
//     << "alpha: " << p.alpha() << "\t" << "L: " << p.l() <<"\t" 
//     << "shell: (" << lmn[0] << "," << lmn[1] << "," << lmn[2] << ")" << "\t" <<"N: " << p.N()<< std::endl;

    
//     return os;
// }


// std::ostream& operator<<(std::ostream& os, const ContractedGaussian& c) {
//     std::vector<PrimitiveGaussian> primitives = c.primitives();
//     std::vector<double> contraction_coeffs = c.ds(); 
//     os << "shell type: " << c.shell() << std::endl;
//     for (int i = 0; i < primitives.size(); i++) {
//         os << "Primitive " << i << ": " << std::endl;
//         os << "Primitives" << primitives[i];
//         os << "contraction coeff: " <<contraction_coeffs[i] << std::endl;
//     }

//     return os;
// }

// std::vector<std::vector<int>> get_lmn(int L) {
//     std::vector<std::vector<int>> combos; 
//     for (int l = L; l >= 0; l--) {
//         for (int m = L - l; m >= 0; m--) {
//             int n = L - l - m;
//             if (n >= 0) {
//                 combos.push_back({l, m, n});
//             }
//         }
//     }
//     return combos; 
// }
 

// std::ostream& operator<<(std::ostream& os, const Atom& atom) {
//     os << "Atomic Number: " << atom.Z <<std::endl;
//     os << "(" << atom.x << ", " << atom.y << ", " << atom.z << ")";
//     return os;
// }

// std::vector<double> PrimitiveGaussian::R() const {return R_;}
// double PrimitiveGaussian::alpha() const {return alpha_;}
// double PrimitiveGaussian::N() const {return N_;}
// int PrimitiveGaussian::l() const {return l_;}
// void PrimitiveGaussian::set_alpha(double new_alpha) {alpha_ = new_alpha;}
// void PrimitiveGaussian::set_l(int l) {l_ = l; }
// void PrimitiveGaussian::set_R(std::vector<double> R) {R_ = R;}
// void PrimitiveGaussian::set_lmn(std::vector<int> lmn) {lmn_ = lmn;}
// std::vector<int> PrimitiveGaussian::lmn() const {return lmn_;}
// void PrimitiveGaussian::set_N(double N) {N_ = N;}



// double get_N() {
//     double integral = overlap1d(R_[0], R_[0], alpha_, alpha_, lmn_[0], lmn_[0]) *
//                     overlap1d(R_[1], R_[1], alpha_, alpha_, lmn_[1], lmn_[1]) *
//                     overlap1d(R_[2], R_[2], alpha_, alpha_,lmn_[2], lmn_[2]);
    
//     N_ = 1.0 / std::sqrt(integral);
//     return N_;

// }
// ContractedGaussian::ContractedGaussian() {}
// ContractedGaussian::ContractedGaussian() {}
// void generateAOs(vector, int atomic_num) {
//     std::vector<PrimitiveGaussian> primitives; 
// }
// // ContractedGaussian::ContractedGaussian(std::vector<int> lmn, std::string shell): lmn_(lmn), shell_(shell){}
// // void ContractedGaussian::addPrimitives(const std::vector<PrimitiveGaussian>& p, const std::vector<double>& contraction_coeffs) {
// //     primitives_ = p;
// //     ds_ = contraction_coeffs;

// // }

// std::vector<double> ContractedGaussian::getN() {
//     std::vector<double> norm_constants; 
//     for (PrimitiveGaussian & p: primitives_) {
//         double N = p.normalize_primitive();
//         p.set_N(N);
//         norm_constants.push_back(N);
//     }
//     Ns_ = norm_constants; 
//     return norm_constants; 
    
// }

// double ContractedGaussian::contracted_overlap(const ContractedGaussian & other) const {
//     double overlap = 0.0;
//     for (int i = 0; i < 3; i++) { 
//         for (int j = i; j < 3; j++) {
//             double primitive_overlap = primitives_[i].overlap3d(other.primitives_[j]);
//             overlap += Ns_[i] * other.Ns_[j]*ds_[i] * other.ds_[j] * primitive_overlap;
//         }
//     }
//     return overlap; 

// }


// std::vector<double> ContractedGaussian::alphas() const {return alphas_;}
// std::vector<int> ContractedGaussian::Ls() const {return Ls_;} // l+m+n = L;
// std::vector<double> ContractedGaussian::R() const {return R_;}
// std::vector<double> ContractedGaussian::ds() const {return ds_;}
// std::string ContractedGaussian::shell() const {return shell_;} 
// void ContractedGaussian::set_R(const std::vector<double>& center) {R_ = center;}
// const std::vector<PrimitiveGaussian>& ContractedGaussian::primitives() const {
//         return primitives_;
// }

