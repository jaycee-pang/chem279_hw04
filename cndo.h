#pragma #once 
#include <iostream>
#include <armadillo> 
#include <map> 
#include <vector> 
#include <cmath> 
/*
General algorithm
1. Guess initial P 
2. Build Fock matrix with guessed P 
        - need function to input a guessed/temp P 
        - have Fock matrix otherwise set up 
            ie. calculate gamma, have I and A
3. Solve eigenvalue problems to update MO coefficients and eigenvalues 
    (for alpha and beta): FC = CE 
4. Make new density matrices P = CC.T 
    occupy p and q lowest energy MOs for alpha and beta 
5. Compare new density matrices, converged to some threshold? 
    arma::aprox_equal(Pa, Pa_old, "absdiff", tol)
6. If converged, calculate total energy. Else, return to step 2 
*/

#include "molecule.h"

extern std::map <std::string, double> semi_emp_s;
extern std::map <std::string, double> semi_emp_p;
extern std::map <std::string, double> semi_empirical; 
extern std::map <std::string, int> atomic_bonding;
extern const double au_to_EV; // 27.211 eV/au
class CNDO {
    private:
        // Fock matrices alpha and beta 
        arma::mat Fa_; arma::mat Fb_; 
        Molecule molecule; 
        int natoms_; int n_electrons_; int N_; // # basis functions 
        // alpha and beta density matrix P 
        arma::mat Pa_; arma::mat Pb_; 
        // total density matrix 
        arma::mat P_; 

        arma::mat Gamma_; // gamma ab matrix gamma_AB is NA xNA (# atoms)
        arma::mat Ga_; arma::mat Gb_;
        // mol coeffs alpha and beta 
        arma::mat Ca_, Cb_; // C beta matrix  // C alpha matrix
        // eigenvalues 
        arma::vec Ea_; arma::vec Eb_;
        arma::mat H_; // core Hamiltonian
        arma::mat S_; // overlap matrix 

        int q_, p_; // num alpha and beta electrons
        // params for the SCF algorithm 
        double tol; int max_it; 
        double Ec_, Eel_, Etot_; // Nuclear repulsion, electronic, total energy 

        
        

    public: 
        CNDO(Molecule& molecule, double tolerance, int max_it);
        int find_ooc(); 
        void build_F(); 

        void buildH(); 
        // maybe move G calculations to molecule class since it only has to do with atoms 
        void build_G(); 
        void update_P(); 
        void G_ab(); 
        void calc_Ec(); 
        double calculate_E(); 
        void SCF(); 

        double gamma_ab(const AO &ao1, const AO& ao2); 

        double integral6d(double sigma_A, double sigma_B, arma::vec Ra, arma::vec Rb);
        void print_CNDO_info() const; 
        void print_SCF() const; 
        // void update_F(); 

};


