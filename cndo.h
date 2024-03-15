#pragma #once 
#include <iostream>
#include <armadillo> 
#include <map> 
#include <vector> 
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
#include 
class CNDO {
    private:
        arma::mat Fa_; // fock matrix alpha

        Molecule molecule; 
        arma::mat Pa_; // alpha density matrix 
        arma::mat Pb_; // beta density matrix 
        arma::mat P_; // density matrix --> alpha or beta should this class be interactive or for both 
        arma::mat Gamma_; // gamma ab matrix
        arma::mat Ca_, Cb_; // C beta matrix  // C alpha matrix
 
        arma::mat Fb_; // beta Fock matrix 
        int q, p; // num alpha and beta electrons

    public: 
        CNDO()
        void build_F(); 
        void build_G() {

        }
        

};
void SCF(); // driver code 
