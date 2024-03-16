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

std::map <std::string, double> semi_emp_s = {{"H1s", 7.716}, {"C2s", 14.051}, {"N2s", 19.316},
                                            {"O2s", 25.390}, {"F2s", 32.272}};
std::map <std::string, double> semi_emp_p = {{"C2p", 5.572}, {"N2p", 7.275},
                                            {"O2p", 9.111}, {"F2p", 11.080}};
std::map <std::string, double> semi_empirical = {{"H1s", 7.716}, {"C2s", 14.051}, {"N2s", 19.316},
                                            {"O2s", 25.390}, {"F2s", 32.272}, {"C2p", 5.572}, {"N2p", 7.275},
                                            {"O2p", 9.111}, {"F2p", 11.080}}; 
std::map <std::string, int> atomic_bonding = {{"H", 9}, {"C", 21}, {"N", 25}, {"O", 31}, {"F", 39} };
const double au_to_EV = 27.211; // 27.211 eV/au
class CNDO {
    private:
        arma::mat Fa_; // fock matrix alpha
        Molecule molecule; 
        arma::mat Pa_; // alpha density matrix 
        arma::mat Pb_; // beta density matrix 
        arma::mat P_; // density matrix --> alpha or beta should this class be interactive or for both 
        arma::mat Gamma_; // gamma ab matrix gamma_AB is NA xNA (# atoms)
        arma::mat Ca_, Cb_; // C beta matrix  // C alpha matrix
        arma::mat H_; // core Hamiltonian
        arma::mat Fa_; arma::mat Fb_; // beta Fock matrix 
        int q, p; // num alpha and beta electrons
        double tol; 
        int max_it; 
        int N_; // # basis functions 
        int natoms_; 
        int n_electrons;
        arma::mat S_; // overlap matrix 

    public: 
        CNDO(Molecule& molecule, double tolerance, int max_it) : molecule(molecule), tol(tolerance), max_it(max_it) {
            N_ = molecule.N();
            natoms_ = molecule.natoms(); 
            n_electrons = molecule.n_electrons();
            // make sure S is made, To-Do: bool if S is made in molecule class 
            molecule.make_overlap_matrix();
            S_ = molecule.S_overlap(); 
            // p = n_electrons / 2 + m; 
            // q = n_electrons /2 -m; 
            Gamma_.resize(natoms_, natoms_); // G only relies on atoms A and B 


        }
        void build_F() {
            int AO_mu = 0; 
            for (int i = 0; i < natoms_; i++) {
                const Atom& atomA = molecule.get_atom(i);
                const std::vector<AO> AAOs = molecule.atom_AOs(atomA.element); 
                int ZA = atomA.valence_e; 
                // for same atom AA; 
                for (const AO& ao_A: AAOs) {
                    std::string shell_type = ao_A.shell();
                    double IA = semi_empirical[shell_type]; 
                    double betaA = atomic_bonding[atomA.element];
                    // density elements for alpha and beta for A 
                    double pAA_alpha = Pa_(AO_mu, AO_mu); 
                    double pAA_beta = Pb(AO_mu, AO_mu);
                    double G_AA = Gamma_(i,i); // like gamma_AA 
                    
                    Fa_(AO_mu, AO_mu) = -IA + ((pAA_alpha - ZA) - (pAA_alpha - 0.5))  * G_AA; // sum (pBB-ZB)*GAB;
                    Fb_(AO_mu, AO_mu) = -IA + ((pAA_alpha - ZA) - (pAA_beta - 0.5))  * G_AA;

                    // for B!=A 
                    for (int j = 0; j < natoms_; j++) {
                        
                        if (i!=j) {
                            const Atom& atomB = molecule.get_atom(j);
                            const std::vector<AO> BAOs = molecule.atom_AOs(atomB.element); 
                            double G_AB = G_(i,j); 
                            double betaB = atomic_bonding[atomB.element];
                            double betaAB = 0.5 * (betaA + betaB);

                            for (const AO& ao_B : BAOs) {
                                for (int AO_nu = 0; AO_nu < N_; AO_nu++) {
                                    // Adjust condition if necessary
                                    if (AO_mu != AO_nu) {
                                        Fa_(AO_mu, AO_nu) += betaAB * S(AO_mu, AO_nu) - Pa_(AO_mu, AO_nu) * G_AB; 
                                        Fb_(AO_mu, AO_nu) += betaAB * S(AO_mu, AO_nu) - Pb_(AO_mu, AO_nu) * G_AB;
                                    }
                                }
                            }
                        }

                        
                    }
                    AO_mu++;
                }

        }
        // core Hamiltonian
        void buildH() {
            // need ZA atomic numbers (valence e numbers)
            // beta_AB*S_munu
            H_.resize(N_, N_);
            int AO_mu = 0; 
            for (int mu = 0; mu < natoms_; mu++) {
                const Atom& atomA = molecule.get_atom(mu); 
                int ZA = atomA.valence_e; 
                double G_AA = Gamma_(mu,mu); 
                const std::vector<AO> AAOs = molecule.atom_AOs(atomA.element); 
                for (const AO& ao_A : AAOs) {
                    std::string shell_type = ao_A.shell();
                    double IA = semi_empirical[shell_type]; 
                    double ZB_GAB = 0.0; // sum 
                    int AO_nu = 0; 
                    for (int nu = 0; nu < natoms_; nu++) {
                        const Atom& atomB = molecule.get_atom(nu); 
                        const std::vector<AO> BAOs = molecule.atom_AOs(atomB.element); 
                        int ZB = atomA.valence_e; 
                        double betaA = atomic_bonding[atomA.element];
                        double betaB = atomic_bonding[atomB.element]; 
                        double betaAB = 0.5 * (betaA + betaB);
                        if (mu!=nu) {
                            // if atomB != atomA
                            ZB_GAB += ZB*Gamma_(mu,nu);
                        }
                        for (const AO& ao_B :BAOs) {
                            if (AO_mu != AO_nu) {
                                H_(AO_mu, AO_nu) = betaAB * S_(AO_mu, AO_nu); 
                            }
                            AO_nu ++; 
                        }
                    } 
                    H_(AO_mu, AO_mu) = -IA - (ZA - 0.5) * G_AA - ZB_GAB;
                    AO_mu++; 
                }
            }

        }
        // maybe move G calculations to molecule class since it only has to do with atoms 
     
        void build_G() {
            // use elements of the Pa Pb total density matrix to update
            /* 
            for i < num atoms 
                for j < atoms 
                    get AO i, get AO j
                    evlauate the 2 electron integral s shell electrons 
            */
            Gamma_.zeros(natoms_, natoms_);
            for (int i = 0; i < natoms_; i ++) {
                for (int j = 0; j < natoms_; j++) {
                    Atom& atom_i = molecule.get_atom(i);
                    Atom& atom_j = molecule.get_atom(j);
                    AO& ao_i = get_AO(atom_i.element + "s");
                    AO& ao_j = get_AO(atom_j.element+ "s");
                    Gamma_(i,j) = gamma_ab(ao_i, ao_j);

                }
           }

        }
        void calculate_V() {

        }

        // driver code 
        void SCF() {
            // Fa, Fb 
            // build core hamiltonian
            for (int i = 0; i < max_iter; i++) {
                // build Fock matrix with iniitla guess P 
                Fa_; 
                Fb_;
                // solve eigenvalue problesmm to get MO coeff and eigenvalues FC=Ceps
                // make new density matrices 
                //      Pa_ vs Pa_old, Pb_ vs Pb_old
                // check for convergence 
                if (arma::aprox_equal(Pa_, Pa, "absdiff", tol))

            }
            

            

        } 
        

};

double gamma_ab(const AO &ao1, const AO& ao2) {
    arma::uvec lmn1 = ao1.lmn(); arma::uvec lmn2 = ao2.lmn(); 
    arma::uvec lmns = {0, 0, 0}; // s-shell only 
    if (lmn1!= lmns || lmn2!= lmns) {
        throw std::invalid_argument("The two-electron integral evaluation occurs between s-shell only."); 
    }
    
    arma::vec da(ao1.ds().size); arma::vec db(ao2.ds().size); 
    arma::vec a_alphas = ao1.alphas(); arma::vec b_alphas = ao2.alphas();
    assert(a_alphas.size() == b_alphas.size())
    arma::vec Ra = ao1.R(); arma::vec Rb = ao2.R();
    double sum = 0.0; 
    // iterate over k, k', l, l' to get all primitives
    for (size_t k = 0; i < a_alphas.size(); k++) {
        for (size_t kp = 0; kp < a_alphas.size(); kp++) {
            double sigma_A = 1.0/(a_alphas(k) + a_alphas(kp));
            for (size_t l=0; l < b_alphas.size(); l++) {
                for (size_t lp =0; lp < b_alphas.size(); lp++) {
                    double sigma_B = 1.0 / (b_alphas(l)+b_alphas(lp)); 
                    double boys = integral6d(sigma_A, sigma_B, Ra, Rb); 
                }
            }
            sum+= a_alphas(k) * a_alphas(kp)*b_alphas(l) b_alphas(lp) * boys; 
        }
    }
    return sum; 




}

double integral6d(double sigma_A, double sigma_B, arma::vec Ra, arma::vec Rb) {
    double U_A = std::pow(sigma_A*M_PI, 1.5); 
    double U_B = std::pow(sigma_B*M_PI, 1.5);
    double U = U_A*U_B;
    double V2 = 1/(sigma_A+sigma_B);
    double R = arma::dot(Ra - Rb, Ra - Rb);
    double T = V2 * std::pow(R,2);
    double boys; 
    if (Ra == Rb) {
        boys= U*std::sqrt(2*V2)*std::sqrt(2/M_PI);
    }
    else {
        boys = U*std::sqrt(1/(R*R)) *std::erf(std::sqrt(T));
    }
    
    
    return boys; 

}

