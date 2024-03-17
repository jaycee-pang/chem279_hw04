#include "cndo.h"
std::map <std::string, double> semi_emp_s = {{"H1s", 7.716}, {"C2s", 14.051}, {"N2s", 19.316},
                                            {"O2s", 25.390}, {"F2s", 32.272}};
std::map <std::string, double> semi_emp_p = {{"C2p", 5.572}, {"N2p", 7.275},
                                            {"O2p", 9.111}, {"F2p", 11.080}};
std::map <std::string, double> semi_empirical = {{"H1s", 7.716}, {"C2s", 14.051}, {"N2s", 19.316},
                                            {"O2s", 25.390}, {"F2s", 32.272}, {"C2p", 5.572}, {"N2p", 7.275},
                                            {"O2p", 9.111}, {"F2p", 11.080}}; 
std::map <std::string, int> atomic_bonding = {{"H", 9}, {"C", 21}, {"N", 25}, {"O", 31}, {"F", 39} };
const double au_to_EV = 27.211;
CNDO::CNDO(Molecule& molecule, double tolerance, int max_it) : molecule(molecule), tol(tolerance), max_it(max_it) {
    N_ = molecule.N();
    natoms_ = molecule.natoms(); 
    n_electrons_ = molecule.num_electrons();
    // make sure S is made, To-Do: bool if S is made in molecule class 
    molecule.make_overlap_matrix();
    S_ = molecule.S_overlap(); 
    std::cout << "total electrons: " << n_electrons_ << std::endl;
    q_ = n_electrons_ /2; 
    p_ = (n_electrons_-q_); 
    std::cout << "p: " << p_ << ", q: " << q_ << std::endl;
       // find_ooc();
    Gamma_.resize(natoms_, natoms_); // G only relies on atoms A and B 
    Ga_.zeros(); 
    Gb_.zeros();
    
    Gamma_ *= au_to_EV; 
    Pa_.resize(N_, N_); Pb_.resize(N_, N_);
    Ga_.resize(N_, N_); Gb_.resize(N_, N_);
    build_G(); 
    
    Eel_ = 0.0; Etot_ = 0.0; 
    calc_Ec(); 
    // std::cout << "Ec: " << Ec_ << std::endl;
    Ea_.resize(N_); Eb_.resize(N_);
    H_.resize(N_, N_);
    H_.zeros(N_, N_);
    buildH(); 
    H_.print();


}
int CNDO::find_ooc() {
   
    int p = 0; int q = 0; 
    int charge = molecule.charge(); 
    for (const Atom& atoms : molecule.atoms()) {
        int single = (charge < 0) ? -charge: 0; 
        p += single; 
        q += (charge < 0 ) ? 0 : single; 
    }
    p_ = p;
    q_=q; 

} 
double CNDO::gamma_ab(const AO &ao1, const AO& ao2) {
    arma::uvec lmn1 = ao1.lmn(); arma::uvec lmn2 = ao2.lmn(); 
    arma::uvec lmns = {0, 0, 0}; // s-shell only 
    // if (lmn1!= lmns || lmn2!= lmns) {
    //     throw std::invalid_argument("The two-electron integral evaluation occurs between s-shell only."); 
    // }
    
    arma::vec da(ao1.ds().size()); arma::vec db(ao2.ds().size()); 
    arma::vec a_alphas = ao1.alphas(); arma::vec b_alphas = ao2.alphas();
    assert(a_alphas.size() == b_alphas.size());
    arma::vec Ra = ao1.R(); arma::vec Rb = ao2.R();
    double sum = 0.0; 
    // iterate over k, k', l, l' to get all primitives
    for (size_t k = 0; k < a_alphas.size(); k++) {
        for (size_t kp = 0; kp < a_alphas.size(); kp++) {
            double sigma_A = 1.0/(a_alphas(k) + a_alphas(kp));
            for (size_t l=0; l < b_alphas.size(); l++) {
                for (size_t lp =0; lp < b_alphas.size(); lp++) {
                    double sigma_B = 1.0 / (b_alphas(l)+b_alphas(lp)); 
                    double boys = integral6d(sigma_A, sigma_B, Ra, Rb); 
                    sum+= a_alphas(k) * a_alphas(kp)*b_alphas(l)* b_alphas(lp) * boys; 
                }
            }
            
        }
    }
    return sum; 

}

double CNDO::integral6d(double sigma_A, double sigma_B, arma::vec Ra, arma::vec Rb) {
    double U_A = std::pow(sigma_A*M_PI, 1.5); 
    double U_B = std::pow(sigma_B*M_PI, 1.5);
    double U = U_A*U_B;
    double V2 = 1/(sigma_A+sigma_B);
    // double R = arma::dot(Ra - Rb, Ra - Rb);
    double R = arma::norm(Ra-Rb, 2);
    double T = V2 * std::pow(R,2);
    double boys; 
    
    if (R==0) {
        boys= U*std::sqrt(2*V2)*std::sqrt(2/M_PI);
    }
    else {
        boys = U*std::sqrt(1/(R*R)) *std::erf(std::sqrt(T));
    }
    
    
    return boys; 

}
// maybe move G calculations to molecule class since it only has to do with atoms 
void CNDO::build_G() {
    // use elements of the Pa Pb total density matrix to update
    /* 
    for i < num atoms 
        for j < atoms 
            get AO i, get AO j
            evlauate the 2 electron integral s shell electrons 
    */

    // Gamma_.zeros();
    // for (int i = 0; i < natoms_; i++) {
    //     for (int j = 0; j < natoms_; j++) {
    //         double gamma_ij = compute_gamma_ij(i, j); 
    //         Gamma_(i, j) = gamma_ij;
    //     }
    // }
    
    
    for (int i = 0; i < natoms_; i ++) {
        for (int j = 0; j < natoms_; j++) {
            const Atom& atom_i = molecule.get_atom(i);
            const Atom& atom_j = molecule.get_atom(j);
            std::string atomi_shell = atom_i.element + (atom_i.Z == 1 ? "1s" : "2s");
            // std::cout << "atom i " << atomi_shell << std::endl;
            const AO& ao_i = molecule.get_AO(atomi_shell);
            std::string atomj_shell = atom_j.element + (atom_j.Z == 1 ? "1s" : "2s");
            const AO& ao_j = molecule.get_AO(atomj_shell);
            // std::cout << "atom j " << atomj_shell <<std::endl;
            Gamma_(i,j) = gamma_ab(ao_i, ao_j);
            // std::cout << atomi_shell << ", " << atomj_shell << std::endl;

        }
    }
    // std::cout << Gamma_ << std::endl;

}
 
void CNDO::build_F() {
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
            double pAA_beta = Pb_(AO_mu, AO_mu);
            double G_AA = Gamma_(i,i); // like gamma_AA 
            
            Fa_(AO_mu, AO_mu) = -IA + ((pAA_alpha - ZA) - (pAA_alpha - 0.5))  * G_AA; 
            Fb_(AO_mu, AO_mu) = -IA + ((pAA_alpha - ZA) - (pAA_beta - 0.5))  * G_AA;

            // for B!=A 
            for (int j = 0; j < natoms_; j++) {
                
                if (i!=j) {
                    const Atom& atomB = molecule.get_atom(j);
                    const std::vector<AO> BAOs = molecule.atom_AOs(atomB.element); 
                    double G_AB = Gamma_(i,j); 
                    double betaB = atomic_bonding[atomB.element];
                    double betaAB = 0.5 * (betaA + betaB);

                    for (const AO& ao_B : BAOs) {
                        for (int AO_nu = 0; AO_nu < N_; AO_nu++) {
                            double p_BB = Pb_(AO_nu, AO_nu);
                            double ZB = atomB.valence_e;
                            double gamma_AB = Gamma_(i, j);
                            if (AO_mu != AO_nu) {
                                
                                Fa_(AO_mu, AO_nu) += betaAB * S_(AO_mu, AO_nu) - Pa_(AO_mu, AO_nu) * G_AB; 
                                Fb_(AO_mu, AO_nu) += betaAB * S_(AO_mu, AO_nu) - Pb_(AO_mu, AO_nu) * G_AB;
                            }
                            if (AO_mu == AO_nu) {
                                // add the last term in 
                                Fa_(AO_mu, AO_nu) += (p_BB - ZB) * gamma_AB;
                                Fb_(AO_mu, AO_nu) += (p_BB - ZB) * gamma_AB;
                            }
                        }
                    }
                }

            }  
            
            AO_mu++;
        }
    }


}
// core Hamiltonian
void CNDO::buildH() {
    // need ZA atomic numbers (valence e numbers)
    // beta_AB*S_munu
    H_.resize(N_, N_);
    // 

    for (int mu = 0; mu < N_; mu++) {
        for (int nu = 0; nu < N_; nu++) {
            double H_mu_nu = 0.0;
            for (int i = 0; i < natoms_; i++) {
                const Atom& atomA = molecule.get_atom(i);
                // atomA.print_atom(); 
               
                const std::vector<AO>& AAOs = molecule.atom_AOs(atomA.element);
                // std::cout << "AAos size: " << AAOs.size() << std::endl;
                for (const AO& ao_A : AAOs) {
                    // ao_A.print_AO();
                    double IA = semi_empirical[ao_A.shell()];
                    
                    int AO_mu = i * atomA.n_basis + mu; 
                    // ZA
                    H_mu_nu += -0.5 * (IA + atomA.valence_e - 0.5) * Gamma_(i, i) * S_(mu, nu); // IA , gamma_AA
                    for (int j = 0; j < natoms_; j++) {
                        if (i != j) {
                            const Atom& atomB = molecule.get_atom(j);
                            double ZB = atomB.valence_e;
                            double G_AB = Gamma_(i, j);
                            
                            double betaA = atomic_bonding[atomA.element];
                            double betaB = atomic_bonding[atomB.element];
                            double beta_AB = 0.5 * (betaA + betaB);
                            H_mu_nu += 0.5 * (betaA + betaB) * S_(mu, nu) * (1 - (mu == nu)) * G_AB; // beta_AB * G_AB
                        }
                    }
                }
            }
            H_(mu, nu) = H_mu_nu;
        }
    }


}
       
        

void CNDO::update_P() {
    // update Pa 
    Pa_ = Ca_.cols(0, p_-1) *trans(Ca_.cols(0, p_-1)); 
    // update Pb if q > 0 
    if (q_>0) {
        Pb_ = Cb_.cols(0, q_-1) * trans(Cb_.cols(0, q_-1)); 
    }
    else {
        Pb_.zeros(); 
    }
}
void CNDO::G_ab() {
    
    arma::vec P_tot = arma::zeros(natoms_);
    int AO_mu = 0; 
    // get total density 
    for (int mu = 0; mu < natoms_; mu++) {
        const Atom& atom = molecule.get_atom(mu);
        const std::vector<AO> atomAOs = molecule.atom_AOs(atom.element); 
        for (const AO& ao: atomAOs) {
            P_tot += Pa_(AO_mu, AO_mu) + Pb_(AO_mu, AO_mu); 
            AO_mu++;
        }
    }
    AO_mu = 0;
    for (int mu = 0; mu < natoms_; mu++) {
        double G_AA = Gamma_(mu, mu);
        const Atom& atom_u = molecule.get_atom(mu);
        const std::vector<AO> atomuAOs = molecule.atom_AOs(atom_u.element); 
        for (const AO& ao: atomuAOs) {
            Ga_(AO_mu, AO_mu) = (P_tot(mu) - Pa_(AO_mu, AO_mu)) * G_AA;
            Gb_(AO_mu, AO_mu) = (P_tot(mu) - Pb_(AO_mu, AO_mu)) * G_AA;

            int AO_nu = 0;
            for (int nu = 0; nu < natoms_; nu++) {
                double G_AB = Gamma_(mu, nu);
                const Atom& atom_v = molecule.get_atom(nu);
                const std::vector<AO> atomvAOs = molecule.atom_AOs(atom_v.element); 
                if (mu != nu){
                    Ga_(AO_mu, AO_mu) += P_tot[nu] * G_AB;
                    Gb_(AO_mu, AO_mu) += P_tot[nu] * G_AB;
                }
                for (const AO& ao: atomvAOs) {
                    if (mu != nu) {
                        Ga_(AO_mu, AO_nu) = -G_AB * Pa_(AO_mu, AO_nu);
                        Gb_(AO_mu, AO_nu) = -G_AB * Pb_(AO_mu, AO_nu);
                    }
                    AO_nu++;
                }
            }
            AO_mu++; 
        }
    }

}

void CNDO::calc_Ec() {
    Ec_  = 0.0; 
    for (int i = 0; i < natoms_; i++) {
        const Atom& atom_i = molecule.get_atom(i); 
        for (int j = i+1; j<natoms_; j++) {
            
            const Atom& atom_j = molecule.get_atom(j);
            arma::vec Ra = {atom_i.x, atom_i.y, atom_i.z}; 
            arma::vec Rb = {atom_j.x, atom_j.y, atom_j.z}; 
            // arma::vec R = arma::dot(Ra - Rb, Ra - Rb);
            double R = arma::norm(Ra - Rb, 2);
            int ZA = atom_i.valence_e; 
            int ZB = atom_j.valence_e; 
            Ec_ += ZA * ZB / R; 
        
        }
    }
    Ec_ /= 2.0; // e piars
    Ec_ *= au_to_EV; 

}

double CNDO::calculate_E() {
    // 
    // Etot_ = 0.0; 
    Eel_ += 0.5 * arma::dot(Pa_, H_ + Fa_) + 0.5 * arma::dot(Pb_, H_ + Fb_);
    Etot_ += Ec_ + Eel_;
 
    // Etot_ = Eel_ + Ec_; 
    return Etot_; 
    
}
// void CNDO::update_F() {
//     update_P();
//     Fa_ = H_ + 2 * Ga_ - Gb_;
//     // F_beta = H + 2 * Gb - Ga
//     Fb_ = H_ + 2 * Gb_ - Ga_; 


// }

// driver code 
void CNDO::SCF() {
    std::cout << "SCF:\nmax iterations: " << max_it << ", tolerance: " << tol << std::endl; 
    // initial guess for P, guess all 0
    Pa_.zeros();
    Pb_.zeros();
    print_CNDO_info();
    // std::cout << "initial guess: \nPa: " << Pa_ << "\nPb:" << Pb_ << std::endl; 
    // Fa, Fb 
    // place holders 
    Fa_ = H_ + 2 * Ga_ - Gb_;
    Fb_ = H_ + 2 * Gb_ - Ga_;
    // Fa_ = H_ + Ga_; 
    // Fb_ = H_ + Gb_; 
    build_F();
    // std::cout << "Fa" << Fa_ << std::endl;
    
    // build core hamiltonian
    for (int i = 0; i < max_it; i++) {
        // G_ab(); // this does for current P 
        // build Fock matrix with iniitla guess P 
        // build_F(); 
        G_ab();
        build_F();
        
        
        std::cout << "Iteration: " << i << std::endl;
        

        arma::eig_sym(Ea_, Ca_, Fa_); 
        arma::eig_sym(Eb_, Cb_, Fb_);
        arma::mat Pa_old = Pa_; 
        arma::mat Pb_old = Pb_; 
        update_P();

        // solve eigenvalue problesmm to get MO coeff and eigenvalues FC=Ceps
        // make new density matrices 
        //      Pa_ vs Pa_old, Pb_ vs Pb_old
        // check for convergence 
        print_SCF(); 
        if (arma::approx_equal(Pa_, Pa_old, "absdiff", tol) && arma::approx_equal(Pa_, Pb_old, "absdiff", tol)){
            std::cout << "SCF converged at iteration: " << i << std::endl;
            break; 
        }

    }
    Etot_ = calculate_E(); 
    std::cout << "Molecule " << molecule.name() << std::endl;
    std::cout << "Nuclear Repulsion Energy is "<< Eel_<< " eV" << std::endl;
    std::cout << "Electron Energy is " << Ec_ << " eV" << std::endl;
    std::cout << "Total energy is " << Etot_ << " eV" << std::endl; 
    


} 
void CNDO::print_CNDO_info() const {
    std::cout << "gamma:\n" << Gamma_ << std::endl;
    std::cout << "overlap:\n" << S_ << std::endl;
    std::cout << "p: " << p_ << ", q: " << q_ << std::endl;
    std::cout << "core Hamiltonian:\n" << H_ << std::endl; 
    
    
}
void CNDO::print_SCF() const {
    std::cout << "Fa:\n" << Fa_ << "Fb:\n" << Fb_<< std::endl; 
    std::cout << "Ca:\n" << Ca_ << "Cb:\n" << Cb_ << std::endl;
    std::cout << "Ea:\n" << Ea_ << "Eb:\n" << Eb_ << std::endl;
    std::cout << "Pa:\n" <<Pa_ << "Pb:\n" << Pb_ << std::endl;
    std::cout << "Ga:\n" <<Ga_ << "Gb:\n" << Gb_ << std::endl;


}

