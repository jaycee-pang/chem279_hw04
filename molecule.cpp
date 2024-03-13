#include "molecule.h"
#include "AO.h"
// class AO;


Atom::Atom(int Z, double x, double y, double z): Z(Z), x(x), y(y), z(z),element(element_reverse[Z]) {}
void Atom::print_atom() const {
    std::cout << "atom : " << element << "atomic number: " << Z << ", ("<< x << "," << y<< ","<<z<< ")" << std::endl;
}


Molecule::Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms) : name_(name), natoms_(n_atoms),
                                                                                charge_(charge), atoms_(atoms) {
    
    
    int a = 0; int b = 0; 
    for (const Atom&atom : atoms) {
        getBasis_data(atom.element);
        if (atom.Z == 6 || atom.Z == 7 || atom.Z == 8 || atom.Z == 9) {
            a++; 
        }
        else if(atom.Z == 1) {
            b++;
        }
     
    }
    // epairs should be 1 for H, 2 pairs for C,N,O,F 
    // change for radicals? 
    int e_pairs = 2*a + (b/2); 
    if (e_pairs != int(e_pairs)) {
        throw std::runtime_error("The number of electron pairs should be even.");
    } 
    else {
        n_ = e_pairs;
        std::cout << "n pairs" << n_ << std::endl;
    } 
    AOs_= generateAOs(); 
    N_ = 4*a + b; 
    assert (N_ == AOs_.size());
    S_.resize(N_, N_);

}
        
std::vector<AO> Molecule::generateAOs() {
    std::vector<AO> AOs; 
    arma::uvec lmn; 
    for (const Atom& atom: atoms_) {
        if (atom.Z == 1) {
            AO ao("H1s", arma::vec{atom.x, atom.y, atom.z}, H_alphas, H1s_coeffs, arma::uvec{0,0,0}); 
            AOs.push_back(ao);
        }
        else if (atom.Z == 6) {
            AO ao2s("C2s", arma::vec{atom.x, atom.y, atom.z}, C_alphas, C2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p("C2p",{atom.x, atom.y, atom.z}, C_alphas, C2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }

        }
        else if (atom.Z == 7) {
            AO ao2s("N2s", arma::vec{atom.x, atom.y, atom.z}, N_alphas, N2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p("C2p",{atom.x, atom.y, atom.z}, N_alphas, N2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }

        }
        else if (atom.Z == 8) {
            AO ao2s("O2s", arma::vec{atom.x, atom.y, atom.z}, O_alphas, O2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p("O2p",{atom.x, atom.y, atom.z}, O_alphas, O2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }


        }
        else if (atom.Z == 9) {
            AO ao2s("F2s", arma::vec{atom.x, atom.y, atom.z}, F_alphas, F2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p("F2p",{atom.x, atom.y, atom.z}, F_alphas, F2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }


        }
        else {
            throw std::invalid_argument("No STO3G basis sets for this atom. Only compatable for H, C, N, O, F.");
        }

    }
    return AOs;
}



void Molecule::molecule_info() const {
    std::cout << name_ << ":" << std::endl; 
    for (const Atom& atom: atoms_) {
        atom.print_atom();
    }
    for (const AO& ao: AOs_) {
        ao.print_AO(); 
    }

}

Molecule read_mol(const std::string& molecule_name) {
    std::string filename = molecule_name+".txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error opening the file " + filename);
    } 
    try {
        
        int natoms; 
        int charge; 
        std::vector<Atom> atoms;  
        std::string line; 
  
        if (std::getline(file, line)) {
            std::istringstream iss(line);
            if (!(iss >> natoms >> charge)) {
                throw std::runtime_error("Error reading molecule atoms and charge");
            }
        } 
        else {
            throw std::runtime_error("Empty file or unable to read the first line.");
        }
        int Z; 
        double x, y, z;
        while (file >> Z >> x >> y >> z) {

            atoms.emplace_back(Z, x, y, z);
        }
        if (atoms.size() != natoms) {
            throw std::runtime_error("Atom information not consistent.");
        }
        Molecule mol(filename, natoms, charge, atoms); // will throw an error if not even # electron pairs 
        file.close();
        return mol;
    }

    catch (const std::exception& e) {
        std::cerr << "Error! " << e.what() << std::endl;
        file.close();

    }
    
}




void make_overlap_matrix(std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix) {
    int dim = MoleculeAOs.size();
    // overlap_matrix(dim, dim);
    for (int i = 0; i < dim; i++) { 
        for (int j = 0; j <= i; j++) {
            double overlap_elm = evaluate_contracted_overlap(MoleculeAOs[i], MoleculeAOs[j]); 
            overlap_matrix(i,j) = overlap_elm;
            overlap_matrix(j,i) = overlap_elm; 
            // overlap += Ns_[i] * other.Ns_[j]*ds_[i] * other.ds_[j] * primitive_overlap;
        }
    }
    // return overlap_matrix; 

}
