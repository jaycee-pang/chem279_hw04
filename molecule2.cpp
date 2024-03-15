#include "molecule.h"


Atom::Atom(int Z, double x, double y, double z): Z(Z), x(x), y(y), z(z),element(element_reverse[Z]) {}
void Atom::print_atom() const {
    std::cout << "atom : " << element << "atomic number: " << Z << ", ("<< x << "," << y<< ","<<z<< ")" << std::endl;
}
// #ifndef MOLECULE_H
// #define MOLECULE_H
// #include "AO.h"
class AO;
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
        arma::vec R = {atom.x, atom.y, atom.z};
        if (atom.Z == 1) {
            AO ao(std::string("H1s"), R, H_alphas, H1s_coeffs, arma::uvec{0,0,0}); 
            AOs.push_back(ao);
        }
        else if (atom.Z == 6) {
            AO ao2s(std::string("C2s"), R, C_alphas, C2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p(std::string("C2p"),R, C_alphas, C2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }

        }
        else if (atom.Z == 7) {
            AO ao2s(std::string("N2s"), R, N_alphas, N2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p(std::string("C2p"),R, N_alphas, N2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }

        }
        else if (atom.Z == 8) {
            AO ao2s(std::string("O2s"), R, O_alphas, O2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p(std::string("O2p"),R, O_alphas, O2p_coeffs, lmn);
                AOs.push_back(ao2p);
            
            }


        }
        else if (atom.Z == 9) {
            AO ao2s(std::string("F2s"), R, F_alphas, F2s_coeffs, arma::uvec{0,0,0});
            AOs.push_back(ao2s);
            // get the correct lmn combo
            for(size_t j = 0; j < 3; j++){
                lmn.zeros();
                lmn(j) = 1;
                AO ao2p(std::string("F2p"),R, F_alphas, F2p_coeffs, lmn);
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

void Molecule::make_overlap_matrix(std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix) {
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





void getBasis_data(const std::string&atom_type) {
    std::string filename;
    
    
    if (atom_type == "C") {
        filename = "C_STO3G.txt"; 
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        C_alphas = arma::vec(alphas);
        C2s_coeffs = arma::vec(coeffs_2s);
        C2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            C2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
        
    }
    else if(atom_type == "H") {
        filename = "H_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_1s;
        // std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_1s; 
            if (!(iss >> alpha >> coeff_1s)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_1s.push_back(coeff_1s);
            // coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        H_alphas = arma::vec(alphas);
        H1s_coeffs = arma::vec(coeffs_1s);


        file.close();
    }
    else if (atom_type == "O") {
        filename = "O_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        O_alphas = arma::vec(alphas);
        O2s_coeffs = arma::vec(coeffs_2s);
        O2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            O2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else if(atom_type == "N") {
        filename = "N_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        N_alphas = arma::vec(alphas);
        N2s_coeffs = arma::vec(coeffs_2s);
        N2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            N2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else if(atom_type == "F") {
        filename = "F_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        N_alphas = arma::vec(alphas);
        N2s_coeffs = arma::vec(coeffs_2s);
        N2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            N2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else {
        throw std::runtime_error("No basis set for this atom type ");
    }

}
// #endif 