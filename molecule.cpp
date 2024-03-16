#include "molecule.h"
#include "AO.h"
const std::map<std::string, int> valence_map = {{"H",1}, {"C",4}, {"N", 5},{"O",6}, {"F", 7}};

const std::map<std::string, int> element_map = {{"H",1}, {"C",6}, {"N", 7}, {"O",8}, {"F", 9}};
const std::map<int, std::string> element_reverse = {{1, "H"}, {6, "C"}, {7, "N"}, {8, "O"}, {9, "F"}};

Atom::Atom(int Z, double x, double y, double z): Z(Z), x(x), y(y), z(z) {
    if (Z == 1 || Z == 3) {
        valence_e = 1; 

    }
    else if (Z == 6 || Z == 7 || Z == 8 || Z == 9) {
        valence_e = Z - 2; 
    }
    auto it = element_reverse.find(Z); 
    if (it != element_reverse.end()) {
        element = it->second;
    }
    else {
        throw std::invalid_argument("Atom type not supported.");
    }
}
void Atom::print_atom() const {
    std::cout << "atom: " << element << ", atomic number: " << Z << ", ("<< x << "," << y<< ","<<z<< ")" << std::endl;
}

arma::vec C_alphas; arma::vec H_alphas; arma::vec O_alphas;
arma::vec N_alphas; arma::vec F_alphas;
arma::vec H1s_coeffs; 
arma::vec C2s_coeffs; arma::vec C2p_coeffs; 
arma::vec O2s_coeffs; arma::vec O2p_coeffs;
arma::vec N2s_coeffs; arma::vec N2p_coeffs;
arma::vec F2s_coeffs; arma::vec F2p_coeffs; 
// #include "AO.h"
// class AO;
// Molecule mol(molecule_name, natoms, charge, atoms);
Molecule::Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms) : name_(name), natoms_(n_atoms),
                                                                                charge_(charge), atoms_(atoms) {
    
    getBasis_data();

    int a = 0; int b = 0; 
    for (const Atom&atom : atoms_) {
        // getBasis_data();
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
    }  
     
    N_ = 4*a + b; 
    // AOs_.resize(N_);
    // for (const Atom& atom: atoms_) {
    //     atom.print_atom();
    // }
    AOs_ = generateAOs();
    assert (N_ == AOs_.size());
    S_.resize(N_, N_);
    H_.resize(N_, N_);
    C_.resize(N_, N_);


}
int Molecule::N() const {return N_;}
int Molecule::num_electrons() const { return n_electrons; }
int Molecule::natoms() const {return natoms_;}
const std::vector<Atom>& Molecule::atoms() const { return atoms_;}
const Atom& Molecule::get_atom(int i) const {return atoms_[i];}
const std::vector<AO>& Molecule::AOs() const {return AOs_;}
// for getting a certain type of basis function ie. s-type, or p-type 
const AO& Molecule::get_AO(std::string shell_type) const {

    // for (const AO& ao: AOs_) {
    //     if (ao.shell().back() == shell_type) {
    //         return ao;
    //     }
    // }
    for (const AO& ao: AOs_) {
        if (ao.shell() == shell_type) {
            return ao; 
        }
    }
    throw std::runtime_error("No AO with " + shell_type + " found.");

}

void getBasis_data() {
    std::string filename;
    std::vector<double> alphas;
    std::vector<double> coeffs_2s;
    std::vector<double> coeffs_p; // same for px, py, pz
    std::vector<std::string> STO3Gfiles = {"basis/C_STO3G.txt",
                                            "basis/H_STO3G.txt",
                                            "basis/N_STO3G.txt",
                                            "basis/O_STO3G.txt",
                                            "basis/F_STO3G.txt"};
    
    std::vector<std::string> atom_filenames = {"C", "H", "N", "O", "F"};

    for (size_t i = 0; i < STO3Gfiles.size(); ++i) {
        std::string filename = STO3Gfiles[i];

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }

        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<double> coeffs_p; // same for px, py, pz

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_p;
            if (!(iss >> alpha >> coeff_2s)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            
            if (atom_filenames[i] != "H") {
                if (!(iss >> coeff_p)) {
                    std::cerr << "Error reading line: " << line << std::endl;
                    continue;
                }
                coeffs_p.push_back(coeff_p);
            }
        }

        if (atom_filenames[i] == "C") {
            C_alphas = arma::vec(alphas);
            C2s_coeffs = arma::vec(coeffs_2s);
            C2p_coeffs = arma::vec(coeffs_p);
        } 
        else if (atom_filenames[i] == "H") {
            H_alphas = arma::vec(alphas);
            H1s_coeffs = arma::vec(coeffs_2s);
        } 
        else if (atom_filenames[i] == "O") {
            O_alphas = arma::vec(alphas);
            O2s_coeffs = arma::vec(coeffs_2s);
            O2p_coeffs = arma::vec(coeffs_p);
        } 
        else if (atom_filenames[i] == "N") {
            N_alphas = arma::vec(alphas);
            N2s_coeffs = arma::vec(coeffs_2s);
            N2p_coeffs = arma::vec(coeffs_p);
        }
        else if (atom_filenames[i] == "F") {
            F_alphas = arma::vec(alphas);
            F2s_coeffs = arma::vec(coeffs_2s);
            F2p_coeffs = arma::vec(coeffs_p);

        }
        else {
            throw std::invalid_argument("No STO3G basis set for this atom type."); 
        }

        file.close();
    }
}

std::vector<AO> Molecule::generateAOs() {
    std::vector<AO> AOs; 
    arma::uvec lmn(3);
    for (const Atom& atom : atoms_) {
        arma::vec R = { atom.x, atom.y, atom.z };
        switch (atom.Z) {
            case 1: 
                AOs.emplace_back("H1s", R, H_alphas, H1s_coeffs, arma::uvec{0, 0, 0});
                break;
            case 6: 
                AOs.emplace_back("C2s", R, C_alphas, C2s_coeffs, arma::uvec{0, 0, 0});
                for (size_t j = 0; j < 3; ++j) {
                    lmn.zeros();
                    lmn(j) = 1;
                    AOs.emplace_back("C2p", R, C_alphas, C2p_coeffs, lmn);
                }
                break;
            case 7:
                AOs.emplace_back("N2s", R, N_alphas, N2s_coeffs, arma::uvec{0, 0, 0});
                for (size_t j = 0; j < 3; ++j) {
                    lmn.zeros();
                    lmn(j) = 1;
                    AOs.emplace_back("N2p", R, N_alphas, N2p_coeffs, lmn);
                }
                break;
            case 8:
                AOs.emplace_back("O2s", R, O_alphas, O2s_coeffs, arma::uvec{0, 0, 0});
                for (size_t j = 0; j < 3; ++j) {
                    lmn.zeros();
                    lmn(j) = 1;
                    AOs.emplace_back("O2p", R, O_alphas, O2p_coeffs, lmn);
                }
                break;
            case 9: 
                AOs.emplace_back("F2s", R, F_alphas, F2s_coeffs, arma::uvec{0, 0, 0});
                for (size_t j = 0; j < 3; ++j) {
                    lmn.zeros();
                    lmn(j) = 1;
                    AOs.emplace_back("F2p", R, F_alphas, F2p_coeffs, lmn);
                }
                break;
            default:
                throw std::invalid_argument("No STO3G basis sets for this atom. Only compatible for H, C, N, O, F.");
        }
    }
    return AOs;
}

void Molecule::molecule_info() const {
    std::cout << name_ << ", # atoms: "<< natoms_ << ", basis functions: " << N_ << ", e- pairs: " << n_ <<std::endl; 
    // int natoms_; // number of atoms in the mol
    //     int N_; // num basis functions 
    //     int n_; // num AOs electron pairs = total/2 
    for (const Atom& atom: atoms_) {
        atom.print_atom();
    }
    for (const AO& ao: AOs_) {
        ao.print_AO(); 
    }

}


void Molecule::make_overlap_matrix() { // std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix
    int dim = AOs_.size();
    // overlap_matrix(dim, dim);
    for (int i = 0; i < dim; i++) { 
        for (int j = 0; j <= i; j++) {
            double overlap_elm = evaluate_contracted_overlap(AOs_[i], AOs_[j]); 
            std::cout << "AO i : " <<std::endl; AOs_[i].print_AO();
            S_(i,j) = overlap_elm;
            S_(j,i) = overlap_elm; 
        }
    }
    S_.print();
    // return overlap_matrix; 

}
// void make_H_mat(); 
// void H_mat(); 
const arma::mat& Molecule::S_overlap() const {
    return S_; 
}
int Molecule::routine() {
    // want this to be the driver function 
    return 0; 
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
        

        Molecule mol(molecule_name, natoms, charge, atoms); // will throw an error if not even # electron pairs 
        for (const Atom& atom: atoms) {
            atom.print_atom();
        }
        file.close();
        return mol;
    }

    catch (const std::exception& e) {
        std::cerr << "Error! " << e.what() << std::endl;
        file.close();

    }
    
}


// #endif 


// std::vector<AO> Molecule::generateAOs() {
//     std::vector<AO> AOs; 
//     arma::uvec lmn(3); 
//     for (const Atom& atom: atoms_) {
//         arma::vec R = {atom.x, atom.y, atom.z};
//         switch (atom.Z) {

        
//             if (atom.Z == 1) {
//                 AO ao(std::string("H1s"), R, H_alphas, H1s_coeffs, arma::uvec{0,0,0}); 
//                 AOs.push_back(ao);
//             }
//         else if (atom.Z == 6) {
//             AO ao2s(std::string("C2s"), R, C_alphas, C2s_coeffs, arma::uvec{0,0,0});
//             AOs.push_back(ao2s);
//             // get the correct lmn combo
//             for(size_t j = 0; j < 3; j++){
//                 arma::uvec lmn;
//                 lmn.zeros();
//                 lmn(j) = 1;
//                 AO ao2p(std::string("C2p"),R, C_alphas, C2p_coeffs, lmn);
//                 AOs.push_back(ao2p);
            
//             }

//         }
//         else if (atom.Z == 7) {
//             AO ao2s(std::string("N2s"), R, N_alphas, N2s_coeffs, arma::uvec{0,0,0});
//             AOs.push_back(ao2s);
//             // get the correct lmn combo
//             for(size_t j = 0; j < 3; j++){
//                 arma::uvec lmn;
//                 lmn.zeros();
//                 lmn(j) = 1;
//                 AO ao2p(std::string("C2p"),R, N_alphas, N2p_coeffs, lmn);
//                 AOs.push_back(ao2p);
            
//             }

//         }
//         else if (atom.Z == 8) {
//             AO ao2s(std::string("O2s"), R, O_alphas, O2s_coeffs, arma::uvec{0,0,0});
//             AOs.push_back(ao2s);
//             // get the correct lmn combo
//             for(size_t j = 0; j < 3; j++){
//                 arma::uvec lmn;
//                 lmn.zeros();
//                 lmn(j) = 1;
//                 AO ao2p(std::string("O2p"),R, O_alphas, O2p_coeffs, lmn);
//                 AOs.push_back(ao2p);
            
//             }


//         }
//         else if (atom.Z == 9) {
//             AO ao2s(std::string("F2s"), R, F_alphas, F2s_coeffs, arma::uvec{0,0,0});
//             AOs.push_back(ao2s);
//             // get the correct lmn combo
//             for(size_t j = 0; j < 3; j++){
//                 arma::uvec lmn;
//                 lmn.zeros();
//                 lmn(j) = 1;
//                 AO ao2p(std::string("F2p"),R, F_alphas, F2p_coeffs, lmn);
//                 AOs.push_back(ao2p);
            
//             }


//         }
//         else {
//             throw std::invalid_argument("No STO3G basis sets for this atom. Only compatable for H, C, N, O, F.");
//         }
//         }

//     }
//     return AOs;
// }