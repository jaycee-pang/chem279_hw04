#include <iostream> 
#include <armadillo> 
#include <map> 
std::map<std::string, int> valence_map = {{"H",1}, {"C",4}, {"N", 5},
                                         {"O",6}, {"F", 7}};

std::map<std::string, int> element_map = {{"H",1}, {"C",4}, {"N", 5},
                                         {"O",6}, {"F", 7}};
const double H1s = -13.6;
const double C2s = -21.4;
const double C2px = -11.4;
const double C2py = -11.4; 
const double C2pz = -11.4; 
const double K = 1.85; 
const double convertBohr = 0.52917706; // 1 Bohr = 0.52917706 Angstrom 
struct Atom {
    int N;  // number of basis functions  N = 4a+b
    int Z;  // atomic numer 
    double x, y, z;     // xyz coords 
    Atom(int Z, double x, double y, double z): Z(Z), x(x), y(y), z(z) {}

}; 

std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    os << "Atomic Number: " << atom.Z <<std::endl;
    os << "(" << atom.x << ", " << atom.y << ", " << atom.z << ")";
    return os;
}


class Molecule{
    private: 
        int natoms_; // number of atoms in the mol
        int N_; // num basis functions 
        int n_; // num AOs electron pairs = total/2 
        int charge_; // 
        std::string name_; // txt file molecule name 
        std::vector<Atom> atoms_; // vector of Atom objects (Z, x,y,z)
        arma::vec Zs_; // vector of the atomic numbers for atoms
        std::vector<std::string> atom_elements; // string atom names 
        std::vector<AO> AOs_; 
        
    public: 
        Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms) : name_(name), natoms_(n_atoms),
                                                                                     charge_(charge), atoms_(atoms) {
            getBasis_data();
            AOs_= generateAOs(); 
                                                                                     }
        std::vector<AO> generateAOs() {
            std::vector<AO> AOs; 
            for (const Atom& atom: atoms_) {
                if (atom.Z == 1) {
                    AO ao("H1s", arma::vec{atom.x, atom.y, atom.z}, H_alphas, H1s_coeffs, arma::vec{0,0,0}); 
                    AOs.push_back(ao);
                }
                else if (atom.Z == 6) {
                    AO ao2s("C2s", arma::vec{atom.x, atom.y, atom.z}, C_alphas, C2s_coeffs, arma::vec{0,0,0});
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
                    AO ao2s("N2s", arma::vec{atom.x, atom.y, atom.z}, N_alphas, N2s_coeffs, arma::vec{0,0,0});
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
                    AO ao2s("O2s", arma::vec{atom.x, atom.y, atom.z}, O_alphas, O2s_coeffs, arma::vec{0,0,0});
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
                    AO ao2s("F2s", arma::vec{atom.x, atom.y, atom.z}, F_alphas, F2s_coeffs, arma::vec{0,0,0});
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
                    throw invalid_argument("No STO3G basis sets for this atom. Only compatable for H, C, N, O, F.");
                }
                
                // const std::string& shell, const arma::vec& R, const arma::vec& alphas, const arma::vec& d_coeff, const arma::uvec& lmn)

            }
            return AOs;
        }

};


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

