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
        
    public: 
        Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms) : name_(name), natoms_(n_atoms),
                                                                                     charge_(charge), atoms_(atoms) {}

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

void makeAOs() {

}
void makeMOs() {

}