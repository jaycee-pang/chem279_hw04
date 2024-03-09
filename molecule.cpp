#include <iostream> 
#include <armadillo> 
const double H1s = -13.6;
const double C2s = -21.4;
const double C2px = -11.4;
const double C2py = -11.4; 
const double C2pz = -11.4; 
const double K = 1.85; 
const double convertBohr = 0.52917706; // 1 Bohr = 0.52917706 Angstrom 
struct Atom {
    int N;  // number of basis functions  C: 4, H: 1
    int Z;  // atomic numer 
    double x, y, z;     // xyz coords 
    Atom(int Z, double x, double y, double z);

}; 

class Molecule{
    private: 
        int natoms_; 
        int N_; // num basis functions 
        int n_; // num AOs electron pairs 
        int charge_; 
        std::string name_; 
        
        arma::vec Zs_;
        std::vector<std::string> atom_elements; 
        
    public: 

};
void makeAOs() {

}
void makeMOs() {

}