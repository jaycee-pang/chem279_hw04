#pragma once
#include <iostream> 
#include <armadillo> 
#include <map> 
#include <cassert> 
#include "AO.h"
std::map<std::string, int> valence_map = {{"H",1}, {"C",4}, {"N", 5},{"O",6}, {"F", 7}};

std::map<std::string, int> element_map = {{"H",1}, {"C",6}, {"N", 7}, {"O",8}, {"F", 9}};
std::map<int, std::string> element_reverse = {{1, "H"}, {6, "C"}, {7, "N"}, {8, "O"}, {9, "F"}};


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
    std::string element;
    Atom(int Z, double x, double y, double z);
    void print_atom() const; 

}; 




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
        arma::mat S_; 
        int n_electrons;
   
        
    public: 


        Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms);
    
        // arma::vec C_alphas, H_alphas, O_alphas, N_alphas, F_alphas; 
        // arma::mat C2s_coeffs, C2p_coeffs, H1s_coeffs,O2s_coeffs, O2p_coeffs,
        //                     N2s_coeffs, N2p_coeffs, F2s_coeffs, F2p_coeffs; 
        std::vector<AO> generateAOs(); 
    
        // void overlap_matrix();
        void molecule_info() const;
        // void getBasis_data(const std::string&atom_type);

};
// arma::mat overlap_matrix(const std::vector<AO> &MoleculeAOs); 
void make_overlap_matrix(std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix);
Molecule read_mol(const std::string& molecule_name); 