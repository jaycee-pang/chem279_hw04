#pragma once
#include "molecule.h"
#include <iostream> 
#include <armadillo> 
#include <map> 
#include <cassert> 
#include <algorithm> 
#include <iterator> 
#include "AO.h"
// class AO;
extern const std::map<std::string, int> valence_map;
extern const std::map<std::string, int> element_map;
extern const std::map<int, std::string> element_reverse;



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
    int valence_e;
    int n_basis; 

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
        std::vector<AO> AOs_; // basis functions 
        arma::mat S_;  // overlap matrix 
        arma::mat H_; // H matrix
        arma::mat C_; // coeff matrix 
        int n_electrons;
        
    public: 
        Molecule(std::string name, int n_atoms, int charge, std::vector<Atom> atoms);
        std::vector<AO> generateAOs(); 
        // void overlap_matrix();
        void molecule_info() const;
        void make_overlap_matrix(); 
        const arma::mat& S_overlap() const;
        int routine(); 
        int N() const;
        int num_electrons() const;
        int natoms() const; 
        const std::vector<Atom>& atoms() const; 
        const Atom& get_atom(int i) const;
        const std::vector<AO>& AOs() const; 
        const AO& get_AO(std::string label) const; 
        const std::vector<AO> atom_AOs(std::string atom_element) const; 
        int charge() const; 
        std::string name() const; 
        // void make_H_mat(); 
        // void H_mat(); 
   

};
// void make_overlap_matrix(std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix);
void getBasis_data(); 
// arma::mat overlap_matrix(const std::vector<AO> &MoleculeAOs); 
// void make_overlap_matrix(std::vector<AO> &MoleculeAOs, arma::mat &overlap_matrix);
Molecule read_mol(const std::string& molecule_name); 
// #endif