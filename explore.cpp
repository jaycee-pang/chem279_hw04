#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "molecule.h"
#include "cndo.h"
int main(void) {

    getBasis_data();
    double R_NN = 1.098; 
    Atom N1(7, 0.0,0.0,0.0);
    Atom N2(7, 0.0, 0.0, R_NN); 
    std::vector <Atom> atoms = {N1, N2};
    // assume neutral 
    int charge = 0; 
    Molecule molecule("N2", atoms.size(), 0 , atoms);
    
    molecule.molecule_info();
    // molecule.S_.print();
    molecule.make_overlap_matrix();
    arma::mat overlap_mat = molecule.S_overlap(); 
    double tolerance = 10e-6; 
    int max_iterations = 50; 
    CNDO cndo(molecule, tolerance, max_iterations);
    cndo.SCF();
    return 0;
}