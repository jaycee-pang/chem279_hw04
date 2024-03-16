#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "molecule.h"
#include "cndo.h"
int main(int argc, char* argv[]) {

    getBasis_data();
    if (argc !=2)
    {
        printf("usage ./hw3 filename, for example hw3 example.txt\n");
        return 1;
    }
    std::string fname(argv[1]);
    std::string path_file = "sample_input/"+fname;
    Molecule molecule = read_mol(path_file);
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