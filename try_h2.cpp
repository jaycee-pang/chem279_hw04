#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "molecule.h"

// #include "AO_Zhe.h"
// #include "EH_Zhe.h"
int main(int argc, char* argv[]) {

    getBasis_data();
    // arma::mat H_STO3G, C_STO3G;
    // H_STO3G.load("H_STO3G.txt");
    // C_STO3G.load("C_STO3G.txt");
    // std::vector<AO> MoleculeAOs;
    // int num_ele;
    
    if (argc !=2)
    {
        printf("usage ./hw3 filename, for example hw3 example.txt\n");
        return 1;
    }
    std::string fname(argv[1]);
    // try {
    //     num_ele = GenerateAOs(MoleculeAOs, fname, H_STO3G, C_STO3G);
    // }
    // catch (std::invalid_argument &e)
    // {
    //     std::cerr << e.what() << std::endl;
    //     return EXIT_FAILURE;
    // }
    // PrintAOs(MoleculeAOs);

    // int dim = MoleculeAOs.size();
    // arma::mat OV_mat(dim, dim);
    // Eval_OV_mat(MoleculeAOs, OV_mat);
    // OV_mat.print("Overlap matrix:");

    // arma::mat H_mat(dim, dim);
    // Generate_Hmat(OV_mat, MoleculeAOs, H_mat);
    // H_mat.print("Hamiltonian matrix");

    // arma::mat C_mat(dim, dim);
    // arma::vec energy_vec(dim);
    // double Energy = Solve_EH(OV_mat, H_mat, C_mat, energy_vec, num_ele);
    // C_mat.print("MO coefficients (C matrix):");
    // // check the MOs are orthonormal
    // arma::mat MO_ov = C_mat.t()*OV_mat*C_mat;
    // MO_ov.print("MO overlap matrix:");
    // // energy_vec.print("energy_vec");
    // printf("The molecule in file %s has energy %f\n", argv[1], Energy);

    // return EXIT_SUCCESS;
    return 0;
}