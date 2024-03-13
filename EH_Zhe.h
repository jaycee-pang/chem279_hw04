#if !defined EH_Zhe_H
#define EH_Zhe_H
#include <armadillo>
#include <cassert>
#include "AO_Zhe.h"

void Generate_Hmat(arma::mat &OV_mat, std::vector<AO> &AOs, arma::mat &H_mat);

double Solve_EH(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat, arma::vec &energy_vec, int num_ele);


#endif // EH_Zhe_H