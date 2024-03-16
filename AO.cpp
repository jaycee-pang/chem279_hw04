#include "AO.h"
// #ifndef AO_H
// #define AO_H
// #include "molecule.h"
// class Molecule; 
// class AO;

AO::AO(const std::string& shell, const arma::vec& R, const arma::vec& alphas, const arma::vec& d_coeff, const arma::uvec& lmn)
        : shell_type(shell), R_(R), alphas_(alphas), ds_(d_coeff), lmn_(lmn) {
            // normalize and combine these N with the contraction coefficient term 
            for (size_t k = 0; k < alphas_.n_elem; k++) {
                double self_overlap = overlap3d(R_, R_, alphas_(k), alphas_(k), lmn_, lmn_); 
                ds_(k) /= std::sqrt(self_overlap); 
            }

        }
  
  
void AO::print_AO() const {
    // std::cout << atom_.print_atom()<< ", " << shell_type<< ", " << lmn_.print() << std::endl;
    std::cout << shell_type << ", angular momentum: (" << lmn_(0)<< ", " << lmn_(1) << ", " << lmn_(2) << ")"  <<std::endl;
    std::cout << "alphas: ";
    for (int i = 0; i < alphas_.n_elem; i++) {
        std::cout << alphas_(i) << " ";
    }
    std::cout << "\ncontraction coefficients: ";
    for (int i = 0; i < ds_.n_elem; i++) {
        std::cout << ds_(i) << " ";
    }
    std::cout<<std::endl;
}




std::string AO::shell() const {return shell_type;} 
arma::vec AO::R() const{return R_;}
arma::vec AO::alphas() const {return alphas_;};
arma::vec AO::ds() const {return ds_;}; 
arma::uvec AO::lmn() const {return lmn_;}
double overlap1d(double xa, double xb, double alpha, double beta, int la, int lb) {
    double dx = (xa - xb);
    
    double exponential_term = exp(-(alpha*beta*dx*dx)/(alpha+ beta))* sqrt(M_PI / (alpha+ beta));
    double Rp = (alpha * xa + beta * xb)/ (alpha + beta); // new center 
    double result = 0.0; 
    
    for (int i = 0; i <= la; i++) 
        for (int j = 0; j <= lb; j++) {
            if ((i +j)%2 == 1) 
                continue; 
            double double_fact = double_factorial(i +j-1);
            double coeffs = choose(la, i) * choose(lb, j);
            double num = std::pow(Rp-xa, (la -i)) * std::pow(Rp -xb, (lb-j));
            double denom = std::pow(2*(alpha+beta), (double(i+j)/2.0));
            double term = coeffs * double_fact * num / denom; 
            result += term;
            
        }
    
    result *= exponential_term; //*result; // *= exponential_term; 
    return result; 
}
double overlap3d(const arma::vec& Ra, const arma::vec& Rb, double alphas, double betas, const arma::uvec& la, const arma::uvec& lb)  {
    double Sabx = overlap1d(Ra(0), Rb(0), alphas, betas, la(0), lb(0));
    double Saby = overlap1d(Ra(1), Rb(1), alphas, betas, la(1), lb(1));
    double Sabz = overlap1d(Ra(2), Rb(2), alphas, betas, la(2), lb(2));
    return Sabx * Saby* Sabz; 
}

double evaluate_contracted_overlap(const AO& ao1, const AO& ao2) {
    double sum = 0.0; 
    arma::vec alphas = ao1.alphas(); arma::vec betas = ao2.alphas(); 
    arma::uvec la = ao1.lmn(); arma::uvec lb = ao2.lmn();
    arma::vec da = ao1.ds(); arma::vec db = ao2.ds(); 
    arma::vec Ra = ao1.R(); arma::vec Rb = ao2.R(); 
    for (size_t k = 0; k < alphas.n_elem; k++ ) {
        double alpha_k = alphas(k); 
        for (size_t j = 0; j < alphas.n_elem; j++) {
            double beta_j = betas(j);
            double overlap = overlap3d(Ra, Rb, alpha_k, beta_j, la, lb); 
            sum += da(k) * db(j) * overlap; 
        }
        
    }
    return sum; 
}
// #endif 
double choose(int n, int k) { // n! / k!*(n-k!)!
    if(n < k || k < 0) {
        throw std::invalid_argument("Combination: the number of elements should be bigger than selection numbers AND two numbers should be positive\n");
    }
    double result = 1.0;
    double n_d = n;
    
    if (k > n - k) {
        k = n - k;
    }

    for(int j = 1; j <= k; j++){
        result *= n_d;
        result /= j;
        n_d --;
    }
    
    return result;

    
}



/**
 * Calculate the factorial of n
 *
 * @param n: int to calculate double factorial of n 
 * check base case, else recursively call with previous ints
 * @return n!
 */
int double_factorial(int n) {
    if (n <= 0 || n ==1 )
        return 1;
    else
        return n * double_factorial(n - 2);
}