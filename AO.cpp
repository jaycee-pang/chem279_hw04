#include "AO.h"
#include "molecule.h"
AO::AO() {}
AO::AO(std::string& shell, arma::vec& R, arma::vec& alphas, arma::vec& d_coeff,arma::uvec& lmn)
        : shell_type(shell), R_(R), alphas_(alphas), ds_(d_coeff), lmn_(lmn) {
            // normalize and combine these N with the contraction coefficient term 
            for (size_t k = 0; k < alphas_.n_elem; k++) {
                double self_overlap = overlap3d(R_, R_, alphas_(k), alphas_(k), lmn_, lmn_); 
                ds_(k) /= std::sqrt(self_overlap); 
            }

        }
double AO::overlap1d(double xa, double xb, double alpha, double beta, int la, int lb) const {
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
      
overlap3d(const arma::vec& Ra, const arma::vec& Rb, double alphas, double betas, const arma::uvec& la, const arma::uvec& lb)  {
    double Sabx = overlap1d(Ra(0), Rb(0), alphas, betas, la(0), lb(0));
    double Saby = overlap1d(Ra(1), Rb(1), alphas, betas, la(1), lb(1));
    double Sabz = overlap1d(Ra(2), Rb(2), alphas, betas, la(2), lb(2));
    return Sabx * Saby* Sabz; 
}

evaluate_contracted_overlap(const AO& ao1, const AO& ao2) {
    double sum = 0.0; 
    arma::vec alphas = ao1.alphas();
    arma::vec betas = ao2.alphas(); 
    arma::uvec la = ao1.lmn();
    arma::uvec lb = ao2.lmn();
    arma::vec da = ao1.ds();
    arma::vec db = ao2.ds(); 
    arma::vec Ra = ao1.R(); 
    arma::vec Rb = ao2.R(); 
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
void AO::print_AO() const {
        // std::cout << atom_.print_atom()<< ", " << shell_type<< ", " << lmn_.print() << std::endl;
        std::cout << shell_type << ", " << lmn_.print() < <std::endl;
        std::cout <<" alphas: " << alphas_.print() << ",\ncontraction coefficents: " << ds_.print()<< std::endl;
}




std::string AO::shell() const {return shell_type;} 
arma::vec AO::R() const{return R_;}
arma::vec AO::alphas() const {return alphas_;};
arma::vec AO::ds() const {return ds_;}; 
arma::uvec AO::lmn() const {return lmn_;}

void getBasis_data(const std::string&atom_type) {
    std::string filename;
    
    
    if (atom_type == "C") {
        filename = "C_STO3G.txt"; 
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        C_alphas = arma::vec(alphas);
        C2s_coeffs = arma::vec(coeffs_2s);
        C2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            C2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
        
    }
    else if(atom_type == "H") {
        filename = "H_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_1s;
        // std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_1s; 
            if (!(iss >> alpha >> coeff_1s)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_1s.push_back(coeff_1s);
            // coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        H_alphas = arma::vec(alphas);
        H1s_coeffs = arma::vec(coeffs_1s);


        file.close();
    }
    else if (atom_type == "O") {
        filename = "O_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        O_alphas = arma::vec(alphas);
        O2s_coeffs = arma::vec(coeffs_2s);
        O2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            O2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else if(atom_type == "N") {
        filename = "N_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        N_alphas = arma::vec(alphas);
        N2s_coeffs = arma::vec(coeffs_2s);
        N2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            N2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else if(atom_type == "F") {
        filename = "F_STO3G.txt";
        std::vector<double> alphas;
        std::vector<double> coeffs_2s;
        std::vector<std::vector<double>> coeffs_p;
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Error opening the file " + filename);
        }
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double alpha, coeff_2s, coeff_px, coeff_py, coeff_pz;
            if (!(iss >> alpha >> coeff_2s >> coeff_px >> coeff_py >> coeff_pz)) {
                std::cerr << "Error reading line: " << line << std::endl;
                continue;
            }
            alphas.push_back(alpha);
            coeffs_2s.push_back(coeff_2s);
            coeffs_p.push_back({coeff_px, coeff_py, coeff_pz});
        }
        N_alphas = arma::vec(alphas);
        N2s_coeffs = arma::vec(coeffs_2s);
        N2p_coeffs = arma::mat(coeffs_p.size(), 3);
        for (size_t i = 0; i < coeffs_p.size(); ++i) {
            N2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
        }

        file.close();
    }
    else {
        throw std::runtime_error("No basis set for this atom type ");
    }

}
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
