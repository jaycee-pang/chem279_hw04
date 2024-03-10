#include "gaussians.h"
class STO3GBasis {
    private:
        static arma::vec C_alphas, H_alphas, O_alphas, N_alphas, F_alphas; 
        static arma::mat C2s_coeffs, C2p_coeffs, H1s_coeffs,O2s_coeffs, O2p_coeffs,
                            N2s_coeffs, N2p_coeffs, F2s_coeffs, F2p_coeffs; 

    public: 
        static void getBasis(const std::string& atom_type) {
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
                // C2p_coeffs = arma::mat(coeffs_p.size(), 3);
                // for (size_t i = 0; i < coeffs_p.size(); ++i) {
                //     C2p_coeffs.row(i) = arma::vec(coeffs_p[i]);
                // }

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
                throw std::cerr("No basis set for this atom type "+atom_type);
            }
            
            
            

        }
}
std::vector<std::vector<int>> get_lmn(int L) {
    std::vector<std::vector<int>> combos; 
    for (int l = L; l >= 0; l--) {
        for (int m = L - l; m >= 0; m--) {
            int n = L - l - m;
            if (n >= 0) {
                combos.push_back({l, m, n});
            }
        }
    }
    return combos; 
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
PrimitiveGaussian(arma::vec& center, arma::vec& lmn, int atom) {}

// PrimitiveGaussian::PrimitiveGaussian(): R_({0.0, 0.0, 0.0}), N_(0.0), alpha_(0.0), l_(0), lmn_{0, 0,0} {}
// PrimitiveGaussian::PrimitiveGaussian(double x, double y, double z, double alpha, int l,std::vector<int> lmn) : R_({x, y, z}), alpha_(alpha), l_(l), N_(1.0), lmn_(lmn) {}
// PrimitiveGaussian::PrimitiveGaussian(double x, double y, double z, int l, double alpha): R_({x, y, z}), l_(l), alpha_(alpha), N_(1.0){} 
double PrimitiveGaussian::overlap1d(double xa, double xb, double alpha, double beta,int la, int lb) const {
    double dx = (xa - xb);
    
    double exponential_term = exp(-(alpha*beta*dx*dx)/(alpha+ beta))* sqrt(M_PI / (alpha+ beta));
    double Rp = (alpha * xa + beta * xb)/ (alpha + beta); // new center 
    double result = 0.0; 
    
    for (int i = 0; i <= la; i++) {
        for (int j = 0; j <= lb; j++) {
            if ((i +j)%2 == 1) 
                continue; 
            double double_fact = double_factorial(i +j-1);
            double coeffs = choose(la, i) * choose(lb, j);
            double num = std::pow(Rp-xa, (la -i)) * std::pow(Rp -xb, (lb-j));
            double denom = std::pow(2*(alpha+beta), (double(i+j)/2));
            double term = coeffs * double_fact * num / denom; 
            result += term;
            
        }
    }
    result *= exponential_term; //*result; // *= exponential_term; 
    return result; 
}  
int PrimitiveGaussian::dim_func() const { // get number of basis funcitons for this primitive function
    return (l_+1)*(l_+2)/2;
}



double PrimitiveGaussian::normalize_primitive() {
    double integral = overlap1d(R_[0], R_[0], alpha_, alpha_, lmn_[0], lmn_[0]) *
                    overlap1d(R_[1], R_[1], alpha_, alpha_, lmn_[1], lmn_[1]) *
                    overlap1d(R_[2], R_[2], alpha_, alpha_,lmn_[2], lmn_[2]);
    
    N_ = 1.0 / std::sqrt(integral);
    return N_;

}
double PrimitiveGaussian::overlap3d(const PrimitiveGaussian& other) const {
    double Sabx = overlap1d(R_[0], other.R_[0], alpha_, other.alpha_, lmn_[0], other.lmn_[0]);
    double Saby = overlap1d(R_[1], other.R_[1], alpha_, other.alpha_, lmn_[1], other.lmn_[1]);
    double Sabz = overlap1d(R_[2], other.R_[2], alpha_, other.alpha_, lmn_[2], other.lmn_[2]);
    return Sabx * Saby* Sabz; 
}

std::vector<double> PrimitiveGaussian::R() const {return R_;}
double PrimitiveGaussian::alpha() const {return alpha_;}
double PrimitiveGaussian::N() const {return N_;}
int PrimitiveGaussian::l() const {return l_;}
void PrimitiveGaussian::set_alpha(double new_alpha) {alpha_ = new_alpha;}
void PrimitiveGaussian::set_l(int l) {l_ = l; }
void PrimitiveGaussian::set_R(std::vector<double> R) {R_ = R;}
void PrimitiveGaussian::set_lmn(std::vector<int> lmn) {lmn_ = lmn;}
std::vector<int> PrimitiveGaussian::lmn() const {return lmn_;}
void PrimitiveGaussian::set_N(double N) {N_ = N;}
std::ostream& operator<<(std::ostream& os, const PrimitiveGaussian& p) {
    std::vector<double> center = p.R(); 
    std::vector<int> lmn = p.lmn(); 
    
    os << "Center: " << "(" << center[0] << ", " << center[1] << ", " << center[2] << ")\t" 
    << "alpha: " << p.alpha() << "\t" << "L: " << p.l() <<"\t" 
    << "shell: (" << lmn[0] << "," << lmn[1] << "," << lmn[2] << ")" << "\t" <<"N: " << p.N()<< std::endl;

    
    return os;
}




ContractedGaussian::ContractedGaussian() {}
ContractedGaussian::ContractedGaussian() {}
void generateAOs(vector, int atomic_num) {
    std::vector<PrimitiveGaussian> primitives; 
}
// ContractedGaussian::ContractedGaussian(std::vector<int> lmn, std::string shell): lmn_(lmn), shell_(shell){}
// void ContractedGaussian::addPrimitives(const std::vector<PrimitiveGaussian>& p, const std::vector<double>& contraction_coeffs) {
//     primitives_ = p;
//     ds_ = contraction_coeffs;

// }

std::vector<double> ContractedGaussian::getN() {
    std::vector<double> norm_constants; 
    for (PrimitiveGaussian & p: primitives_) {
        double N = p.normalize_primitive();
        p.set_N(N);
        norm_constants.push_back(N);
    }
    Ns_ = norm_constants; 
    return norm_constants; 
    
}

double ContractedGaussian::contracted_overlap(const ContractedGaussian & other) const {
    double overlap = 0.0;
    for (int i = 0; i < 3; i++) { 
        for (int j = i; j < 3; j++) {
            double primitive_overlap = primitives_[i].overlap3d(other.primitives_[j]);
            overlap += Ns_[i] * other.Ns_[j]*ds_[i] * other.ds_[j] * primitive_overlap;
        }
    }
    return overlap; 

}


std::vector<double> ContractedGaussian::alphas() const {return alphas_;}
std::vector<int> ContractedGaussian::Ls() const {return Ls_;} // l+m+n = L;
std::vector<double> ContractedGaussian::R() const {return R_;}
std::vector<double> ContractedGaussian::ds() const {return ds_;}
std::string ContractedGaussian::shell() const {return shell_;} 
void ContractedGaussian::set_R(const std::vector<double>& center) {R_ = center;}
const std::vector<PrimitiveGaussian>& ContractedGaussian::primitives() const {
        return primitives_;
}


std::ostream& operator<<(std::ostream& os, const ContractedGaussian& c) {
    std::vector<PrimitiveGaussian> primitives = c.primitives();
    std::vector<double> contraction_coeffs = c.ds(); 
    os << "shell type: " << c.shell() << std::endl;
    for (int i = 0; i < primitives.size(); i++) {
        os << "Primitive " << i << ": " << std::endl;
        os << "Primitives" << primitives[i];
        os << "contraction coeff: " <<contraction_coeffs[i] << std::endl;
    }

    return os;
}

// std::vector<std::vector<double>> PrimitiveGaussian::overlap_Sab(const PrimitiveGaussian & other) const {
//     std::vector<std::vector<double>> overlap;
//     int dim1 = dim_func();
//     int dim2 = other.dim_func();
//     std::vector<double> Rb = other.R(); 
//     int lb = other.l(); 
//     int la = l_;
//     double beta = other.alpha();
//     overlap.resize(dim1, std::vector<double>(dim2));
//     int a_index = 0; 
//     for (int la_x = la; la_x >= 0; la_x--) {
//         for (int la_y = la - la_x; la_y >= 0; la_y--) {
//             int b_index = 0;
//             for (int lb_x = lb; lb_x >= 0; lb_x--) {
//                 for (int lb_y = lb - lb_x; lb_y >= 0; lb_y--) {
//                     overlap[a_index][b_index] = overlap_1d(R_[0], Rb[0], alpha_, beta, la_x, lb_x) *
//                                                 overlap_1d(R_[1], Rb[1], alpha_, beta, la_y, lb_y) *
//                                                 overlap_1d(R_[2], Rb[2], alpha_, beta, la - la_x - la_y, lb -lb_x - lb_y); // *
//                                                 // N_ * other.N();
//                     b_index += 1;    
//                 }
//             }
//             a_index += 1;
//         }
//     }  
//     return overlap;
// }