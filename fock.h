#include <iostream> 
#include <armadillo> 
#include <map> 
std::map<std::string, int> valence_e = {{"H",1}, {"C",4}, {"N", 5},
                                         {"O",6}, {"F", 7}};
class Fock {
    private:
        arma::mat F_; 
    public: 
};