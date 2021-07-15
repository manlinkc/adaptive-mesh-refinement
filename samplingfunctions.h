#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>

// Declare functions used in main code
std::array<double,3> getSample(const double &pstar,const double &ustar,const double &S,const double &gamma,std::array<double,3> &WL, std::array<double,3> &WR);
void ReadParams(int &torotest, double &rho_L, double &u_L, double &p_L, double &rho_R, double &u_R, double &p_R, double &pstar, double &ustar, double &gamma, double&t, int &ncells, std::string &filename2);

