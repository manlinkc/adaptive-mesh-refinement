#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<sstream>
#include "samplingfunctions.h"

std::array<double,3> getSample(const double &pstar,const double &ustar,const double &S,const double &gamma,std::array<double,3> &WL, std::array<double,3> &WR){
    
    double rho_exact, u_exact, p_exact;
    
    if (S <= ustar){
        //Sampling point lies to the left of the contact discontinuity
        
        double rho_L, u_L, p_L, a_L;
        rho_L = WL[0];
        u_L = WL[1];
        p_L = WL[2];
        
        a_L = sqrt((gamma*p_L)/rho_L); //left sound speed
        
        if (pstar <= p_L){
            //Left rarefraction
            
            double SHL = u_L - a_L;
            
            if (S <= SHL){
                // Sampled point is left data state
                rho_exact = rho_L;
                u_exact = u_L;
                p_exact = p_L;
            }

            else{
                double G, astar_L;
                G = (gamma - 1.0)/(2.0*gamma);
                astar_L = a_L*std::pow( (pstar/p_L), G); //pow might not work
                
                double STL = ustar - astar_L;
                
                if(STL <= S){
                    // Sampled point is Star Left State i.e WL star fan
                    // WL star fan = (rho star fan, ustar, pstar)
                    double G1 = 1.0/gamma;
                    
                    rho_exact = rho_L*std::pow( (pstar/p_L), G1);
                    u_exact = ustar;
                    p_exact = pstar;
                }
                else{
                    // Sampled point is inside left fan i.e. WL fan
                    double G2 = 2.0/(gamma + 1.0);
                    double G3 = (gamma - 1.0)/(gamma + 1.0);
                    double G4 = 2.0/(gamma - 1.0);
                    double G5 = 1.0/G4;
                    
                    double terms = G2 + (G3/a_L)*(u_L - S);
                    
                    rho_exact = rho_L*std::pow(terms,G4);
                    u_exact = G2*(a_L + G5*u_L + S);
                    p_exact = p_L*std::pow(terms, (G4*gamma));
                }
            }
        }
        
        else{
                // Left shock
            double G6 = (gamma + 1.0)/(2.0*gamma);
            double G7 = (gamma - 1.0)/(2.0*gamma);
            
            double PML = pstar/p_L;
            double SL = u_L - a_L*sqrt(G6*PML + G7);
        
            
            if (S <= SL){
                // Sampled point is left data state
                rho_exact = rho_L;
                u_exact = u_L;
                p_exact = p_L;
            }
            else{
                // Sampled point is Star Left State
                // WL star shock = (rho star shock, ustar, pstar)
                double G8 = (gamma - 1.0)/(gamma + 1.0);
        
                rho_exact = rho_L*((PML + G8)/(G8*PML + 1.0));
                u_exact = ustar;
                p_exact = pstar;
            }
        }
    }
    
    else{
        // Sampling point lies to the right of the contact discontinuity
        
        double rho_R, u_R, p_R, a_R;
        rho_R = WR[0];
        u_R = WR[1];
        p_R = WR[2];
        
        a_R = sqrt((gamma*p_R)/rho_R); //right sound speed
        
        if (pstar > p_R){
            //Right shock
            
            double PMR = pstar/p_R;
            double G9 = (gamma + 1.0)/(2.0*gamma);
            double G10 = (gamma - 1.0)/(2.0*gamma);
            
            double SR = u_R + a_R*sqrt(G9*PMR + G10);
            
            if (SR <= S){
                // Sampled point is right data state
                
                rho_exact = rho_R;
                u_exact = u_R;
                p_exact = p_R;
            }
            else{
                
                // Sampled point is Star Right state
                // WR star shock = (rho R star shock, ustar, pstar)
                
                double G11 = (gamma - 1.0)/(gamma + 1.0);
                
                rho_exact = rho_R*((PMR + G11)/(G11*PMR + 1.0));
                u_exact = ustar;
                p_exact = pstar;
            }
        }
        else{
            // Right rarefraction
            double SHR = u_R + a_R;
            
            if (SHR <= S){
                
                // Sampled point is right data state
                rho_exact = rho_R;
                u_exact = u_R;
                p_exact = p_R;
            }
            else{
                double G12 = (gamma - 1.0)/(2.0*gamma);
                double astar_R = a_R*std::pow((pstar/p_R), G12);
                double STR = ustar + astar_R;
                
                if(S <= STR){
                    
                    // Sampled point is Star Right state
                    // W R star fan(rho star R fan, ustar, pstar)
                
                    rho_exact = rho_R*std::pow((pstar/p_R), (1.0/gamma));
                    u_exact = ustar;
                    p_exact = pstar;
                }
                else{
                    //Sampled point is inside left fan
                    
                    double G13 = 2.0/(gamma + 1.0);
                    double G14 = (gamma - 1.0)/(gamma + 1.0);
                    double G15 = 2.0/(gamma - 1.0);
                    double G16 = 1.0/G15;
                    
                    double terms2 = G13 - (G14/a_R)*(u_R - S);
                    
                    rho_exact = rho_R*std::pow(terms2, G15);
                    u_exact = G13*(-1.0*a_R + G16*u_R + S);
                    p_exact = p_R*std::pow(terms2, (G15*gamma));
                }
            }
        }
    }
    
    std::array<double, 3> sample_exact;
    sample_exact[0] = rho_exact;
    sample_exact[1] = u_exact;
    sample_exact[2] = p_exact;
    
    return sample_exact;
}

void ReadParams(int &torotest, double &rho_L, double &u_L, double &p_L, double &rho_R, double &u_R, double &p_R, double &pstar, double &ustar, double &gamma, double&t, int &ncells, std::string & filename2){
    
    std::string temp;
    std::ifstream file;

    file.open("inputs.txt");
    
    std::getline(file, temp);
    std::istringstream input0(temp);
    input0 >> torotest;
    
    std::getline(file, temp);
    std::istringstream input1(temp);
    input1 >> rho_L;
    
    std::getline(file, temp);
    std::istringstream input2(temp);
    input2 >> u_L;
    
    std::getline(file, temp);
    std::istringstream input3(temp);
    input3 >> p_L;
    
    std::getline(file, temp);
    std::istringstream input4(temp);
    input4 >> rho_R;
    
    std::getline(file, temp);
    std::istringstream input5(temp);
    input5 >> u_R;
    
    std::getline(file, temp);
    std::istringstream input6(temp);
    input6 >> p_R;
    
    std::getline(file, temp);
    std::istringstream input7(temp);
    input7 >> pstar;
    
    std::getline(file, temp);
    std::istringstream input8(temp);
    input8 >> ustar;
    
    std::getline(file, temp);
    std::istringstream input9(temp);
    input9 >> gamma;
    
    std::getline(file, temp);
    std::istringstream input10(temp);
    input10 >> t;
    
    std::getline(file, temp);
    std::istringstream input11(temp);
    input11 >> ncells;
    
    std::getline(file, temp);
    std::istringstream input12(temp);
    input12 >> filename2;

}

