#include<iostream> 
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<array>
#include<algorithm>
#include<functional>
#include<iostream>
#include<fstream>
#include "samplingfunctions.h"

/*
 File Description: Sampling function to get exact riemann solution to Toro's 5 tests
 Date: Sunday 21st February 2021
 */

int main(void){
    double rho_L, u_L, p_L, rho_R, u_R, p_R, pstar, ustar, gamma, x0, x1, t, dx, x, S;
    int ncells, torotest;
    std::string filename2;
    
    //Read in from inputs.txt
    ReadParams(torotest, rho_L, u_L, p_L, rho_R, u_R, p_R, pstar, ustar, gamma, t, ncells, filename2);

    //Construct domain
    x0 = -0.5;
    x1 = 0.5;
    dx = (x1-x0)/ncells;
    
    //Construct left and right states
    std::array<double,3> WL;
    std::array<double,3> WR;
    WL[0] = rho_L;
    WL[1] = u_L;
    WL[2] = p_L;
    WR[0] = rho_R;
    WR[1] = u_R;
    WR[2] = p_R;
        

    //exact_values.resize(ncells);
    
    std::vector<double> xvalues;
    // Read in approx solution x values
    std::ifstream file2(filename2, std::ifstream::in);
    
    // While the file is open
    while(file2){
        
        double approx_xval;
        std::string temp5;
        
        file2 >> approx_xval;
        file2 >> temp5;

        // If there is nothing in the file then break
        if(!file2)
        {
          break;
        }
        
        // First string corresponds to firstname, add to the vector
        xvalues.push_back(approx_xval);
    }
    file2.close();
    
    int ncells2 = xvalues.size();
    // Vector containing exact rho, u, p for all x
    std::vector<std::array<double, 3> > exact_values(ncells2);
    
    for (int i = 0; i < ncells2; i++){
        //x = x0 + i*dx;
        x = xvalues[i] - 0.5;
        S = x/t;
        
        //std::cout << "S = " << S<< std::endl;
        std::array<double,3> sampled_values = getSample(pstar,ustar,S,gamma, WL, WR);
        
        exact_values[i][0] = sampled_values[0];
        exact_values[i][1] = sampled_values[1];
        exact_values[i][2] = sampled_values[2];

    }

    //--------------------------OUTPUT FINAL DATA-------------------
    //--------------------------------------------------------------
    std::string resolution = std::to_string(ncells2);
    std::string test = std::to_string(torotest);
    
    // Output the final data for plotting
    std::ofstream output("exact" +  test + "_" + resolution +".txt");
    
    // Don't include extra ghost cell into final data output
    for (int i = 0; i < ncells2; i++)
    {
            double x = xvalues[i];
            double internalenergy = exact_values[i][2]/((gamma - 1.0)*exact_values[i][0]);
    
            output << x << " " << exact_values[i][0] << " " << exact_values[i][1] << " " << exact_values[i][2] << " " << internalenergy << std::endl;
    }
    
    //Close output file
    output.close();
    
  
	
return 0;
}


