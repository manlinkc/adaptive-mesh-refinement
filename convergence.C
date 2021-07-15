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

/*
 File Description: Takes into text files and calculates L1error
 Make sure text files containing solution data have no words/comments
 */

int main(void){
    
    //------------------------------------------------------------
    //----------READ IN FILE NAMES CONTAINING SOLUTION DATA-------
    std::string filename1;
    std::string filename2;
    std::ifstream file0("convergenceinputs.txt", std::ifstream::in);
    
    while(file0){
        std::string temp0;
        file0 >> filename1;
        file0 >> temp0;
        file0 >> filename2;
        file0 >> temp0;
        
        // If there is nothing in the file then break
        if(!file0)
        {
          break;
        }
        
    }
    
    file0.close();
    
    
    //std::cout <<"filename1 = " << filename1 <<std::endl;
    //std::cout <<"filename2 = " << filename2 <<std::endl;
    
    std::vector<double> rhoexact;
    std::vector<double> rhoapprox;
    
    //------------------------------------------------------------
    //----------READ IN EXACT SOLUTION----------------------------
    // Read in file
    std::ifstream file1(filename1, std::ifstream::in);
    
    // While the file is open
    while(file1){
        
        std::string temp1,temp2,temp3,temp4;
        double exact_val;
        file1 >> temp1;
        file1 >> exact_val; // only read in density values
        file1 >> temp2;
        file1 >> temp3;
        file1 >> temp4;
        
        // If there is nothing in the file then break
        if(!file1)
        {
          break;
        }
        
        // First string corresponds to firstname, add to the vector
        rhoexact.push_back(exact_val);
    }
    
    file1.close();
    
    //------------------------------------------------------------
    //----------READ IN APPROX SOLUTION----------------------------
    
    // Read in file
    std::ifstream file2(filename2, std::ifstream::in);
    
    // While the file is open
    while(file2){
        
        std::string temp5;
        double approx_val;
        
        file2 >> temp5;
        file2 >> approx_val;

        // If there is nothing in the file then break
        if(!file2)
        {
          break;
        }
        
        // First string corresponds to firstname, add to the vector
        rhoapprox.push_back(approx_val);
    }
    file2.close();
    
    //------------------------------------------------------------
    //----------CALCULATE L1 ERROR--------------------------------
    
    int ncells = rhoexact.size();
    double x0 = 0.0;
    double x1 = 1.0;
    double dx = (x1 - x0)/ncells;
    double L1error = 0.0;
    
    for (int i=0; i<ncells; i++){
        
        L1error += fabs(rhoexact[i] - rhoapprox[i]);
        /*
        std::cout << "rhoexact = " << rhoexact[i] << std::endl;
        std::cout << "rhoapprox = " << rhoapprox[i] << std::endl;
        std::cout << "L1error = " << fabs(rhoexact[i] - rhoapprox[i]) << std::endl;
        std::cout << "Total error = " << L1error <<std::endl;
         */
    }
    
    L1error = dx*L1error;
    //std::cout << "dx = " << dx << std::endl;
    
    std::cout << "Resolution  = " << ncells <<std::endl;
    std::cout << "L1error = " << L1error << std::endl;
    
    
    
return 0;
}


