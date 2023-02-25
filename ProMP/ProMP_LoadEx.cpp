/**
 * @file ProMP_LoadEx.cpp
 * @author Mohamed Al-Khulaqui (mohamed_alk@outlook.com)
 * @brief Example of how to load a MP
 * @version 0.1
 * @date 2023-02-24
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "ProMP.hpp"

int main(int argc,char** argv){
    if (argc<2) {
        std::cerr<<"Missing file directory!\n";
        std::cerr<<"Usage: ProMP_LoadEx PATH_TO_MP_FILE \n";
        return -1;
    }
    std::string filePath = argv[1];
    
    // MP's Data is loaded and can be used
    ProMP_ns::ProMP MP1(filePath);
    
    // Get current parameters
    ProMP_ns::ProMP_Params params = MP1.getParams();
    
    // Printout the current parameters
    std::cout<<params<<std::endl;
    
    // Edit parameters
    params.bfNum = 35;
    params.bfStd = 0.0286;
    //params.demoData = Eigen::MatrixXd::Zero();
    params.phaseRate = 1;
    params.savePath = "";
    
    std::cout<<params<<std::endl;    

    return 0;
}