/**
 * @file ProMP_Test.cpp
 * @author Mohamed Al-Khulaqui (mohamed_alk@outlook.com)
 * @brief Example of how to use ProMP
 * @version 0.1
 * @date 2023-01-28
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "ProMP/ProMP.hpp"
#include "dataHandle.hpp"

int main(int argc, char** argv)
{

   	std::vector<std::string> fileList;
	std::string filePath = "../testData/";

	fileList.push_back(filePath + "example1.csv");
	fileList.push_back(filePath + "example2.csv");
	fileList.push_back(filePath + "example3.csv");
	fileList.push_back(filePath + "example4.csv");
	
	//Data handle object which takes in index of trajectory in .csv file and also the file list
	Dataset testData(0,fileList);
	
	/// get required trajectory from all demonstrations into a single matrix
	ProMP_ns::ProMP_Params options;
	options.demoData = testData.getData();
	options.savePath = filePath;

	/// initialize promp object with number of basis functions and std as arguments.
	ProMP_ns::ProMP baseMp(options);

	/// generate trajectory of required number of points with the generateTrajectory function.
	Eigen::MatrixXd GenTraj = baseMp.generateTraj(300);

	/// Below code is for writing all the data into an output.csv file for visualization.
	std::ofstream myfile;
	myfile.open (filePath+"output.csv");
	
	for (int i=0;i<GenTraj.rows();++i)
	{
		int j;
		for (j=0; j<GenTraj.cols()-1; j++){	
			myfile <<GenTraj(i,j)<<",";
		}
		myfile <<GenTraj(i,j)<<std::endl;
	}
	myfile.close();
	
	/// Load the MP that was just trained
	ProMP_ns::ProMP loadMP(filePath+"MP.csv");
	return 0;
}