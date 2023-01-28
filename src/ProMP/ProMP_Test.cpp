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
	std::string filePath = "/home/xuande/My/Workspaces/catkin_ws/src/trajectory_generation/data/";
	fileList.push_back(filePath + "example1.csv");
	fileList.push_back(filePath + "example2.csv");
	fileList.push_back(filePath + "example3.csv");
	fileList.push_back(filePath + "example4.csv");
   DEBUG("File=",0);
	/**
	* data handle object which takes in index of trajectory in .csv file and also the file list
	*/
	Dataset pelvisData(0,fileList);
   DEBUG("dataset=",0);
	/// get required trajectory from all demonstrations into a single matrix
	Eigen::MatrixXd dataBase1 = pelvisData.getData();	
   DEBUG("dataBase=",0);
	/// initialize promp object with number of basis functions and std as arguments.
	ProMP_ns::ProMP baseMp(dataBase1,5,0.0286, 1.0);
   DEBUG("MP=",0);
	/// generate trajectory of required number of points with the generateTrajectory function.
	Eigen::VectorXd vect = baseMp.generateTraj(100);
   DEBUG("TrajGen= \n",vect);



	/// Below code is for writing all the data into an output.csv file for visualization.
	std::ofstream myfile;
	myfile.open ("./src/trajectory_generation/data/output.csv");
	int j;
	int k;
	int plotlen;
	if (vect.innerSize()>dataBase1.outerSize())
		plotlen =vect.innerSize();
	else
		plotlen =dataBase1.outerSize();
	for (int i=0;i<plotlen;++i){
		if (i<dataBase1.outerSize())
			j=i;
		else 
			j=dataBase1.outerSize()-1;

		if (i<vect.innerSize())
			k=i;
		else 
			k=vect.innerSize()-1;
		myfile << dataBase1(0,j)<<","<<dataBase1(1,j)<<","<< dataBase1(2,j)<<","<< dataBase1(3,j)<<","<<vect(k)<<"\n";
	}
	myfile.close();
   return 0;
}