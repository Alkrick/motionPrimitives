/**
 * @file dataHandle.hpp
 * @author Mohamed Al-Khulaqui (mohamed_alk@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2023-01-13
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef _DATA_HANDLE_H_
#define _DATA_HANDLE_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>

class CSVReader{
private:
	std::string fileName_;
	std::string delimeter_; 
	std::vector<std::vector<std::string> > dataList_;
	std::vector<std::vector<double>> list;
public:
	CSVReader(std::string filename, std::string delm = ","):fileName_(filename), delimeter_(delm){ }
	Eigen::MatrixXd getData(){
	std::ifstream file(fileName_.c_str()); 
	std::string line = "";
	while (getline(file, line)){
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter_));
		dataList_.push_back(vec);
		}
	file.close();	
	for(std::vector<std::string> vec : dataList_){
		std::vector<double> oneLine;
		for(std::string data : vec){
			oneLine.push_back(std::stod(data));
			}
		list.push_back(oneLine);		
		}
	Eigen::MatrixXd data(list.size(),list[0].size());
	for (int row = 0; row < list.size(); ++row){
   		for (int col = 0; col < list[0].size(); ++col){
        	data(row,col) = list[row][col];
   			}
		}
	return data;
	}
	Eigen::VectorXd get1Ddata(int index,Eigen::MatrixXd data_){
	Eigen::VectorXd data(list.size());
	for (int row = 0; row < data.innerSize(); ++row){   		
        	data(row) = data_(row,index) ;		
		}
	return data;
	}

};

class Dataset{
private:
	/// approximate length
	int reqTrajLen_;
	/// column number of desired trajectory in each csv file.
	int index_;
	/// file names with trajectories for each demonstration.
	std::vector<std::string> fileList_;
	/// std::vector of Eigen::VectorXd of trajectories
	std::vector<Eigen::VectorXd> trajList_;
	/// Eigen::Matrix of trajectories
	Eigen::MatrixXd data_;
public:
	/** \brief				 handle data from multiple demonstrations
	*	\param	reqTrajLen   required length of trajectory
	*	\param	index    	 column of required trajectory in each csv file
	*	\param	fileList  	 vector of .csv file names	
	*/
	Dataset(int index,std::vector<std::string> fileList):index_(index),fileList_(fileList){
		std::vector<int> trjLens;
		for(int i=0;i<fileList_.size();++i){
			CSVReader reader(fileList_[i]);
			Eigen::VectorXd data = reader.get1Ddata(index_,reader.getData());
			trjLens.push_back(data.innerSize());
			trajList_.push_back(data);
			}
		int minLen = *std::min_element(trjLens.begin(), trjLens.end());
		reqTrajLen_=minLen;
		double delta;
		Eigen::MatrixXd data1(fileList_.size(),reqTrajLen_);
		for(int i=0;i<trajList_.size();++i){			
			delta = (trajList_[i].innerSize() - 1)*1.0/(reqTrajLen_-1);
			for(int j=0;j<reqTrajLen_;++j){
				data1(i,j)=trajList_[i]((int)(j*delta));
				}
			}
		data_= data1;
		}
	/** \brief	get required data
	*/
	Eigen::MatrixXd getData(){
		return data_;
		}
};

#endif

