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
	Eigen::MatrixXd data_;
	std::vector<std::string> paramList_;
	std::vector<std::vector<std::string> > dataList_;
	std::vector<std::vector<double>> list;
public:
	CSVReader(std::string filename, std::string delm = ","):fileName_(filename), delimeter_(delm){
		std::ifstream file(fileName_.c_str()); 
		std::string line = "";
		while (getline(file, line)){
			std::vector<std::string> vec;
			boost::algorithm::split(vec, line, boost::is_any_of(delimeter_));
			if(vec[0]!="#") dataList_.push_back(vec);
			else paramList_ = vec;
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
		data_ = data;
	}

	Eigen::MatrixXd getData(){
		return data_;
	}

	Eigen::VectorXd get1Ddata(int index,Eigen::MatrixXd data_){
	Eigen::VectorXd data(list.size());
	for (int row = 0; row < data.innerSize(); ++row){   		
        	data(row) = data_(row,index) ;		
		}
	return data;
	}

	std::vector<std::string> getParamList(){
		return paramList_;
	}

};

class Dataset{
private:
	/// @brief approximate length.
	int reqTrajLen_;
	/// @brief column number of desired trajectory in each csv file.
	int index_;
	/// @brief number of demonstrations in dataset
	int demoNum_;
	/// @brief number of joints in demonstration.
	int jointNum_;
	/// @brief file names with trajectories for each demonstration.
	std::vector<std::string> fileList_;
	/// @brief Eigen::Matrix of trajectories.
	std::vector<Eigen::MatrixXd> data_;

public:
	/** \brief				 handle data from multiple demonstrations
	*	\param	reqTrajLen   required length of trajectory
	*	\param	index    	 column of required trajectory in each csv file
	*	\param	fileList  	 vector of .csv file names	
	*/
	Dataset(int index,std::vector<std::string> fileList):index_(index),fileList_(fileList){
		demoNum_ = fileList_.size();
		std::vector<int> trjLens;
		std::vector<Eigen::MatrixXd> demoData;
		for(int i=0;i<demoNum_;++i){
			CSVReader reader(fileList_[i]);
			demoData.push_back(reader.getData());
			trjLens.push_back(demoData[i].rows());
			}
		jointNum_ = demoData[0].cols();
		int minLen = *std::min_element(trjLens.begin(), trjLens.end());
		reqTrajLen_=minLen;
		double delta;
		Eigen::MatrixXd data1(reqTrajLen_,jointNum_);
		for(int k=0; k <demoNum_;k++)
		{
			for(int col=0;col<demoData[k].cols();++col){			

				delta = (demoData[k].rows() - 1)*1.0/(reqTrajLen_-1);
				for(int row=0;row<reqTrajLen_;++row){
					data1(row,col)=demoData[k]((int)(row*delta),col);
				}
			}
			data_.push_back(data1);
		}
	}
	
	/** \brief	get required data
	*/
	std::vector<Eigen::MatrixXd> getData(){
		return data_;
		}
};

#endif

