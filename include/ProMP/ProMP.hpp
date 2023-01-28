/**
 * @file ProMP.hpp
 * @author Mohamed Al-Khulaqui (mohamed_alk@outlook.com)
 * @brief  My implementation of Probablistic Motion Primitives (ProMP) based on 2013 Paper by Alexandros Paraschos et. al
 * @version 0.1
 * @date 2022-12-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef _PROMP_HPP_
#define _PROMP_HPP_

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "dataHandle.hpp"

#define debugging_enabled 1

#define DEBUG(n,x) do { \
  if (debugging_enabled) { std::cerr << n << x << std::endl; } \
} while (0)

const std::string dataPath = "/home/xuande/My/Workspaces/catkin_ws/src/trajectory_generation/data/";
namespace ProMP_ns{
class ProMP{
public:
    /**
     * @brief constructor
     * @param demoData      demo data
     * @param bfNum         number of basis functions
     * @param bfStd         standard deviation
     */
    ProMP(Eigen::MatrixXd& demoData, int bfNum, double bfStd, double phaseRate);

    /**
     * @brief Load an already trained ProMP from file 
     * 
     * @param MPfilepath trained MP filepath
     */
    ProMP(std::string MPfilepath);

    /**
     * @brief Set the desired start value (time is set at 0)
     * 
     * @param start value at start
     * @param std   standard deviation
     */
    void setStart(double start, double std);

    /**
     * @brief Set the desired goal value (time is set at 1)
     * 
     * @param goal value at end
     * @param std  standard deviation
     */
    void setGoal(double goal, double std);

    /**
     * @brief Set the save Filepath for storing trained MP 
     * 
     * @param filepath 
     */
    void setFilepath(std::string filepath);

    /**
     * @brief Add a via point to the trajectory
     * 
     * @param t         time of via point (0~1)
     * @param pointPos  value at via point 
     * @param std       standard deviation
     */
    void addViaPoint(double t, double pointPos, double std);

    /**
     * @brief Generates a trajectory from MP
     * 
     * @param trajLen trajectory length
     * @return Eigen::VectorXd trajectory
     */
    Eigen::VectorXd generateTraj(int trajLen);

private:
    // Number of Basis Functions
    int bfNum_;

    // Standard deviation of basis functions
    double bfStd_;

    // Width of basis function
    double h_;

    // Number of demos
    int demoNum_;

    // Points in trajectory
    int trajLen_;

    // Phase speed < 1 for slow, > 1 for faster
    int phaseRate_;

    // Basis type 0 -> Gaussian, 1 -> Von-Mises,
    const int basisType_ = 0;

    // Start value of demo trajectories
    double startPoint_;

    // End value of demo tajectories
    double endPoint_;

    // Phase vector
    Eigen::VectorXd z_;

    // Demo trajectories
    Eigen::MatrixXd demoTraj_;

    // C vector of Basis functions
    Eigen::VectorXd bf_c_;

    // Basis functions
    Eigen::MatrixXd phi_;

    // Weight Vector Omega
    Eigen::MatrixXd om_;

    // Mean of Omega
    Eigen::VectorXd omMean_;

    // Standard deviation of Omega
    Eigen::MatrixXd omStd_;

    // Mean of trajectories
    Eigen::VectorXd meanTraj_;

    // Via points of Trajectory
    std::vector<Eigen::Vector3d> viaPoints_;

    // Store MP file path
    std::string savePath_;

    /**
     * @brief Learn average weights of basis functions from demonstration trajectories 
     * 
     * @param dataset Demonstration trajectories
     */
    void learnMP(Eigen::MatrixXd& dataset);

    /** @brief	Computes phase
	*	@param phaseSpeed PhaseSpeed < 1 slower than average trajectory 
	*			          PhaseSpeed > 1 faster than average trajectory
	*/
	void computePhase(double phaseSpeed);
    /**
     * @brief Caluclates the Gaussian basis function
     * 
     * @param z_c Difference between Phase variable z and parameter c 
     * @return Eigen::MatrixXd Matrix of Gaussian basis function values
     */
    Eigen::MatrixXd gaussianBasis(Eigen::MatrixXd& z_c);

    /**
     * @brief Caluclates the Von-Mises basis function
     * 
     * @param z_c Difference between Phase variable z and parameter c 
     * @return Eigen::MatrixXd Matrix of Von-Mises basis function values
     */
    Eigen::MatrixXd vonMisesBasis(Eigen::MatrixXd& z_c);

    /**
     * @brief Generate basis functions
     * 
     * @param z Phase variable
     * @param trajLen Required length of trajectory
     * @return Eigen::MatrixXd Matrix of basis functions
     */
    Eigen::MatrixXd generateBF(Eigen::VectorXd& z, int trajLen);

    /**
     * @brief Store the trained MP parameters in file
     * 
     */
    void saveMP();
    /**
     * @brief load MP parameters from file
     * 
     * @param filepath MP param file path
     */
    void loadMP(const std::string filepath);

    /**
     * @brief create linearly spaced vector 
     * 
     * @param x0 start
     * @param xf end
     * @param len num of points
     * @return Eigen::VectorXd vector of linearly spaced points
     */
    Eigen::VectorXd linspace(double x0,double xf,int len);
};

class Phase{
public:
    /**
     * @brief Construct a new Phase object
     * 
     * @param dt            time step
     * @param phaseRate     <1 for slow, >1 for fast
     * @param numTimePoints Num of points in traj
     */
    Phase(double dt, double phaseRate, int numTimePoints)
    :dt_(dt),phaseRate_(phaseRate),numTimePoints_(numTimePoints)
    {
        double phaseStart = -dt;
        phaseEnd_ = phaseStart + numTimePoints_ * dt_;

        Eigen::VectorXd dz = phaseRate_*Eigen::VectorXd::Ones(numTimePoints_);
        Eigen::VectorXd z(dz.size());
		double acc = 0;
		for(int i = 0; i < dz.innerSize(); i++){
	 		acc += dz(i);
	 		z(i) = acc;
		}
        z_ = z*dt;
    }

    /**
     * @brief Get a new phase From different time steps object
     * 
     * @param numTimeSteps new time step num
     * @return Eigen::VectorXd new phase
     */
    Eigen::VectorXd getPhaseFromTime(int numTimePoints){
        Eigen::VectorXd newPhase(1);
		newPhase(0) = numTimePoints/phaseEnd_;
		return newPhase;
    }

    /**
     * @brief Get the Phase object
     * 
     * @return Eigen::VectorXd phase vector
     */
    Eigen::VectorXd getPhase()
    {
        return z_;
    }
private:
    Eigen::VectorXd z_;
    double dt_;
    double phaseRate_;
    double phaseEnd_;
    int numTimePoints_;
};
}
#endif