/**
 * @file ProMP.cpp
 * @author Mohamed Al-Khulaqui (mohamed_alk@outlook.com)
 * @brief 
 * @version 0.1
 * @date 2023-01-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */


#include "ProMP/ProMP.hpp"

ProMP_ns::ProMP::ProMP(ProMP_Params& params)
{
    // Unpack params
    demoTraj_ = params.demoData;
    bfNum_ = params.bfNum;
    bfStd_ = params.bfStd;
    phaseRate_ = params.phaseRate;
    savePath_ = params.savePath;
    
    // Define Variables
    trajLen_ = demoTraj_[0].rows();
    jointNum_ = demoTraj_[0].cols();
    demoNum_ = demoTraj_.size();
    
    h_ = sqrt(2*M_PI)*bfStd_;
    bf_c_ = linspace(0,1,bfNum_);
    
    // Combute phase and phi
    computePhase(phaseRate_);
    phi_ = generateBF(z_,trajLen_);
    PHI_ = computeDiagBlock(phi_);
    
    // Learn from Demonstration
    learnMP(demoTraj_);
    
    // Save 
    saveMP();
}

ProMP_ns::ProMP::ProMP(std::string MPfilepath)
{
    CSVReader reader(MPfilepath);
    Eigen::MatrixXd MPdata = reader.getData();
    std::vector<std::string> param = reader.getParamList();
    bfNum_ = std::stod(param[1]);
    jointNum_ = std::stod(param[2]);
    trajLen_ = std::stod(param[3]);    
    
    omMean_ = MPdata.col(0);
    omStd_ = MPdata.block(0,1,bfNum_*jointNum_,bfNum_*jointNum_);
}

void ProMP_ns::ProMP::setStart(double start, double std)
{
    addViaPoint(0.0, start, std);
}

void ProMP_ns::ProMP::setGoal(double goal, double std)
{
    addViaPoint(1.0, goal, std);
}

void ProMP_ns::ProMP::addViaPoint(double t, double pointPos, double std)
{
    Eigen::Vector3d viaPoint;
    viaPoint(0) = t;
    viaPoint(1) = pointPos;
    viaPoint(2) = std;
    viaPoints_.push_back(viaPoint);
}

void ProMP_ns::ProMP::computePhase(double phaseSpeed)
{
    int timePoints = (int)(trajLen_/phaseSpeed);
    double dt = 1.0/trajLen_;
    Phase phase(dt,phaseSpeed,timePoints);
    z_ = phase.getPhase();
}

Eigen::MatrixXd ProMP_ns::ProMP::computeDiagBlock(Eigen::MatrixXd phi)
{
    int r = phi.rows();
    int c = phi.cols();
    Eigen::MatrixXd block(r*jointNum_,c*jointNum_);
    block.setZero();
	int j =0;
	for (int i = 0; i < jointNum_; i++)
	{
		block.block(i*r,j*c,r,c) = phi;
		j++;
	}
    return block;
}

Eigen::MatrixXd ProMP_ns::ProMP::gaussianBasis(Eigen::MatrixXd &z_c)
{
    assert(bfNum_==z_c.rows());
    assert(trajLen_==z_c.cols());
    
    Eigen::MatrixXd phi(bfNum_,trajLen_);
    for (int i = 0; i < bfNum_; i++)
    {
        for (int j = 0; j < trajLen_; j++)
        {
            auto x = z_c(i,j);
            phi(i,j) = std::exp(-0.5 * pow(x/bfStd_,2))/h_;
        }
    }
    return phi;
}

Eigen::MatrixXd ProMP_ns::ProMP::vonMisesBasis(Eigen::MatrixXd &z_c)
{
    assert(bfNum_==z_c.rows());
    assert(trajLen_==z_c.cols());

    Eigen::MatrixXd phi(bfNum_,trajLen_);
    for (int i = 0; i < bfNum_; i++)
    {
        for (int j = 0; j < trajLen_; j++)
        {
            phi(i,j) = std::exp(std::cos((2*M_PI*z_c(i,j))/bfStd_)/h_);
        }
    }
    return phi;
}

Eigen::MatrixXd ProMP_ns::ProMP::generateBF(Eigen::VectorXd &z, int trajLen)
{
    Eigen::MatrixXd zMat = z.replicate(1,bfNum_).transpose();
    Eigen::MatrixXd cMat = bf_c_.replicate(1,trajLen);

    Eigen::MatrixXd z_c = zMat - cMat;
    Eigen::MatrixXd phi0(bfNum_,trajLen_);
    Eigen::VectorXd bfSum(bfNum_);
    
    switch(basisType_)
    {
        case 0: 
            phi0 = gaussianBasis(z_c).transpose();
            bfSum = phi0.rowwise().sum();
            break;
        
        case 1:
            phi0 = vonMisesBasis(z_c).transpose();
            bfSum = phi0.rowwise().sum();
            break;
        
        default:
            phi0 = gaussianBasis(z_c).transpose();
            bfSum = phi0.rowwise().sum();
            break;
        
    }

    Eigen::MatrixXd phi(trajLen_,bfNum_);
    for (int i = 0; i < trajLen_; i++)
    {
        for (int j = 0; j < bfNum_; j++)
        {
            phi(i,j) = phi0(i,j)/bfSum(i);
        }
    }
    return phi.transpose();
}

void ProMP_ns::ProMP::learnMP(std::vector<Eigen::MatrixXd> &dataset)
{
    
    double sig = 3e-7;
    Eigen::MatrixXd om(demoNum_,bfNum_*jointNum_);
    Eigen::MatrixXd noise = sig * Eigen::MatrixXd::Identity(bfNum_*jointNum_,bfNum_*jointNum_);
    
    auto c1 = (PHI_*PHI_.transpose()).inverse() + noise;
    // Estimating Omega
    for (int i = 0; i < dataset.size(); i++)
    {
        Eigen::MatrixXd traj = dataset[i];
        traj.resize(trajLen_*jointNum_,1);
        auto c2 = PHI_*traj;
        om.row(i) = ((c1)*(c2)).transpose();        
    }
    
    // Mean of Omega
    om_ = om;
    omMean_ = om_.colwise().mean();

    // Standard Deviation of Omega
    Eigen::MatrixXd cen = om_.rowwise() - omMean_.transpose();
    omStd_ = (cen.adjoint()*cen) / double(om_.rows()-1);
}

Eigen::MatrixXd ProMP_ns::ProMP::generateTraj(int desiredtrajLen)
{
    double dt = 1.0/trajLen_;
    int timePoints = (int)(trajLen_/phaseRate_);

    auto new_omMean = omMean_;
    auto new_omStd = omStd_;

    Phase newPhase(dt,phaseRate_,timePoints);
    auto newZ_ = newPhase.getPhase();

    auto phi = generateBF(newZ_,trajLen_);
    auto PHI = computeDiagBlock(phi);

    for (int i = 0; i < viaPoints_.size(); i++)
    {
        auto z0 = newPhase.getPhaseFromTime(viaPoints_[i][0]);
        auto phiT = generateBF(z0,1);

        Eigen::VectorXd point(1);
        point(0) = viaPoints_[i][1];
        Eigen::VectorXd std(1);
        std(0) = viaPoints_[i][2];

        auto aux = std + (phiT.transpose()*new_omStd)*phiT;

        new_omMean = new_omMean + (new_omStd*phiT/aux(0)) * (point - phiT.transpose()*new_omMean);
        new_omStd = new_omStd - (new_omStd*phiT/aux(0)) * (phiT.transpose() * new_omStd);
    }

    auto sampleOm = new_omMean;
    Eigen::MatrixXd y_t = PHI.transpose() * sampleOm;
    y_t.resize(trajLen_,jointNum_);
    

    double dt0 = (trajLen_ - 1.0)/(desiredtrajLen);
    Eigen::MatrixXd genTraj(desiredtrajLen,jointNum_);
    for (int j=0; j < jointNum_; j++){
        for (int i = 0; i < desiredtrajLen; i++){
            genTraj(i,j) = y_t((int)(i*dt0),j);
        }    
    }
    return genTraj;
}

void ProMP_ns::ProMP::saveMP()
{
    std::ofstream MPfile;
    MPfile.open(savePath_+"MP.csv");
    MPfile << std::scientific <<"#,"<< bfNum_ <<","<< jointNum_<<"," << trajLen_ << std::endl;
    for(int i=0;i<jointNum_*bfNum_;i++){
        int j=0;
        MPfile << omMean_(i) << ",";
        for(j=0;j< jointNum_*bfNum_-1;j++){
            MPfile << "\t" << omStd_(i,j) << ","; 
        }
        MPfile << "\t" << omStd_(i,j) << std::endl;
    }
    MPfile.close();
}

void ProMP_ns::ProMP::loadMP(const std::string filepath)
{
}

Eigen::VectorXd ProMP_ns::ProMP::linspace(double x0, double xf, int len)
{
    Eigen::VectorXd f(len);
    double dx = (xf-x0)/(len-1);
    for(int i = 0; i < len; i++)
    {
        f(i)=x0+i*dx;
    }
    return f;
}
