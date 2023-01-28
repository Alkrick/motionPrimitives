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

ProMP_ns::ProMP::ProMP(Eigen::MatrixXd &demoData, int bfNum, double bfStd, double phaseRate)
:demoTraj_(demoData),bfNum_(bfNum),bfStd_(bfStd),phaseRate_(phaseRate)
{
    // Define Variables
    demoNum_ = demoData.rows();
    trajLen_ = demoData.cols();
    DEBUG("demoNum = ",demoNum_);
    DEBUG("trajLen = ",trajLen_);
    h_ = sqrt(2*M_PI)*bfStd_;
    bf_c_ = linspace(0,1,bfNum_);

    // Combute phase and phi
    computePhase(phaseRate_);
    phi_ = generateBF(z_,trajLen_);
    DEBUG("phi=\n",phi_);
    // Learn from Demonstration
    learnMP(demoTraj_);
}

ProMP_ns::ProMP::ProMP(std::string MPfilepath)
{
}

void ProMP_ns::ProMP::setStart(double start, double std)
{
    addViaPoint(0.0, start, std);
}

void ProMP_ns::ProMP::setGoal(double goal, double std)
{
    addViaPoint(1.0, goal, std);
}

void ProMP_ns::ProMP::setFilepath(std::string filepath)
{
}

void ProMP_ns::ProMP::addViaPoint(double t, double pointPos, double std)
{
    Eigen::Vector3d viaPoint;
    viaPoint(0) = t;
    viaPoint(1) = pointPos;
    viaPoint(2) = std;
    viaPoints_.push_back(viaPoint);
}

void ProMP_ns::ProMP::learnMP(Eigen::MatrixXd &dataset)
{
    double sig = 0.0000000000001;
    Eigen::MatrixXd om(demoNum_,bfNum_);
    Eigen::MatrixXd noise = sig * Eigen::MatrixXd::Identity(bfNum_,bfNum_);
    
    // Estimating Omega
    for (int i = 0; i < dataset.rows(); i++)
    {
        auto c1 = (phi_*phi_.transpose()).inverse() + noise;
        auto c2 = phi_*dataset.row(i).transpose();
        om.row(i) = ((c1)*(c2)).transpose();
    }
    
    // Mean of Omega
    om_ = om;
    omMean_ = om_.colwise().mean();
    DEBUG("Omega= \n",omMean_);
    // Standard Deviation of Omega
    Eigen::MatrixXd cen = om_.rowwise() - om_.colwise().mean();
    omStd_ = (cen.adjoint()*cen) / double(om_.rows()-1);
    DEBUG("STD= \n",omStd_);
    // Mean of Trajectory
    meanTraj_ = dataset.colwise().mean();
    startPoint_ = meanTraj_(0);
    endPoint_ = meanTraj_(trajLen_ - 1);

}

void ProMP_ns::ProMP::computePhase(double phaseSpeed)
{
    int timePoints = (int)(trajLen_/phaseSpeed);
    double dt = 1.0/trajLen_;
    Phase phase(dt,phaseSpeed,timePoints);
    z_ = phase.getPhase();
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
            phi(i,j) = std::exp(-0.5 * pow(x/bfStd_,2))/(sqrt(2*3.14)*bfStd_);
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
    DEBUG("zMAT=\n",zMat);
    DEBUG("cMAT=\n",cMat);
    Eigen::MatrixXd z_c = zMat - cMat;
    Eigen::MatrixXd phi0(bfNum_,trajLen_);
    Eigen::VectorXd bfSum(bfNum_);
    
    switch(basisType_)
    {
        case 0: 
            DEBUG("case",0);
            phi0 = gaussianBasis(z_c).transpose();
            DEBUG("sum",0);
            bfSum = phi0.rowwise().sum();
            break;
        
        case 1:
            DEBUG("what",0);
            phi0 = vonMisesBasis(z_c).transpose();
            bfSum = phi0.rowwise().sum();
            break;
        
        default:
            DEBUG("huh",0);
            phi0 = gaussianBasis(z_c).transpose();
            bfSum = phi0.rowwise().sum();
            break;
        
    }

    DEBUG("PHI0=\n",phi0);
    
    Eigen::MatrixXd phi(trajLen_,bfNum_);
    for (int i = 0; i < trajLen_; i++)
    {
        for (int j = 0; j < bfNum_; j++)
        {
            phi(i,j) = phi0(i,j)/bfSum(i);
        }
    }
    DEBUG("PHI=\n",phi.transpose());
    return phi.transpose();
}

Eigen::VectorXd ProMP_ns::ProMP::generateTraj(int desiredtrajLen)
{
    double dt = 1.0/trajLen_;
    int timePoints = (int)(trajLen_/phaseRate_);

    auto new_omMean = omMean_;
    auto new_omStd = omStd_;

    Phase newPhase(dt,phaseRate_,timePoints);
    auto newZ_ = newPhase.getPhase();

    auto phi = generateBF(newZ_,trajLen_);
    DEBUG("gen_phi= \n",phi);

    for (int i = 0; i < viaPoints_.size(); i++)
    {
        DEBUG("wah",0);
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
    auto y_t = phi.transpose() * sampleOm;
    DEBUG("Omega Mean = \n",sampleOm);
    DEBUG("PHI * Omega = \n",y_t);

    double dt0 = (trajLen_ - 1.0)/(desiredtrajLen);
    Eigen::VectorXd genTraj(desiredtrajLen);
    for (int i = 0; i < desiredtrajLen; i++)
    {
        genTraj(i) = y_t((int)(i*dt0));
    }
    return genTraj;
}

void ProMP_ns::ProMP::saveMP()
{
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
