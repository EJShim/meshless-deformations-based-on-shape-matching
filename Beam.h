#pragma once
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkExtractEdges.h>
#include <eigen3/Eigen/Dense>

class Beam{
    public:
    Beam();
    ~Beam();

    protected:    
    vtkSmartPointer<vtkUnstructuredGrid> m_data;
    vtkSmartPointer<vtkExtractEdges> m_edgeExtractor;
    vtkSmartPointer<vtkActor> m_actor;

    //Boundary
    double m_timeStep = 0.000007;
    double m_gravity = -9.8;
    double m_mass = 1.0;

    double m_youngsModulus = 30000;
    double m_poissonsRatio = 0.3;
    Eigen::Matrix3d m_Identity;
    Eigen::MatrixXd m_E;

    //Force
    std::vector<Eigen::Vector3d> m_force;
    std::vector<Eigen::Vector3d> m_velocity;


    //Original Inverse Matrix for All Tetra
    std::vector<Eigen::Matrix3d> m_orgInverseMatrix;


    

    protected:
    void Initialize();
    void InitializeSystem();

    void UpdateForce();

    public:
    vtkSmartPointer<vtkActor> GetActor(){return m_actor;}
    void ComputeFEM();
    void Update();


};