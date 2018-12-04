#pragma once
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <eigen3/Eigen/Dense>
#include <vtkUnsignedCharArray.h>

class Beam{
    public:
    Beam();
    ~Beam();

    protected:    
    vtkSmartPointer<vtkUnstructuredGrid> m_data;
    vtkSmartPointer<vtkUnstructuredGrid> m_giData;

    
    vtkSmartPointer<vtkActor> m_actor;
    vtkSmartPointer<vtkActor> m_giActor;


    vtkSmartPointer<vtkUnsignedCharArray> m_vertexColors;

    //Boundary
    double m_timeStep = 0.001;
    double m_gravity = 0.0;
    double m_mass = 1.0;

    double m_youngsModulus = 30000;
    double m_poissonsRatio = 0.3;
    Eigen::Matrix3d m_Identity;
    Eigen::MatrixXd m_E;

    //Force
    std::vector<Eigen::Vector3d> m_force;
    std::vector<Eigen::Vector3d> m_velocity;
    
    ///For Masehlsess
    Eigen::Vector3d m_iCenterOfMass;
    std::vector<Eigen::Vector3d> m_qi;


    //Original Inverse Matrix for All Tetra
    std::vector<Eigen::Matrix3d> m_orgInverseMatrix;

    //Seelected Point, Boundary Condition
    int m_selectedIdx = -1;


    

    protected:
    void Initialize();
    void InitializeSystem();

    void UpdateForce();

    public:
    vtkSmartPointer<vtkActor> GetActor(){return m_actor;}
    vtkSmartPointer<vtkActor> GetGiActor(){return m_giActor;}

    void ComputeFEM();
    void ComputeMesheless();
    void Update();

    void SetPointPosition(int ID, double x, double y, double z);
};