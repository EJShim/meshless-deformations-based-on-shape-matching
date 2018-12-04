#pragma once
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <eigen3/Eigen/Dense>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyData.h>

class Beam{
    public:
    Beam();
    ~Beam();

    protected:

    //Data
    vtkSmartPointer<vtkPolyData> m_data;    

    //Actor For Rendering
    vtkSmartPointer<vtkActor> m_actor;

    //Colors, for debug
    vtkSmartPointer<vtkUnsignedCharArray> m_vertexColors;

    //Boundary
    double m_timeStep = 0.001;
    double m_gravity = 0.0;
    double m_mass = 1.0;

    //Force
    std::vector<Eigen::Vector3d> m_force;
    std::vector<Eigen::Vector3d> m_velocity;
    
    ///For Masehlsess
    Eigen::Vector3d m_iCenterOfMass;
    std::vector<Eigen::Vector3d> m_qi;
    Eigen::Matrix3d m_Aqq;


    //Seelected Point, Boundary Condition
    int m_selectedIdx = -1;


    

    protected:
    void Initialize();
    void InitializeSystem();

    void UpdateForce();

    public:
    vtkSmartPointer<vtkActor> GetActor(){return m_actor;}    

    void ComputeFEM();
    void ComputeMesheless();
    void Update();

    void SetPointPosition(int ID, double x, double y, double z);
};