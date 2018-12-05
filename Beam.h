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
    vtkSmartPointer<vtkPolyData> m_gData;  

    //Actor For Rendering
    vtkSmartPointer<vtkActor> m_actor;
    vtkSmartPointer<vtkActor> m_gActor;

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
    std::vector<Eigen::VectorXd> m_Qi;
    
    Eigen::Matrix3d m_Aqq;
    Eigen::MatrixXd m_AQQ;


    //Seelected Point, Boundary Condition
    int m_selectedIdx = -1;


    

    protected:
    void Initialize();
    void InitializeSystem();

    void UpdateForce();

    public:
    vtkSmartPointer<vtkActor> GetActor(){return m_actor;}
    vtkSmartPointer<vtkActor> GetDebugActor(){return m_gActor;}

    void ComputeFEM();
    void ComputeMesheless();
    void Update();

    void SetPointPosition(int ID, double x, double y, double z);
    void ApplyForce(int ID, double x, double y, double z);

    double* GetCurrentSelectedPosition(int id);
};