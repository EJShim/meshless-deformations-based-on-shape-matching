#pragma once
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyData.h>

#include <eigen3/Eigen/Dense>

class deformableMesh{
    public:
    deformableMesh(vtkSmartPointer<vtkPolyData> data);
    ~deformableMesh();

    protected:

    //Clusteterd data
    std::vector<std::vector<int>> m_clusterID;

    //Data
    vtkSmartPointer<vtkPolyData> m_data;  
    vtkSmartPointer<vtkPolyData> m_gData;  

    //Actor For Rendering
    vtkSmartPointer<vtkActor> m_actor;
    vtkSmartPointer<vtkActor> m_gActor;

    //Colors, for debug
    vtkSmartPointer<vtkUnsignedCharArray> m_vertexColors;    

    //Boundary
    int m_nCluster = 8;
    double m_timeStep = 0.1;
    double m_gravity = -0.0;
    double m_mass = 1.0;

    //Force
    std::vector<Eigen::Vector3d> m_force;
    std::vector<Eigen::Vector3d> m_velocity;    
    std::vector<Eigen::MatrixXd> m_results;
    std::vector<int> m_avg;
    
    ///For Masehlsess
    // Eigen::Vector3d m_iCenterOfMass;
    // Eigen::Vector3d m_cCenterOfMass;
    // std::vector<Eigen::Vector3d> m_qi;
    // std::vector<Eigen::VectorXd> m_Qi;    
    // Eigen::Matrix3d m_Aqq;
    // Eigen::MatrixXd m_AQQ;



    ////Temp for clustering    
    std::vector<std::vector<Eigen::Vector3d>> m_c_qi;
    std::vector<std::vector<Eigen::VectorXd>> m_c_Qi;
    std::vector<Eigen::Matrix3d> m_c_Aqq;
    std::vector<Eigen::MatrixXd> m_c_AQQ;



    //Seelected Point, Boundary Condition
    int m_selectedIdx = -1;


    

    protected:
    void Initialize(vtkSmartPointer<vtkPolyData> data);
    void InitializeSystem();
    void MakeCluster();

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
    double GetTotalMass();
};