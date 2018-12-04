#include "Beam.h"
#include <vtkCubeSource.h>
#include <vtkProperty.h>
#include <vtkCubeSource.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkPolyDataMapper.h>
#include <vtkMath.h>
#include <time.h>
#include <vtkPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractEdges.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataNormals.h>

Beam::Beam(){
    m_actor = NULL;
    Initialize();
}

Beam::~Beam(){

}

void Beam::Initialize(){
    // Create a cube.
    vtkSmartPointer<vtkRectilinearGridToTetrahedra> formMesh = vtkSmartPointer<vtkRectilinearGridToTetrahedra>::New();
    formMesh->SetInput(2, 2, 2, 1, 1, 1, 1.0);    
    formMesh->SetTetraPerCellTo6();
    formMesh->Update();
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceExtractor = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceExtractor->SetInputData(formMesh->GetOutput());
    surfaceExtractor->Update();


    // Create a sphere
    vtkSmartPointer<vtkCylinderSource> cylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSource->SetCenter(0.0, 0.0, 0.0);
    cylinderSource->SetRadius(5.0);
    cylinderSource->SetHeight(7.0);
    cylinderSource->SetResolution(10);
    cylinderSource->Update();

      // Create a cube.
    vtkSmartPointer<vtkCubeSource> cubeSource = vtkSmartPointer<vtkCubeSource>::New();
    cubeSource->Update();

    m_data = surfaceExtractor->GetOutput();

    


    // Create a mapper and actor.
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(m_data);
    mapper->Update();

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(mapper);
    m_actor->GetProperty()->SetColor(1, 1, 0);
    m_actor->GetProperty()->SetRepresentationToWireframe();


    //Create Debug Actor
    m_gData = vtkSmartPointer<vtkPolyData>::New();
    m_gData->DeepCopy(m_data);

    vtkSmartPointer<vtkPolyDataMapper> groundMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    groundMapper->SetInputData(m_gData);
    groundMapper->Update();

    m_gActor = vtkSmartPointer<vtkActor>::New();
    m_gActor->SetMapper(groundMapper);
    m_gActor->GetProperty()->SetColor(1, 0, 0);
    m_gActor->GetProperty()->SetRepresentationToPoints();
    m_gActor->GetProperty()->SetPointSize(10);
    
    

    InitializeSystem();
}

void Beam::InitializeSystem(){    

    vtkSmartPointer<vtkPoints> pointSet = m_data->GetPoints();
    


    //Initialize Force and Velocity, color
    m_vertexColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    m_vertexColors->SetNumberOfComponents(3);
    m_vertexColors->SetName("Colors");

    m_iCenterOfMass = Eigen::Vector3d(0, 0, 0);
    int nPoints = m_data->GetNumberOfPoints();
    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d velocity(0, 0, 0);
        Eigen::Vector3d force(0, 0, 0);

        m_velocity.push_back(velocity);
        m_force.push_back(force);

        m_iCenterOfMass += Eigen::Vector3d(m_data->GetPoint(idx));
        m_vertexColors->InsertNextTuple3(1, 1, 0);
    }
    // m_data->GetPointData()->SetScalars(m_vertexColors);

    m_iCenterOfMass /= nPoints;

    m_Aqq = Eigen::MatrixXd::Zero(3,3);
    m_AQQ = Eigen::MatrixXd::Zero(9,9);

    for(int idx = 0 ; idx<nPoints ; idx++){
        Eigen::Vector3d qi = Eigen::Vector3d(m_data->GetPoint(idx)) - m_iCenterOfMass;

        Eigen::VectorXd Qi(9);
        Qi << qi[0], qi[1], qi[2], qi[0]*qi[0], qi[1]*qi[1], qi[2]*qi[2], qi[0]*qi[1], qi[1]*qi[2], qi[2]*qi[0];
        


        m_Aqq += m_mass * qi * qi.transpose();
        m_AQQ += m_mass * Qi * Qi.transpose();
        
        m_qi.push_back(qi);
        m_Qi.push_back(Qi);
    }

    m_Aqq = m_Aqq.inverse();
    m_AQQ = m_AQQ.inverse();
}


void Beam::ComputeMesheless(){

    UpdateForce();


    //Get CenterOfMass    
    Eigen::Vector3d CenterOfMass = Eigen::Vector3d(0, 0, 0);    
    int nPoints = m_data->GetNumberOfPoints();
    for(int idx = 0 ; idx < nPoints ; idx++){
        CenterOfMass += Eigen::Vector3d(m_data->GetPoint(idx));
    }
    CenterOfMass /= nPoints;



    //Calculate Apq
    Eigen::Matrix3d Apq = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd APQ = Eigen::MatrixXd::Zero(3,9);
    
    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d pi = Eigen::Vector3d(m_data->GetPoint(idx)) - CenterOfMass;
        // Eigen::VectorXd Pi(9);
        // Pi << pi[0], pi[1], pi[2], pi[0]*pi[0], pi[1]*pi[1], pi[2]*pi[2], pi[0]*pi[1], pi[1]*pi[2], pi[2]*pi[0];
        
        
        
        Apq += m_mass * pi * m_qi[idx].transpose();
        APQ += m_mass * pi * m_Qi[idx].transpose();
    }

    


    //Inverse SQRT???
    Eigen::Matrix3d S = Apq.transpose() * Apq;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S);
    Eigen::Matrix3d sqrt_S = es.operatorInverseSqrt();

    Eigen::Matrix3d R = Apq * sqrt_S;
    Eigen::MatrixXd RR = Eigen::MatrixXd::Zero(3, 9);
    RR.col(0) = R.col(0);
    RR.col(1) = R.col(1);
    RR.col(2) = R.col(2);


    //Matrix A
    Eigen::Matrix3d A = Apq * m_Aqq;
    Eigen::MatrixXd AA = APQ * m_AQQ;


    double alpha = 0.9;
    double beta = 0.9;


    
    Eigen::Matrix2d factor(2, 2);
    factor  << 0.1, -alpha/m_timeStep,
                m_timeStep, 1-alpha;    
    Eigen::MatrixXd current(2, 3);
    Eigen::MatrixXd ground(2, 3);



    //Update Position
    for(int idx = 1 ; idx < nPoints ; idx++){

        Eigen::Vector3d gi =  ( beta*AA + (1-beta)*RR )*m_qi[idx]+CenterOfMass;
        Eigen::Vector3d xi = Eigen::Vector3d(m_data->GetPoint(idx));

        current.row(0) = m_velocity[idx];
        current.row(1) = xi;


        ground.row(0) = (alpha * (gi) / m_timeStep);
        ground.row(1) = alpha * (gi);

        Eigen::MatrixXd results = factor*current + ground;        
        

        //Velocity, *0.9 for temporariy
        m_velocity[idx] = results.row(0);
        Eigen::Vector3d color = results.row(0).normalized() * 255.0;
        m_vertexColors->SetTuple3(idx, abs(color[0]), abs(color[1]), abs(color[2]));

        //Position        
        m_data->GetPoints()->SetPoint(idx, results.row(1)[0], results.row(1)[1], results.row(1)[2]);
        m_gData->GetPoints()->SetPoint(idx, results.row(1)[0], results.row(1)[1], results.row(1)[2]);
    }


    // //Update Vis Info    
    m_data->GetPoints()->Modified();
    m_gData->GetPoints()->Modified();
    m_actor->GetMapper()->Update();
}

void Beam::UpdateForce(){

    Eigen::Vector3d force(0.0, m_gravity*m_mass, 0.0);

    for(int idxPoint = 0 ; idxPoint< m_data->GetNumberOfPoints() ; idxPoint++){        
            m_force[idxPoint] = force;
    }
    
}

void Beam::SetPointPosition(int idx, double x, double y, double z){    
    m_selectedIdx = idx;
    if(idx == -1) return;

    //Temp
    // idx = m_data->GetNumberOfPoints()-1;
    // m_selectedIdx = idx;


    m_velocity[idx] = Eigen::Vector3d(0.0, 0.0, 0.0);
    m_data->GetPoints()->SetPoint(idx, x, y, z);
    m_data->GetPoints()->Modified();

    m_gData->GetPoints()->SetPoint(idx, x, y, z);
    m_gData->GetPoints()->Modified();
}