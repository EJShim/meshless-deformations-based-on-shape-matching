#include "deformableMesh.h"
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkMath.h>
#include <time.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkOBBDicer.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>

deformableMesh::deformableMesh(vtkSmartPointer<vtkPolyData> data){
    m_actor = NULL;
    Initialize(data);
}

deformableMesh::~deformableMesh(){

}

void deformableMesh::Initialize(vtkSmartPointer<vtkPolyData> data){
  
    //Clearn..
    vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanPolyData->SetInputData(data);
    cleanPolyData->Update();

    m_data = cleanPolyData->GetOutput();

    // Create a mapper and actor.
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(m_data);
    mapper->Update();

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(mapper);
    m_actor->GetProperty()->SetColor(1, 0, 0);
    m_actor->GetProperty()->SetRepresentationToWireframe();
    m_actor->GetProperty()->RenderPointsAsSpheresOn();
    m_actor->GetProperty()->SetPointSize(15);


    //Create Debug Actor
    m_gData = vtkSmartPointer<vtkPolyData>::New();
    m_gData->DeepCopy(m_data);

    vtkSmartPointer<vtkPolyDataMapper> groundMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    groundMapper->SetInputData(m_gData);
    groundMapper->Update();

    m_gActor = vtkSmartPointer<vtkActor>::New();
    m_gActor->SetMapper(groundMapper);
    m_gActor->GetProperty()->SetColor(.9, .9, 0);
    // m_gActor->GetProperty()->SetRepresentationToPoints();
    m_gActor->GetProperty()->RenderPointsAsSpheresOn();
    m_gActor->GetProperty()->SetPointSize(15);
    
    

    InitializeSystem();
    MakeCluster();
}

void deformableMesh::InitializeSystem(){    

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
        m_vertexColors->InsertNextTuple3(255, 254, 0);
        
    }
    m_gData->GetPointData()->SetScalars(m_vertexColors);

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


void deformableMesh::MakeCluster(){
//Initialize Clustering
    vtkSmartPointer<vtkOBBDicer> dicer = vtkSmartPointer<vtkOBBDicer>::New();
    dicer->SetInputData(m_data);
    dicer->SetNumberOfPieces(4);
    dicer->SetDiceModeToSpecifiedNumberOfPieces();
    dicer->Update();

    std::cout << "Actual Number of Cluster : " << dicer->GetNumberOfActualPieces() << std::endl;

    vtkSmartPointer<vtkThreshold> selector = vtkSmartPointer<vtkThreshold>::New();
    selector->SetInputArrayToProcess(0, 0, 0, 0, "vtkOBBDicer_GroupIds");
    selector->SetInputConnection(dicer->GetOutputPort());
    selector->AllScalarsOff();


    vtkMath::RandomSeed(8775070); // for reproducibility
    for(int cluster=0 ; cluster<dicer->GetNumberOfActualPieces() ; cluster++){
        selector->ThresholdBetween(cluster, cluster);

        vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
        geometryFilter->SetInputConnection(selector->GetOutputPort());
        geometryFilter->Update();

        vtkSmartPointer<vtkPolyData> subCluster = geometryFilter->GetOutput();
        std::vector<int> idArray;



        //Visualize Different Cluster Color
        float color[3];
        color[0] = vtkMath::Random(0, 255);
        color[1] = vtkMath::Random(0, 255);
        color[2] = vtkMath::Random(0, 255);
        
        for(int idx = 0 ; idx < subCluster->GetNumberOfPoints() ; idx++){
            double* p = subCluster->GetPoint(idx);
            int pointId = m_data->FindPoint(p);
            
            idArray.push_back(pointId);
            m_vertexColors->SetTuple3(pointId, color[0], color[1], color[2]);
        }
        m_clusterID.push_back(idArray);
    }    
}

void deformableMesh::ComputeMesheless(){
    //Get CenterOfMass    
    int nPoints = m_data->GetNumberOfPoints();    


    m_cCenterOfMass = Eigen::Vector3d(0, 0, 0);        
    for(int idx = 0 ; idx < nPoints ; idx++){
        m_cCenterOfMass += Eigen::Vector3d(m_data->GetPoint(idx));
    }
    m_cCenterOfMass /= nPoints;



    //Calculate Apq
    Eigen::Matrix3d Apq = Eigen::MatrixXd::Zero(3,3);
    Eigen::MatrixXd APQ = Eigen::MatrixXd::Zero(3,9);
    
    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d pi = Eigen::Vector3d(m_data->GetPoint(idx)) - m_cCenterOfMass;
        
        Apq += m_mass * pi * m_qi[idx].transpose();
        APQ += m_mass * pi * m_Qi[idx].transpose();
    }

    


    //Inverse SQRT???
    Eigen::Matrix3d S = Apq.transpose() * Apq;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S);
    Eigen::Matrix3d sqrt_S = es.operatorInverseSqrt();

    Eigen::Matrix3d R = Apq * sqrt_S;
    Eigen::MatrixXd RR = Eigen::MatrixXd::Zero(3, 9);
    RR.block<3,3>(0,0) = R;
    


    //Matrix A
    Eigen::Matrix3d A = Apq * m_Aqq;
    Eigen::MatrixXd AA = APQ * m_AQQ;


    double alpha = 0.9;
    double beta = 0.9;
    double damping = 0.05;


    
    Eigen::Matrix2d factor(2, 2);
    factor  << 1.0-damping, -alpha/m_timeStep,
                m_timeStep, 1-alpha;    
    Eigen::MatrixXd current(2, 3);
    Eigen::MatrixXd ground(2, 3);


    // std::cout << AA << std::endl << std::endl;

    //Update Position
    for(int idx = 0 ; idx < nPoints ; idx++){

        Eigen::Vector3d pi = Eigen::Vector3d(m_data->GetPoint(idx)) - m_cCenterOfMass;
        

        Eigen::Vector3d gi =  ( beta*AA + (1-beta)*RR )*m_Qi[idx]+ m_cCenterOfMass;
        Eigen::Vector3d xi = Eigen::Vector3d(m_data->GetPoint(idx));

        current.row(0) = m_velocity[idx];
        current.row(1) = xi;        


        ground.row(0) = (alpha * (gi) / m_timeStep) + m_timeStep*m_force[idx]/m_mass;
        ground.row(1) = alpha * (gi);

        Eigen::MatrixXd results = factor*current + ground;        
        

        //Velocity, *0.9 for temporariy
        m_velocity[idx] = results.row(0);
        Eigen::Vector3d color = results.row(0).normalized() * 255.0;
        // m_vertexColors->SetTuple3(idx, abs(color[0]), abs(color[1]), abs(color[2]));

        //Position        
        m_data->GetPoints()->SetPoint(idx, results.row(1)[0], results.row(1)[1], results.row(1)[2]);
        m_gData->GetPoints()->SetPoint(idx, gi[0], gi[1], gi[2]);
    }


    // //Update Vis Info    
    m_data->GetPoints()->Modified();
    m_gData->GetPoints()->Modified();
    m_actor->GetMapper()->Update();


    UpdateForce();
}

void deformableMesh::UpdateForce(){

    Eigen::Vector3d force(0.0, m_gravity*m_mass, 0.0);

    for(int idxPoint = 0 ; idxPoint< m_data->GetNumberOfPoints() ; idxPoint++){        
            m_force[idxPoint] = force;
    }
    
}

void deformableMesh::ApplyForce(int idx, double x, double y, double z){
    m_selectedIdx = idx;
    if(idx == -1) return;

    m_force[idx] = GetTotalMass()*m_timeStep * Eigen::Vector3d(x, y, z);

}


double* deformableMesh::GetCurrentSelectedPosition(int idx){
    return m_gData->GetPoint(idx);
}

void deformableMesh::SetPointPosition(int idx, double x, double y, double z){    
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

double deformableMesh::GetTotalMass(){
    return m_mass * m_data->GetNumberOfPoints();
}