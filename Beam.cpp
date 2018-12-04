#include "Beam.h"
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkPolyDataMapper.h>
#include <vtkMath.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <time.h>
#include <vtkPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractEdges.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>

Beam::Beam(){
    m_actor = NULL;
    Initialize();
}

Beam::~Beam(){

}

void Beam::Initialize(){
    // Create a cube.
    vtkSmartPointer<vtkRectilinearGridToTetrahedra> formMesh = vtkSmartPointer<vtkRectilinearGridToTetrahedra>::New();
    formMesh->SetInput(2, 2, 10, 1, 1, 1, 0.1);    
    formMesh->SetTetraPerCellTo6();
    formMesh->Update();

    m_data = formMesh->GetOutput();
    m_giData = vtkSmartPointer<vtkUnstructuredGrid>::New();
    

    vtkSmartPointer<vtkExtractEdges> edgeExtractor = vtkSmartPointer<vtkExtractEdges>::New();
    edgeExtractor->SetInputData(m_data);
    edgeExtractor->Update();    

    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(edgeExtractor->GetOutputPort());
    tubes->SetRadius(0.1);
    tubes->SetNumberOfSides(6);


    //Show Surface
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceExtractor = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceExtractor->SetInputData(m_data);
    surfaceExtractor->Update();


    // Create a mapper and actor.
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(mapper);
    m_actor->GetProperty()->SetColor(1, 1, 0);

    //Show Temp data
    m_giData->DeepCopy(m_data);
    

    
    vtkSmartPointer<vtkDataSetMapper> giMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    giMapper->SetInputData(m_giData);

    m_giActor =  vtkSmartPointer<vtkActor>::New();
    m_giActor->SetMapper(giMapper);
    m_giActor->GetProperty()->SetColor(1, 0, 0);
    m_giActor->GetProperty()->SetPointSize(50);
    m_giActor->GetProperty()->SetLineWidth(20.0);
    m_giActor->GetProperty()->SetRepresentationToPoints();


    InitializeSystem();
}

void Beam::InitializeSystem(){    

    vtkSmartPointer<vtkPoints> pointSet = m_data->GetPoints();
    Eigen::Matrix3d orgMatrix_V;

    for(int idxTetra = 0 ; idxTetra < m_data->GetNumberOfCells() ; idxTetra++){
        vtkSmartPointer<vtkCell> tetra = m_data->GetCell(idxTetra);

        //Four Points of Current Tetrahedra
        Eigen::Vector3d pt0 ( pointSet->GetPoint(tetra->GetPointId(0)) );
        Eigen::Vector3d pt1 ( pointSet->GetPoint(tetra->GetPointId(1)) );
        Eigen::Vector3d pt2 ( pointSet->GetPoint(tetra->GetPointId(2)) );
        Eigen::Vector3d pt3 ( pointSet->GetPoint(tetra->GetPointId(3)) );

        //Original Matrix
        Eigen::Matrix3d orgMatrix_V;
        orgMatrix_V.col(0) = pt1 - pt0;
        orgMatrix_V.col(1) = pt2 - pt0;
        orgMatrix_V.col(2) = pt3 - pt0;
        
        //What is this?
        Eigen::PartialPivLU<Eigen::Matrix3d> lu(orgMatrix_V);        
        m_orgInverseMatrix.push_back(lu.inverse());
    }

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
    m_giData->GetPointData()->SetScalars(m_vertexColors);

    m_iCenterOfMass /= nPoints;

    for(int idx = 0 ; idx<nPoints ; idx++){
        m_qi.push_back( Eigen::Vector3d(m_data->GetPoint(idx)) - m_iCenterOfMass  );
    }

    //Initialize Factors
    m_Identity = Eigen::Matrix3d::Identity();
    m_E = Eigen::MatrixXd::Zero(6, 6);
    m_E << 1.0 - m_poissonsRatio, m_poissonsRatio, m_poissonsRatio, 0.0, 0.0, 0.0,
		m_poissonsRatio, 1.0 - m_poissonsRatio, m_poissonsRatio, 0.0, 0.0, 0.0,
		m_poissonsRatio, m_poissonsRatio, 1.0 - m_poissonsRatio, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0 - 2.0 * m_poissonsRatio, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 1.0 - 2.0 * m_poissonsRatio, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 1.0 - 2.0 * m_poissonsRatio;
    m_E = (m_youngsModulus / ((1.0f + m_poissonsRatio)*(1.0 - 2.0 * m_poissonsRatio)))*m_E;
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
    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d pi = Eigen::Vector3d(m_data->GetPoint(idx)) - CenterOfMass;
        Apq += m_mass * pi * m_qi[idx].transpose();
    }


    //Inverse SQRT???
    Eigen::Matrix3d S = Apq.transpose() * Apq;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S);
    Eigen::Matrix3d sqrt_S = es.operatorInverseSqrt();

    Eigen::Matrix3d R = Apq * sqrt_S;

    
    double alpha = 0.3;


    
    Eigen::Matrix2d factor(2, 2);
    factor  << 1, -alpha/m_timeStep,
                m_timeStep, 1-alpha;    
    Eigen::MatrixXd current(2, 3);
    Eigen::MatrixXd ground(2, 3);



    //Update Position
    for(int idx = 9 ; idx < nPoints ; idx++){

        Eigen::Vector3d gi =  R*m_qi[idx]+CenterOfMass;
        Eigen::Vector3d xi = Eigen::Vector3d(m_data->GetPoint(idx));

        current.row(0) = m_velocity[idx];
        current.row(1) = xi;

        // std::cout << gi.transpose() << "," << current.row(1) << std::endl;

        ground.row(0) = (alpha * (gi) / m_timeStep) + m_timeStep*(m_force[idx] / m_mass);
        ground.row(1) = alpha * (gi);

        Eigen::MatrixXd results = factor*current + ground;        
        

        //Velocity
        m_velocity[idx] = results.row(0);
        Eigen::Vector3d color = results.row(0).normalized() * 255.0;
        m_vertexColors->SetTuple3(idx, abs(color[0]), abs(color[1]), abs(color[2]));

        //Position        
        m_data->GetPoints()->SetPoint(idx, results.row(1)[0], results.row(1)[1], results.row(1)[2]);
        m_giData->GetPoints()->SetPoint(idx, gi[0], gi[1], gi[2]);
    }


    // //Update Vis Info    
    m_data->GetPoints()->Modified();
    m_giData->GetPoints()->Modified();    
}


void Beam::ComputeFEM(){
    UpdateForce();

    

    vtkSmartPointer<vtkPoints> pointSet = m_data->GetPoints();
    

    for(int idxTetra = 0 ; idxTetra < m_data->GetNumberOfCells() ; idxTetra++){
        vtkSmartPointer<vtkCell> tetra = m_data->GetCell(idxTetra);
        
        int id0 = tetra->GetPointId(0);
        int id1 = tetra->GetPointId(1);
        int id2 = tetra->GetPointId(2);
        int id3 = tetra->GetPointId(3);
        //Four Points of Current Tetrahedra
        Eigen::Vector3d pt0 ( pointSet->GetPoint(id0) );
        Eigen::Vector3d pt1 ( pointSet->GetPoint(id1) );
        Eigen::Vector3d pt2 ( pointSet->GetPoint(id2) );
        Eigen::Vector3d pt3 ( pointSet->GetPoint(id3) );

        //Deform Matrix
        Eigen::Matrix3d deformMatrix;
        deformMatrix.col(0) = pt1 - pt0;
        deformMatrix.col(1) = pt2 - pt0;
        deformMatrix.col(2) = pt3 - pt0;


        Eigen::Matrix3d PMatrix = deformMatrix * m_orgInverseMatrix[idxTetra];
        Eigen::Matrix3d strainMatrix = 0.5 * (PMatrix + PMatrix.transpose()) - m_Identity;
        Eigen::VectorXd strainVector(6);
        strainVector << strainMatrix(0, 0), strainMatrix(1, 1), strainMatrix(2, 2), strainMatrix(0, 1), strainMatrix(1, 2), strainMatrix(2, 0);
        Eigen::VectorXd stressVector(6);
        stressVector = m_E * strainVector;
        Eigen::Matrix3d stressMatrix;
        stressMatrix << stressVector(0, 0), stressVector(3, 0), stressVector(5, 0),
                        stressVector(3, 0), stressVector(1, 0), stressVector(4, 0),
                        stressVector(5, 0), stressVector(4, 0), stressVector(2, 0);


        //Get Direction Vector of Current Tet
        Eigen::Vector3d p10(pt1 - pt0);
        Eigen::Vector3d p20(pt2 - pt0);
        Eigen::Vector3d p30(pt3 - pt0);
        Eigen::Vector3d p21(pt2 - pt1);
        Eigen::Vector3d p31(pt3 - pt1);

        //Calculate Force
        Eigen::Vector3d force012 = stressMatrix * (p10.cross(p20))*(-1.0/3.0);
        Eigen::Vector3d force013 = stressMatrix * (p30.cross(p10))*(-1.0/3.0);
        Eigen::Vector3d force023 = stressMatrix * (p20.cross(p30))*(-1.0/3.0);
        Eigen::Vector3d force123 = stressMatrix * (p31.cross(p21))*(-1.0/3.0);

        //Update Force
        m_force[id0] += force012 + force013 + force023;
        m_force[id1] += force012 + force013 + force123;
        m_force[id2] += force012 + force023 + force123;
        m_force[id3] += force013 + force023 + force123;

                        
    }

}

void Beam::UpdateForce(){

    Eigen::Vector3d force(0.0, m_gravity*m_mass, 0.0);

    for(int idxPoint = 0 ; idxPoint< m_data->GetNumberOfPoints() ; idxPoint++){        
            m_force[idxPoint] = force;
    }
    
}

void Beam::Update(){
    vtkMath::RandomSeed(time(NULL));

    for(int idx = 9 ; idx < m_data->GetNumberOfPoints() ; idx++){
        Eigen::Vector3d acceleration = m_force[idx] / m_mass;
        m_velocity[idx] += acceleration * m_timeStep;
        //Change Position
        Eigen::Vector3d position ( m_data->GetPoint(idx) );
        position += m_velocity[idx];

        m_data->GetPoints()->SetPoint(idx, position[0], position[1], position[2]);
    }

    m_data->GetPoints()->Modified();
}

void Beam::SetPointPosition(int idx, double x, double y, double z){    
    m_selectedIdx = idx;
    if(idx == -1) return;

    //Temp
    idx = m_data->GetNumberOfPoints()-1;
    m_selectedIdx = idx;


    m_velocity[idx] = Eigen::Vector3d(0.0, 0.0, 0.0);
    m_data->GetPoints()->SetPoint(idx, x, y, z);
    m_data->GetPoints()->Modified();
}