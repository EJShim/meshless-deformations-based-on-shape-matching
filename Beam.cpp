#include "Beam.h"
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkMath.h>
#include <time.h>



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

    m_edgeExtractor = vtkSmartPointer<vtkExtractEdges>::New();
    m_edgeExtractor->SetInputData(m_data);
    m_edgeExtractor->Update();

    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputConnection(m_edgeExtractor->GetOutputPort());
    tubes->SetRadius(0.1);
    tubes->SetNumberOfSides(6);


    // Create a mapper and actor.
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());

    m_actor = vtkSmartPointer<vtkActor>::New();
    m_actor->SetMapper(mapper);
    m_actor->GetProperty()->SetColor(1, 1, 0);


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

    //Initialize Force and Velocity

    m_iCenterOfMass = Eigen::Vector3d(0, 0, 0);
    int nPoints = m_data->GetNumberOfPoints();
    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d velocity(0, 0, 0);
        Eigen::Vector3d force(0, 0, 0);

        m_velocity.push_back(velocity);
        m_force.push_back(force);

        m_iCenterOfMass += Eigen::Vector3d(m_data->GetPoint(idx));
    }

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

    std::cout << sqrt_S << std::endl<<std::endl;

    // std::cout << R << std::endl<< std::endl;

    // std::vector<Eigen::Vector3d> gi;


    double alpha = 0.3;


    Eigen::MatrixXd results(2, 3);    
    Eigen::Matrix2d factor(2, 2);
    factor  << 1, -alpha/m_timeStep,
                m_timeStep, 1-alpha;
    Eigen::MatrixXd current(2, 3);
    Eigen::MatrixXd ground(2, 3);


    for(int idx = 0 ; idx < nPoints ; idx++){
        Eigen::Vector3d gi =  R*m_qi[idx]+m_iCenterOfMass ;

        
        current.col(0) = m_velocity[idx];
        current.col(1) = Eigen::Vector3d(m_data->GetPoint(idx));

        ground.col(0) = alpha * gi / m_timeStep;
        ground.col(1) = alpha * gi;


        results = factor*current + ground;

        //Velocity
        m_velocity[idx] = results.col(0);

        //Position
        Eigen::Vector3d position = results.col(1);

        
        m_data->GetPoints()->SetPoint(idx, position[0], position[1], position[2]);
    }

    m_edgeExtractor->Modified();


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

    m_edgeExtractor->Modified();
}