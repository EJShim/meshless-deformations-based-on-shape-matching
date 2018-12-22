#include "Mainwindow.h"
#include <QTimer>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <QCoreApplication>
#include <QDir>

#include <QTimer>

Mainwindow::Mainwindow()
{	
	
    InitObjects();
	
    //Initialize Renderer
    setCentralWidget(InitCentralWidget());

    //Set Event Loop
	QTimer* timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(Tick()));
	timer->start(10);
}


Mainwindow::~Mainwindow()
{
}

void Mainwindow::InitObjects(){

    std::string stlPath = QCoreApplication::applicationDirPath().toStdString() + "/bunny.stl";
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(stlPath.c_str());
    reader->Update();

    // TestObject 1. deformableMesh
    vtkSmartPointer<vtkRectilinearGridToTetrahedra> formMesh = vtkSmartPointer<vtkRectilinearGridToTetrahedra>::New();
    formMesh->SetInput(20, 20, 100, 10.0, 10.0, 10.0, 0.1);        
    formMesh->Update();
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceExtractor = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceExtractor->SetInputData(formMesh->GetOutput());
    surfaceExtractor->Update();

    // TestObject2. Sphere - has backface problem
    vtkSmartPointer<vtkCylinderSource> cylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
    cylinderSource->SetCenter(0.0, 0.0, 0.0);
    cylinderSource->SetRadius(5.0);
    cylinderSource->SetHeight(7.0);
    cylinderSource->SetResolution(10);
    cylinderSource->Update();

    m_currentObject = new deformableMesh(surfaceExtractor->GetOutput());


    //Initialize Arrow Actor
    m_arrowData = vtkSmartPointer<vtkArrowSource>::New();
    m_arrowData->SetShaftRadius(0.01);
    m_arrowData->SetTipRadius(0.02);
    m_arrowData->SetTipLength(0.1);
    m_arrowData->Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(m_arrowData->GetOutputPort());

    m_arrowActor = vtkSmartPointer<vtkActor>::New();
    m_arrowActor->SetMapper(mapper);
    m_arrowActor->GetProperty()->SetColor(0.0, 0.1, 0.8);

}


void Mainwindow::Tick(){

    
    int pickID = m_interactorStyle->GetID();

    
    double* position = m_interactorStyle->GetPosition();
    
    if(pickID == -1){
        m_renderer->RemoveActor(m_arrowActor);
    }else{
        
        m_renderer->AddActor(m_arrowActor);

        double* startPosition = m_currentObject->GetCurrentSelectedPosition(pickID);
        double start[3] = {startPosition[0], startPosition[1], startPosition[2]};
        double end[3] = {position[0], position[1], position[2]};
        UpdateArrow(start, end);

        double force[3];
        vtkMath::Subtract(end, start, force);
        m_currentObject->ApplyForce(pickID, force[0], force[1], force[2] );
    }

    m_currentObject->ComputeMesheless();
    m_renderer->ResetCameraClippingRange();
	m_renderer->GetRenderWindow()->Render();	

}

void Mainwindow::UpdateArrow(double startPoint[3], double endPoint[3]){
    // Compute a basis
    double normalizedX[3];
    double normalizedY[3];
    double normalizedZ[3];

    // The X axis is a vector from start to end
    vtkMath::Subtract(endPoint, startPoint, normalizedX);
    double length = vtkMath::Norm(normalizedX);
    vtkMath::Normalize(normalizedX);

    // The Z axis is an arbitrary vector cross X
    double arbitrary[3];
    arbitrary[0] = vtkMath::Random(-10,10);
    arbitrary[1] = vtkMath::Random(-10,10);
    arbitrary[2] = vtkMath::Random(-10,10);
    vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
    vtkMath::Normalize(normalizedZ);

    // The Y axis is Z cross X
    vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();

    // Create the direction cosine matrix
    matrix->Identity();
    for (unsigned int i = 0; i < 3; i++)
    {
        matrix->SetElement(i, 0, normalizedX[i]);
        matrix->SetElement(i, 1, normalizedY[i]);
        matrix->SetElement(i, 2, normalizedZ[i]);
    }    

    // Apply the transforms
    vtkSmartPointer<vtkTransform> transform =  vtkSmartPointer<vtkTransform>::New();
    transform->Translate(startPoint);
    transform->Concatenate(matrix);
    transform->Scale(length, length, length);

    m_arrowActor->SetUserTransform(transform);
}

QWidget* Mainwindow::InitCentralWidget(){
    m_rendererWidget = new QVTKOpenGLWidget();

    // Visualize
	m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_renderer->SetBackground(.1, .1, .1);
	// m_renderer->AddActor(m_currentObject->GetActor());    
    m_renderer->AddActor(m_currentObject->GetDebugActor());
	m_renderer->ResetCamera();	
	m_renderer->ResetCameraClippingRange();


    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
	renderWindow->AddRenderer(m_renderer);	
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    m_interactorStyle = vtkSmartPointer<PickInteractor>::New();
	renderWindowInteractor->SetInteractorStyle(m_interactorStyle);
	renderWindowInteractor->SetRenderWindow(renderWindow);	

    m_rendererWidget->SetRenderWindow(renderWindow);
    renderWindow->Render();


    #ifdef __APPLE__
        //Force to use GL>3.2,, mac default is 2.1
        QSurfaceFormat::setDefaultFormat(m_rendererWidget->defaultFormat());
        m_rendererWidget->setFormat(m_rendererWidget->defaultFormat());        
    #endif


    return m_rendererWidget;
}