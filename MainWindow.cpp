#include "Mainwindow.h"
#include <QTimer>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <QTimer>

Mainwindow::Mainwindow()
{	
	
    //Initialize Toolbar
	
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


void Mainwindow::Tick(){
    //Calculate Based on FERM
	// m_currentObject.ComputeFEM();
	

    m_currentObject.ComputeMesheless();
    // m_currentObject.Update();



	m_renderer->GetRenderWindow()->Render();	

}

QWidget* Mainwindow::InitCentralWidget(){
    m_rendererWidget = new QVTKOpenGLWidget();

    // Visualize
	m_renderer = vtkSmartPointer<vtkRenderer>::New();
    m_renderer->SetBackground(.1, .1, .1);
	m_renderer->AddActor(m_currentObject.GetActor());
	m_renderer->ResetCamera();	
	m_renderer->ResetCameraClippingRange();


    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
	renderWindow->AddRenderer(m_renderer);	
	// vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	// vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	// renderWindowInteractor->SetInteractorStyle(style);
	// renderWindowInteractor->SetRenderWindow(renderWindow);	

    m_rendererWidget->SetRenderWindow(renderWindow);
    renderWindow->Render();


    #ifdef __APPLE__
        //Force to use GL>3.2,, mac default is 2.1
        QSurfaceFormat::setDefaultFormat(m_rendererWidget->defaultFormat());
        m_rendererWidget->setFormat(m_rendererWidget->defaultFormat());        
    #endif


    return m_rendererWidget;
}