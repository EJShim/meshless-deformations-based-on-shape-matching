#include "Mainwindow.h"
#include <QTimer>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <QTimer>

Mainwindow::Mainwindow()
{	

	// Visualize
	m_renderer = vtkSmartPointer<vtkRenderer>::New();
	m_renderer->AddActor(m_currentObject.GetActor());
	m_renderer->SetBackground(.1, .1, .1);
	m_renderer->ResetCamera();	
	m_renderer->ResetCameraClippingRange();


	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(m_renderer);
	renderWindow->SetSize(1000, 1000);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	renderWindowInteractor->SetInteractorStyle(style);
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	
	

	QTimer* timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(Tick()));
	timer->start(10);
	renderWindowInteractor->Start();
}


Mainwindow::~Mainwindow()
{
}


void Mainwindow::Tick(){

	m_currentObject.ComputeFEM();
	m_currentObject.Update();
	m_renderer->GetRenderWindow()->Render();	

}