#pragma once
#include <Beam.h>

#include <iostream>
#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <QVTKOpenGLWidget.h>
#include <PickInteractor.h>

class Mainwindow : public QMainWindow
{
	Q_OBJECT

private:

    QVTKOpenGLWidget* m_rendererWidget;
	vtkSmartPointer<vtkRenderer> m_renderer;
	vtkSmartPointer<PickInteractor> m_interactorStyle;
	Beam m_currentObject;

protected:
    QWidget* InitCentralWidget();

public:
	Mainwindow();
	virtual ~Mainwindow();

public Q_SLOTS:
	void Tick();


};