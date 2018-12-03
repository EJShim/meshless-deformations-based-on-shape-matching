#pragma once
#include <Beam.h>

#include <iostream>
#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkActor.h>
#include <vtkRenderer.h>


class Mainwindow : public QMainWindow
{
	Q_OBJECT

private:
	vtkSmartPointer<vtkPoints> beamPoints;
	vtkSmartPointer<vtkActor> beamActor;
	vtkSmartPointer<vtkRenderer> m_renderer;
	Beam m_currentObject;

public:
	Mainwindow();
	virtual ~Mainwindow();

public Q_SLOTS:
	void Tick();


};