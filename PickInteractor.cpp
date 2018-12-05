#include <PickInteractor.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAbstractPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
PickInteractor::PickInteractor(){
    m_pointPicker = vtkSmartPointer<vtkPointPicker>::New();
    m_id = -1;
    m_position[0] = 0.0;
    m_position[1] = 0.0;
    m_position[2] = 0.0;

    m_startPosition[0] = 0.0;
    m_startPosition[1] = 0.0;
    m_startPosition[2] = 0.0;
}

PickInteractor::~PickInteractor(){

}

void PickInteractor::OnLeftButtonDown(){    

    m_pointPicker->Pick(this->Interactor->GetEventPosition()[0], this->Interactor->GetEventPosition()[1], 0, this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

    m_pointPicker->GetPickPosition(m_position);
    m_pointPicker->GetPickPosition(m_startPosition);
    m_id = m_pointPicker->GetPointId();
    

    if(m_id == -1){
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }
}

void PickInteractor::OnMouseMove(){
    if(m_id == -1) {
        vtkInteractorStyleTrackballCamera::OnMouseMove();
        return;
    }

    m_pointPicker->Pick(this->Interactor->GetEventPosition()[0], this->Interactor->GetEventPosition()[1], 0,
this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

      m_pointPicker->GetPickPosition(m_position);

    // std::cout << "ID: " << m_id << ", Picked value: " << m_position[0] << " " << m_position[1] << " " << m_position[2] << std::endl;

}

void PickInteractor::OnLeftButtonUp(){
    if(m_id == -1) {
        vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
        return;
    }

    m_id = -1;
    m_position[0] = 0.0;
    m_position[1] = 0.0;
    m_position[2] = 0.0;

    m_startPosition[0] = 0.0;
    m_startPosition[1] = 0.0;
    m_startPosition[2] = 0.0;
}