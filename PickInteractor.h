#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSmartPointer.h>
#include <vtkPointPicker.h>


class PickInteractor : public vtkInteractorStyleTrackballCamera{
    public:
        static PickInteractor* New(){return new PickInteractor;}
        vtkTypeMacro(PickInteractor, vtkInteractorStyleTrackballCamera);

        PickInteractor();
        ~PickInteractor();

        virtual void OnLeftButtonDown();
        virtual void OnMouseMove();
        virtual void OnLeftButtonUp();

    protected:
    ///Picker
    vtkSmartPointer<vtkPointPicker> m_pointPicker;

    ///Picked ID
    int m_id;

    ///Picked Position
   double m_position[3];

   public:
   int GetID(){return m_id;}
   double* GetPosition(){return m_position;}
};