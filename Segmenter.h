#ifndef Segmenter_H
#define Segmenter_H

#include <QMainWindow>
#include <vtkSmartPointer.h>
#include <vtkImageViewer2.h>
#include <vtkDICOMImageReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include <vtkActor.h>
#include <vtkCommand.h>
#include <vtkSeedWidget.h>
#include <vtkSeedRepresentation.h>
#include <vtkPointHandleRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkProperty.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <QVTKOpenGLWidget.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextMapper.h>
#include <vtkPolyDataMapper.h>
// some standard vtk headers
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkImageViewer2.h>

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include <itkGDCMImageIO.h>
#include "itkGDCMSeriesFileNames.h"
#include "itkImageToVTKImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"

VTK_MODULE_INIT(vtkRenderingOpenGL2)

static int seedCoords[3][3];
static int z_value;
static int z_slice;

typedef float InputPixelType;//Pixel Type
const unsigned int InputDimension = 3;//Dimension of image
using OutputPixelType = float;

typedef itk::Image< InputPixelType, InputDimension > InputImageType;
typedef itk::Image< OutputPixelType, InputDimension > OutputImageType;

/*************** ITK filters **********************/
typedef itk::ImageSeriesReader< InputImageType > ReaderType;
typedef itk::GDCMImageIO                        ImageIOType;
typedef itk::GDCMSeriesFileNames                NamesGeneratorType;

using CastingFilterType = itk::CastImageFilter< InputImageType, OutputImageType >;
typedef itk::BinaryThresholdImageFilter <InputImageType, InputImageType> BinaryThresholdImageFilterType;
using SmoothingFilterType = itk::CurvatureAnisotropicDiffusionImageFilter<InputImageType, InputImageType>;
using GradientFilterType = itk::GradientMagnitudeRecursiveGaussianImageFilter<InputImageType, InputImageType >;
using SigmoidFilterType = itk::SigmoidImageFilter<InputImageType, InputImageType >;
using FastMarchingFilterType = itk::FastMarchingImageFilter< InputImageType, InputImageType >;
using NodeContainer = FastMarchingFilterType::NodeContainer;
using NodeType = FastMarchingFilterType::NodeType;

typedef itk::CurvatureFlowImageFilter< InputImageType, InputImageType >
CurvatureFlowImageFilterType;
typedef itk::BinaryBallStructuringElement< OutputImageType::PixelType, 3 > StructuringElementType;
using BinaryMorphologicalOpeningImageFilterType = itk::BinaryMorphologicalOpeningImageFilter<OutputImageType, OutputImageType, StructuringElementType>;

namespace Ui {
    class Segmenter;
}

class Segmenter : public QMainWindow {
    Q_OBJECT

public:
    explicit Segmenter(QWidget *parent = 0);
    ~Segmenter();

private slots:
    void openDICOMFolder();
    void drawDICOMSeries(std::string folderDICOM);
    void on_buttonOpenFolder_clicked();
    void on_sliderSlices_sliderMoved(int posicion);
    void on_SeedPoint_clicked();
    //int on_buttonLiverSegmentation_clicked();
    int on_LiverSegmentation_clicked();
    int on_TumorSegmentation_clicked();
   // int on_buttonVesselSegmentation_clicked();
    void on_ShowModel_clicked();

private:
    Ui::Segmenter *ui;
    vtkSmartPointer<vtkDICOMImageReader> reader;
    vtkSmartPointer<vtkImageViewer2> imageViewer;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindow> renderWindow2;

    ReaderType::Pointer readeritk;
    ImageIOType::Pointer gdcmIO;
    NamesGeneratorType::Pointer namesGenerator;

    CastingFilterType::Pointer caster;
    CastingFilterType::Pointer caster1;
    CastingFilterType::Pointer caster2;
    CastingFilterType::Pointer caster3;

    SmoothingFilterType::Pointer smoothing1;
    GradientFilterType::Pointer  gradientMagnitude1;
    SigmoidFilterType::Pointer sigmoid1;
    FastMarchingFilterType::Pointer  fastMarching1;
    NodeContainer::Pointer seeds1;
    BinaryThresholdImageFilterType::Pointer liverThresholdFilter;

    SmoothingFilterType::Pointer smoothing2;
    GradientFilterType::Pointer gradientMagnitude2;
    SigmoidFilterType::Pointer sigmoid2;
    FastMarchingFilterType::Pointer fastMarching2;
    NodeContainer::Pointer seeds2;
    BinaryThresholdImageFilterType::Pointer tumourThresholdFilter;
   
    CurvatureFlowImageFilterType::Pointer smoothing;
    BinaryThresholdImageFilterType::Pointer vesselThresholdFilter;
    BinaryMorphologicalOpeningImageFilterType::Pointer openingFilter;
 
    vtkSmartPointer<vtkImageData> volume1;
    vtkSmartPointer<vtkImageData> volume2;
    vtkSmartPointer<vtkImageData> volume3;

    vtkSmartPointer<vtkRenderWindow> vtkrenderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> interactor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();

    vtkSmartPointer<vtkRenderer> leftRenderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> centerRenderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderer> rightRenderer =
        vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkPolyDataMapper> rightmapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> rightactor =
        vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
        vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    vtkSmartPointer<vtkPolyDataMapper> leftmapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> leftactor =
        vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkPolyDataMapper> centermapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> centeractor =
        vtkSmartPointer<vtkActor>::New();
    vtkSmartPointer<vtkMarchingCubes> surface1 =
        vtkSmartPointer<vtkMarchingCubes>::New();
    vtkSmartPointer<vtkMarchingCubes> surface2 =
        vtkSmartPointer<vtkMarchingCubes>::New();
    vtkSmartPointer<vtkMarchingCubes> surface3 =
        vtkSmartPointer<vtkMarchingCubes>::New();
   
    //vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindowOpenGL2;

    int minSlice;
    int maxSlice;

    void VesselSegmentation(void);

};

#endif // Segmenter_H
