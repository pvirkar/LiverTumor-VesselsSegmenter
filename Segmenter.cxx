#include "Segmenter.h"
#include "ui_Segmenter.h"

#include <QFileDialog>
#include <QDebug>


class vtkSeedCallback : public vtkCommand
{
public:
    static vtkSeedCallback* New()
    {
        return new vtkSeedCallback;
    }

    vtkSeedCallback():m_nCount(0)
    {}

    virtual void Execute(vtkObject* caller, unsigned long event, void* calldata)
    {
        vtkSeedWidget* pSeedWidget = reinterpret_cast<vtkSeedWidget*>(caller);

        if (m_nCount >= 3)
            pSeedWidget->CompleteInteraction();

        if (event == vtkCommand::PlacePointEvent)
        {
            std::cout << "Point placed, total of: "
                << this->SeedRepresentation->GetNumberOfSeeds() << std::endl;
            m_nCount++;
        }
        if (event == vtkCommand::InteractionEvent)
        {
            if (calldata)
            {
                std::cout << "Interacting with seed : "
                    << *(static_cast<int*>(calldata)) << std::endl;
            }
        }


        std::cout << "List of seeds (Display coordinates):" << std::endl;
        for (vtkIdType i = 0; i < this->SeedRepresentation->GetNumberOfSeeds(); i++)
        {
            double pos[3];
            this->SeedRepresentation->GetSeedDisplayPosition(i, pos);
            std::cout << "(" << pos[0] << " " << pos[1] << " " << pos[2] << ")" << std::endl;

            seedCoords[i][0] = pos[0];
            seedCoords[i][1] = pos[1];
           // seedCoords[i][2] = z_slice;
            
            std::cout << "x " << seedCoords[i][0] << "& " << "y " << seedCoords[i][1] << std::endl;
           
            
        }

    }

    void SetRepresentation(vtkSmartPointer<vtkSeedRepresentation> rep) { this->SeedRepresentation = rep; }
private:
    int m_nCount;	// to limit the seed point
    vtkSmartPointer<vtkSeedRepresentation> SeedRepresentation;
};
//vtkStandardNewMacro(vtkSeedCallback)


// helper class to format slice status message
class StatusMessage {
public:
    static std::string Format(int slice, int maxSlice) {
        std::stringstream tmp;
        tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
        z_slice = slice + 1;
        for (int i = 0; i < 3; ++i)
        {
            seedCoords[i][2] = z_slice;
        }

        return tmp.str();
    }
};


// Define own interaction style
class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
    static myVtkInteractorStyleImage* New();
    vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
    vtkImageViewer2* _ImageViewer;
    vtkTextMapper* _StatusMapper;
    int _Slice;
    int _MinSlice;
    int _MaxSlice;

public:
    void SetImageViewer(vtkImageViewer2* imageViewer) {
        _ImageViewer = imageViewer;
        _MinSlice = imageViewer->GetSliceMin();
        _MaxSlice = imageViewer->GetSliceMax();
        _Slice = _MinSlice;
        cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
    }

    void SetStatusMapper(vtkTextMapper* statusMapper) {
        _StatusMapper = statusMapper;
    }


protected:
    void MoveSliceForward() {
        if (_Slice < _MaxSlice) {
            _Slice += 1;
            cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
            _ImageViewer->SetSlice(_Slice);
            std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
            _StatusMapper->SetInput(msg.c_str());
            _ImageViewer->Render();
        }
    }

    void MoveSliceBackward() {
        if (_Slice > _MinSlice) {
            _Slice -= 1;
            cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
            _ImageViewer->SetSlice(_Slice);
            std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
            _StatusMapper->SetInput(msg.c_str());
            _ImageViewer->Render();
        }
    }


    virtual void OnMouseWheelForward() {
        //std::cout << "Scrolled mouse wheel forward." << std::endl;
        MoveSliceForward();
       
    }


    virtual void OnMouseWheelBackward() {
        //std::cout << "Scrolled mouse wheel backward." << std::endl;
        if (_Slice > _MinSlice) {
            MoveSliceBackward();
        }
       
    }
};

vtkStandardNewMacro(myVtkInteractorStyleImage);


Segmenter::Segmenter(QWidget *parent) : QMainWindow(parent), ui(new Ui::Segmenter)
{
    ui->setupUi(this);
   
    QFont fontB("Arial", 12, QFont::Bold);
    ui->label->setFont(fontB);
    ui->label->setText("Click on 'Open DICOM Folder' to load DICOM series.");

    reader = vtkSmartPointer<vtkDICOMImageReader>::New();
    imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow2 = vtkSmartPointer<vtkRenderWindow>::New();

    readeritk = ReaderType::New();
    gdcmIO = ImageIOType::New();
    namesGenerator = NamesGeneratorType::New();

    caster = CastingFilterType::New();
    caster1 = CastingFilterType::New();
    caster2 = CastingFilterType::New();
    caster3 = CastingFilterType::New();
 
    smoothing1 = SmoothingFilterType::New(); 
    gradientMagnitude1 = GradientFilterType::New();   
    sigmoid1 = SigmoidFilterType::New();
    fastMarching1 = FastMarchingFilterType::New();
    seeds1 = NodeContainer::New();
    liverThresholdFilter = BinaryThresholdImageFilterType::New();

    smoothing2 = SmoothingFilterType::New();
    gradientMagnitude2 = GradientFilterType::New();
    sigmoid2 = SigmoidFilterType::New();
    fastMarching2 = FastMarchingFilterType::New();
    seeds2 = NodeContainer::New();
    tumourThresholdFilter = BinaryThresholdImageFilterType::New();
   
    smoothing = CurvatureFlowImageFilterType::New();
    vesselThresholdFilter = BinaryThresholdImageFilterType::New();
    openingFilter = BinaryMorphologicalOpeningImageFilterType::New();   

    volume1 = vtkSmartPointer<vtkImageData>::New();
    volume2 = vtkSmartPointer<vtkImageData>::New();
    volume3 = vtkSmartPointer<vtkImageData>::New();

    leftRenderer = vtkSmartPointer<vtkRenderer>::New();
    rightmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    rightactor = vtkSmartPointer<vtkActor>::New();
    style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    leftmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    leftactor = vtkSmartPointer<vtkActor>::New();
    centermapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    centeractor = vtkSmartPointer<vtkActor>::New();
    surface1 = vtkSmartPointer<vtkMarchingCubes>::New();
    surface2 = vtkSmartPointer<vtkMarchingCubes>::New();
    surface3 = vtkSmartPointer<vtkMarchingCubes>::New();

   // renderWindowOpenGL2 = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkrenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

}

Segmenter::~Segmenter() {
    delete ui;
}

void Segmenter::openDICOMFolder() {
    QString folderNameDICOM = QFileDialog::getExistingDirectory(this, tr("Open DICOM Folder"), QDir::currentPath(), QFileDialog::ShowDirsOnly);
    drawDICOMSeries(folderNameDICOM.toUtf8().constData());
}

void Segmenter::drawDICOMSeries(std::string folderDICOM) {

    reader->SetDirectoryName(folderDICOM.c_str());
    reader->Update();

    namesGenerator->SetInputDirectory(folderDICOM.c_str());
    const itk::ImageSeriesReader< InputImageType >::FileNamesContainer& filenames =
        namesGenerator->GetInputFileNames();
    readeritk->SetFileNames(filenames);

    //Exceptional handling
    try
    {
        readeritk->Update();
    }
    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
    }

    imageViewer->SetInputConnection(reader->GetOutputPort());

    // slice status message
    vtkSmartPointer<vtkTextProperty> sliceTextProp = vtkSmartPointer<vtkTextProperty>::New();
    sliceTextProp->SetFontFamilyToCourier();
    sliceTextProp->SetFontSize(20);
    sliceTextProp->SetVerticalJustificationToBottom();
    sliceTextProp->SetJustificationToLeft();

    vtkSmartPointer<vtkTextMapper> sliceTextMapper = vtkSmartPointer<vtkTextMapper>::New();
    std::string msg = StatusMessage::Format(imageViewer->GetSliceMin(), imageViewer->GetSliceMax());
    sliceTextMapper->SetInput(msg.c_str());
    sliceTextMapper->SetTextProperty(sliceTextProp);

    vtkSmartPointer<vtkActor2D> sliceTextActor = vtkSmartPointer<vtkActor2D>::New();
    sliceTextActor->SetMapper(sliceTextMapper);
    sliceTextActor->SetPosition(15, 10);

    // usage hint message
    vtkSmartPointer<vtkTextProperty> usageTextProp = vtkSmartPointer<vtkTextProperty>::New();
    usageTextProp->SetFontFamilyToCourier();
    usageTextProp->SetFontSize(14);
    usageTextProp->SetVerticalJustificationToTop();
    usageTextProp->SetJustificationToLeft();

    vtkSmartPointer<vtkTextMapper> usageTextMapper = vtkSmartPointer<vtkTextMapper>::New();
    usageTextMapper->SetInput("- Scroll through slices with mouse wheel");
    usageTextMapper->SetTextProperty(usageTextProp);

    vtkSmartPointer<vtkActor2D> usageTextActor = vtkSmartPointer<vtkActor2D>::New();
    usageTextActor->SetMapper(usageTextMapper);
    usageTextActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    usageTextActor->GetPositionCoordinate()->SetValue(0.05, 0.95);

    //ui->widget->SetRenderWindow(viewer->GetRenderWindow());
    //viewer->SetupInteractor(ui->widget->GetInteractor());

    vtkSmartPointer<myVtkInteractorStyleImage> myInteractorStyle =
        vtkSmartPointer<myVtkInteractorStyleImage>::New();

    // make imageviewer2 and sliceTextMapper visible to our interactorstyle
    // to enable slice status message updates when scrolling through the slices
    myInteractorStyle->SetImageViewer(imageViewer);
    myInteractorStyle->SetStatusMapper(sliceTextMapper);

    ui->widget->SetRenderWindow(renderWindow.Get());
    renderWindowInteractor = ui->widget->GetInteractor();
    imageViewer->SetRenderWindow(ui->widget->GetRenderWindow());
    imageViewer->SetupInteractor(renderWindowInteractor);

    renderWindowInteractor->SetInteractorStyle(myInteractorStyle);
    // add slice status message and usage hint message to the renderer
    imageViewer->GetRenderer()->AddActor2D(sliceTextActor);
    imageViewer->GetRenderer()->AddActor2D(usageTextActor);

    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->Render();

    renderWindowInteractor->Start();

    minSlice = imageViewer->GetSliceMin();
    maxSlice = imageViewer->GetSliceMax();

    ui->sliderSlices->setMinimum(minSlice);
    ui->sliderSlices->setMaximum(maxSlice);
    ui->labelSlicesNumber->setText(QString::number(maxSlice - minSlice));
    ui->labelFolderName->setText(QString::fromStdString(folderDICOM));

    QFont fontB("Arial", 12, QFont::Bold);
    QFont font("Arial", 12);
    ui->label->setFont(fontB);
    ui->label_2->setFont(font);
    ui->label->setText("DICOM series is successfully loaded.");
    ui->label_2->setText(" Use mouse to scroll through slices.\n Go to slice where tumor can be seen properly.\n Then click on 'Select Seed Point' box and \n put 1st seed on center of liver \n (make sure its not vessel) \n and now 2nd seed on center of tumour part.\n After selecting seed points,\n click on 'Liver Segmentation' box \n & please wait.");

}

void Segmenter::on_buttonOpenFolder_clicked() {
    openDICOMFolder();
}

void Segmenter::on_sliderSlices_sliderMoved(int position) {
    imageViewer->SetSlice(position);
    imageViewer->Render();
}

void Segmenter::on_SeedPoint_clicked() {

    // An interactor
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Create the representation
    vtkSmartPointer<vtkPointHandleRepresentation2D> handle =
        vtkSmartPointer<vtkPointHandleRepresentation2D>::New();
    handle->GetProperty()->SetColor(1, 0, 0);
    vtkSmartPointer<vtkSeedRepresentation> rep =
        vtkSmartPointer<vtkSeedRepresentation>::New();
    rep->SetHandleRepresentation(handle);

    // Seed widget
    vtkSmartPointer<vtkSeedWidget> seedWidget =
        vtkSmartPointer<vtkSeedWidget>::New();
    seedWidget->SetInteractor(renderWindowInteractor);
    seedWidget->SetRepresentation(rep);

    vtkSmartPointer<vtkSeedCallback> seedCallback =
        vtkSmartPointer<vtkSeedCallback>::New();
    seedCallback->SetRepresentation(rep);
    seedWidget->AddObserver(vtkCommand::PlacePointEvent, seedCallback);
    seedWidget->AddObserver(vtkCommand::InteractionEvent, seedCallback);

    ui->widget->SetRenderWindow(renderWindow);

    renderWindow->Render();
    seedWidget->On();
    renderWindowInteractor->Start();

}

int Segmenter::on_LiverSegmentation_clicked()
{
    ui->label->setText("Processing the Liver segmentation. Please wait..");

    smoothing1->SetTimeStep(0.015);
    smoothing1->SetNumberOfIterations(5);
    smoothing1->SetConductanceParameter(9.0);
    smoothing1->SetInput(readeritk->GetOutput());

    const InputPixelType timeThreshold = 250;

    liverThresholdFilter->SetLowerThreshold(0.0);
    liverThresholdFilter->SetUpperThreshold(timeThreshold);
    liverThresholdFilter->SetInsideValue(255);
    liverThresholdFilter->SetOutsideValue(0);

    const double sigma = 1.5;
    gradientMagnitude1->SetSigma(sigma);

    sigmoid1->SetOutputMinimum(0.0);
    sigmoid1->SetOutputMaximum(1.0);

    const double alpha = -0.5;
    const double beta = 3.0;

    sigmoid1->SetAlpha(alpha);
    sigmoid1->SetBeta(beta);

    double x = ((seedCoords[0][0] - 76) * 1.4);
    double y = ((seedCoords[0][1] - 76) * 1.4);

    std::cout << "x value: " << x << std::endl;
    std::cout << "y value: " << y << std::endl;
    std::cout << "z value: " << seedCoords[0][2] << std::endl;

    InputImageType::IndexType  seedPosition;
    seedPosition[0] = x; 
    seedPosition[1] = (512 - y); 
    seedPosition[2] = seedCoords[0][2]; 

    std::cout << "1st coordinate of liver: " << seedPosition[0] << std::endl;
    std::cout << "2nd coordinate of liver: " << seedPosition[1] << std::endl;
    std::cout << "3rd coordinate of liver: " << seedPosition[2] << std::endl;

    NodeType node;
    constexpr double seedValue = 0.0;
    node.SetValue(seedValue);
    node.SetIndex(seedPosition);
  
    seeds1->Initialize();
    seeds1->InsertElement(0, node);

    fastMarching1->SetTrialPoints(seeds1);

    fastMarching1->SetOutputSize(
        readeritk->GetOutput()->GetBufferedRegion().GetSize());

    const double stoppingTime = 250;
    fastMarching1->SetStoppingValue(stoppingTime);
   
    gradientMagnitude1->SetInput(smoothing1->GetOutput());
    sigmoid1->SetInput(gradientMagnitude1->GetOutput());
    fastMarching1->SetInput(sigmoid1->GetOutput());
    liverThresholdFilter->SetInput(fastMarching1->GetOutput());
    caster->SetInput(liverThresholdFilter->GetOutput());

    typedef itk::ImageToVTKImageFilter<OutputImageType> ConnectorType;
    ConnectorType::Pointer connector1 = ConnectorType::New();
    connector1->SetInput(caster->GetOutput());

    //Exceptional handling
    try
    {
        connector1->Update();
    }
    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
   
    vtkSmartPointer<vtkImageData> volumeLiver = vtkSmartPointer<vtkImageData>::New();
    volumeLiver->DeepCopy(connector1->GetOutput());

    vtkSmartPointer<vtkImageFlip> flipXFilter1 =
        vtkSmartPointer<vtkImageFlip>::New();
    flipXFilter1->SetFilteredAxis(1); // flip x axis
    flipXFilter1->SetInputData(volumeLiver);
    flipXFilter1->Update();

    imageViewer->SetInputConnection(flipXFilter1->GetOutputPort());
    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->Render();

    renderWindowInteractor->Start();

    ui->label->setText("Result of Liver segmentation is completed.");
    ui->label_2->setText("Now click on 'Tumor segmentation' box.");
}

int Segmenter::on_TumorSegmentation_clicked()
{
    ui->label->setText("Processing the Tumor segmentation. Please wait..");

    smoothing2->SetTimeStep(0.015);
    smoothing2->SetNumberOfIterations(5);
    smoothing2->SetConductanceParameter(9.0);
    smoothing2->SetInput(readeritk->GetOutput());

    const InputPixelType timeThreshold = 20;

    tumourThresholdFilter->SetLowerThreshold(0.0);
    tumourThresholdFilter->SetUpperThreshold(timeThreshold);
    tumourThresholdFilter->SetInsideValue(255);
    tumourThresholdFilter->SetOutsideValue(0);

    const double sigma = 1.5;
    gradientMagnitude2->SetSigma(sigma);

    sigmoid2->SetOutputMinimum(0.0);
    sigmoid2->SetOutputMaximum(1.0);

    const double alpha = -0.5;
    const double beta = 3.0;

    sigmoid2->SetAlpha(alpha);
    sigmoid2->SetBeta(beta);

    double x = ((seedCoords[1][0] - 76) * 1.4);
    double y = ((seedCoords[1][1] - 76) * 1.4);

    std::cout << "x value: " << x << std::endl;
    std::cout << "y value: " << y << std::endl;
    std::cout << "z value: " << seedCoords[1][2] << std::endl;

    InputImageType::IndexType  seedPosition;
    seedPosition[0] = x; 
    seedPosition[1] = (512 - y); 
    seedPosition[2] = seedCoords[1][2]; 

    std::cout << "1st coordinate of tumor: " << seedPosition[0] << std::endl;
    std::cout << "2nd coordinate of tumor: " << seedPosition[1] << std::endl;
    std::cout << "3rd coordinate of tumor: " << seedPosition[2] << std::endl;

    NodeType node;
    constexpr double seedValue = 0.0;
    node.SetValue(seedValue);
    node.SetIndex(seedPosition);

    seeds2->Initialize();
    seeds2->InsertElement(0, node);
   
    fastMarching2->SetTrialPoints(seeds2);

    fastMarching2->SetOutputSize(
        readeritk->GetOutput()->GetBufferedRegion().GetSize());

    const double stoppingTime = 20;
    fastMarching2->SetStoppingValue(stoppingTime);

    gradientMagnitude2->SetInput(smoothing2->GetOutput());
    sigmoid2->SetInput(gradientMagnitude2->GetOutput());
    fastMarching2->SetInput(sigmoid2->GetOutput());
    tumourThresholdFilter->SetInput(fastMarching2->GetOutput());

    typedef itk::CastImageFilter< InputImageType, OutputImageType >  CastingFilterType;
    CastingFilterType::Pointer caster2 = CastingFilterType::New();

    typedef itk::ImageToVTKImageFilter<OutputImageType>       ConnectorType;
    ConnectorType::Pointer connector2 = ConnectorType::New();

    caster2->SetInput(tumourThresholdFilter->GetOutput());
    connector2->SetInput(caster2->GetOutput());
    try
    {
        connector2->Update();
    }
    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    vtkSmartPointer<vtkImageData> volumeTumour = vtkSmartPointer<vtkImageData>::New();
    volumeTumour->DeepCopy(connector2->GetOutput());

    vtkSmartPointer<vtkImageFlip> flipXFilter2 =
        vtkSmartPointer<vtkImageFlip>::New();
    flipXFilter2->SetFilteredAxis(1); // flip x axis
    flipXFilter2->SetInputData(volumeTumour);
    flipXFilter2->Update();

    imageViewer->SetInputConnection(flipXFilter2->GetOutputPort());
    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->Render();

    ui->label->setText("Result of Tumor segmentation is completed.");
    ui->label_2->setText("Now click on Show 3D Model");
}

void Segmenter::VesselSegmentation(void)
{
    ui->label->setText("Processing the Vessel segmentation. Please wait..");

    smoothing->SetInput(readeritk->GetOutput());
    smoothing->SetNumberOfIterations(5);
    smoothing->SetTimeStep(0.030);
    smoothing->Update();

    const InputPixelType lowerThreshold = 100;
    const InputPixelType upperThreshold = 160;

    vesselThresholdFilter->SetInput(readeritk->GetOutput());

    vesselThresholdFilter->SetLowerThreshold(lowerThreshold);
    vesselThresholdFilter->SetUpperThreshold(upperThreshold);
    vesselThresholdFilter->SetOutsideValue(0);
    vesselThresholdFilter->SetInsideValue(255);
    vesselThresholdFilter->Update();

    StructuringElementType structuringElement;
    StructuringElementType::SizeType newRadius;
    newRadius.Fill(1.0);
    structuringElement.SetRadius(newRadius);
    structuringElement.CreateStructuringElement();
   
    InputPixelType foreground = 255.0;

    openingFilter->SetInput(vesselThresholdFilter->GetOutput());
    openingFilter->SetKernel(structuringElement);
    openingFilter->SetForegroundValue(foreground);
    openingFilter->Update();

    caster3->SetInput(openingFilter->GetOutput());

    typedef itk::ImageToVTKImageFilter<OutputImageType> ConnectorType;
    ConnectorType::Pointer connector3 = ConnectorType::New();
    connector3->SetInput(caster3->GetOutput());

    //Exceptional handling
    try
    {
        connector3->Update();
    }
    catch (itk::ExceptionObject& e)
    {
        std::cerr << "exception in file reader " << std::endl;
        std::cerr << e << std::endl;
        //return EXIT_FAILURE;
    }
    
    vtkSmartPointer<vtkImageData> volumeVessel = vtkSmartPointer<vtkImageData>::New();
    volumeVessel->DeepCopy(connector3->GetOutput());

    vtkSmartPointer<vtkImageFlip> flipXFilter3 =
        vtkSmartPointer<vtkImageFlip>::New();
    flipXFilter3->SetFilteredAxis(1); // flip x axis
    flipXFilter3->SetInputData(volumeVessel);
    flipXFilter3->Update();

    imageViewer->SetInputConnection(flipXFilter3->GetOutputPort());
    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->Render();

    ui->label_2->setText("Result of Vessel segmentation is completed.");
}

void Segmenter::on_ShowModel_clicked()
{
    VesselSegmentation();

    ui->label->setText("Please wait! The 3D Model is loading...");

    typedef itk::ImageToVTKImageFilter<OutputImageType> ConnectorType;

    caster1->SetInput(liverThresholdFilter->GetOutput());
    ConnectorType::Pointer connector1 = ConnectorType::New();
    connector1->SetInput(caster1->GetOutput());
    connector1->Update();

    caster2->SetInput(tumourThresholdFilter->GetOutput());
    ConnectorType::Pointer connector2 = ConnectorType::New();
    connector2->SetInput(caster2->GetOutput());
    connector2->Update();

    caster3->SetInput(openingFilter->GetOutput());
    ConnectorType::Pointer connector3 = ConnectorType::New();
    connector3->SetInput(caster3->GetOutput());
    connector3->Update();

    volume1->DeepCopy(connector1->GetOutput());
    volume2->DeepCopy(connector2->GetOutput());
    volume3->DeepCopy(connector3->GetOutput());

    vtkSmartPointer<vtkImageFlip> flipXFilter1 =
      vtkSmartPointer<vtkImageFlip>::New();
    flipXFilter1->SetFilteredAxis(1); // flip x axis
    flipXFilter1->SetInputData(volume3);
    flipXFilter1->Update();

    double Value = 1.3;

    surface1->SetInputData(volume1);
    surface1->ComputeNormalsOn();
    surface1->SetValue(0, Value);

    surface2->SetInputData(volume2);
    surface2->ComputeNormalsOn();
    surface2->SetValue(0, Value);

    surface3->SetInputData(volume3);
    surface3->ComputeNormalsOn();
    surface3->SetValue(0, Value);

    leftmapper->SetInputConnection(surface1->GetOutputPort());
    leftmapper->ScalarVisibilityOff();
    leftactor->SetMapper(leftmapper);
    leftactor->GetProperty()->SetColor(1, 0.547237, 0.319073);
    leftactor->GetProperty()->SetEdgeVisibility(.9);
    leftactor->GetProperty()->SetEdgeColor(1, 0, 0);
    leftactor->GetProperty()->SetShading(.1);
    leftactor->GetProperty()->SetOpacity(0.3);

    centermapper->SetInputConnection(surface2->GetOutputPort());
    centermapper->ScalarVisibilityOff();
    centeractor->SetMapper(centermapper);
    centeractor->GetProperty()->SetColor(.9, .9, .9);
    centeractor->GetProperty()->SetOpacity(0.5);

    rightmapper->SetInputConnection(surface3->GetOutputPort()); 
    rightmapper->ScalarVisibilityOff(); 
    rightactor->SetMapper(rightmapper);
    rightactor->GetProperty()->SetColor(.8, .3, .3);

    vtkrenderWindow->SetSize(900, 900);

    // And one interactor
    interactor->SetRenderWindow(vtkrenderWindow);
    interactor->SetInteractorStyle(style);

    // Setup both renderers
    vtkrenderWindow->AddRenderer(leftRenderer);
    leftRenderer->SetBackground(0.329412, 0.34902, 0.427451);

    leftRenderer->AddActor(leftactor);
    leftRenderer->AddActor(rightactor);
    leftRenderer->AddActor(centeractor);
    leftRenderer->ResetCamera();

   // vtkrenderWindow->Render();
   // interactor->Start();

    ui->widget_2->SetRenderWindow(renderWindow2.Get());
    ui->widget_2->GetRenderWindow()->AddRenderer(leftRenderer);
    ui->widget_2->GetInteractor();
    ui->widget_2->update();

    QFont fontBold("Arial", 12, QFont::Bold);
    ui->label->setFont(fontBold);
    ui->label->setText(" 3D Model loaded successfully!! ");
    ui->label_2->setText("Results for 3D model are ready.");

}



