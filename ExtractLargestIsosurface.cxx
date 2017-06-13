/*============================================================================================
The following tool was devloped by Ethan Hayns and Claudeson Azurin at the request of IU Health, 
Dr Shiaofen Fang, and Huang Li. The funding for this project was made possible by the NSF and DoD.
Outside libraries were used to complete this project, mainly the VTK library. 
=============================================================================================*/
#include "vtkSmartPointer.h"
#include "vtkStructuredPointsReader.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkTransform.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkNIFTIImageHeader.h"
#include "vtkNIFTIImageReader.h"
#include "vtkNIFTIImageWriter.h"
#include <vtkContourFilter.h>
#include "vtkIndent.h"
#include <vtkDiscreteMarchingCubes.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkPath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>


/*============================================================================================
The following code normalizes the expression value of the data point including the minimum and
maximum of the data set. To avoid negative values, the equations presented are of a piecwise
manner defines the color for each sphere
=============================================================================================*/
static float arr[3]; 
float* colorVal(float X,float min, float max, float divisor){
 //arr[1]=0;x
max=max/divisor;
min=min/divisor;
float med=(max+min)/2;
float theta=(3.1415926535)*exp((-med)+(2*(med)));
float beta=theta*10;
 X=X/(divisor);
 if(X>med){ 
  arr[0]=ceil(((theta)*(X-med))+(1));
  arr[2]=ceil(((-beta)/(med-min)*X));
 }
 else {
  arr[0]=ceil(((beta)/(med-min)*(X-med))+(1));
  arr[2]=ceil(((-theta)*(X-med))+(1));
 }
 return arr;
}
/*===========================================================================================*/

/*============================================================================================
The following code block takes in a normalized brain image from a deceased donor and an image of
a patient from IU Health. These two were needed to normalize the data points so that they fit ont
the patient's brain
==============================================================================================*/
int main(int argc, char* argv[])
{
argv[0] = "brain.nii";
argv[1] = "brain.nii";
argv[2] = "83";
argv[4] = "T1.nii";
argv[5] = "T1.nii";

  const char* fileName = argv[1];
  //float threshold = atof(argv[2]);
  const char* fileName2=argv[4];
  int extractLargest = 1;
  if (argc == 4)
    {
    extractLargest = atoi(argv[3]);
    }
/*==============================================================================================*/

ifstream file("Probes.csv");
bool correctProbe = false;
std::string value;
std::string probeID;
std::string geneLookUp;
cout<<"Input gene ID \n";
cin >>geneLookUp;
int finder = 0;
std::map < std::string, std::string > probeInfo;
std::map < int, std::string > structName;
float geneValue[946];
int i;
for (i=0;i<947;i++){
  geneValue[i]=0;
}
while (file.good()){
//First IF statement facilitates movement to a new line.
if (finder % 7 == 6){
getline(file, value, '\n');
}
//Second IF statement grabs the probe_ID from the matching gene
if(finder % 7 == 0){
probeID = value;
//cout <<value<<"\n";
}
else{
getline(file, value, ',');
if (finder % 7 == 3){
if (geneLookUp == value){

probeInfo.insert(std::pair<std::string, std::string>(probeID, value));
}
}
}
 
finder++;
}
file.close();
std::map<std::string, std::string>::iterator it = probeInfo.begin();
//cout << "probeInfo contains:\n";
for (it = probeInfo.begin(); it != probeInfo.end(); ++it){
//cout << it->first << " => " << it->second << '\n';
}
//cout << "number of probes: " << probeInfo.size() << endl;
//cout << endl;
//cout << "Searching MicroArray, Progress: " << endl;

long int row = 0;
int progress = 0;
finder = 0;
file.open("MicroarrayExpression.csv");

while (file.good()){
//First IF statement facilitates movement to a new line.
if (finder % 947 == 946){
getline(file, value, '\n');
if (correctProbe == true){
geneValue[945] += std::stof(value);
cout<<value<<"\n";
}
 
correctProbe = false;
finder++;
}
//Second IF statement checks to see if the probe ID matches for the gene we're looking for
if (finder % 947 == 0){
getline(file, value, ',');
probeID = value;
for (it = probeInfo.begin(); it != probeInfo.end(); ++it){
if (it->first == probeID){
correctProbe = true;
 cout<< "Found correct ProbeID: " << value << endl;

}
}
if (!correctProbe){
correctProbe = false;
file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//getline(file, value, '\n');
finder += 946;
}
//Percent progress display
row++;
if (row % 586 == 0){
//cout << progress << "%" << endl;
progress++;
}
finder++;
}
 
//Else statement, used when correctProbe to grab the expression data
else{
getline(file, value, ',');
if (correctProbe == true){
value.resize(7);
geneValue[(finder % 947)-1] += std::stof(value);
}
finder++;

}

}
file.close();

/*=============================================================================================
In the following block of code, we search the SampleAnnot file to determine the position of the 
expression values
==============================================================================================*/
file.open("SampleAnnot.csv");
finder = 0;
int xCoord[947];
int yCoord[947];
int zCoord[947];
row = -1;
while (file.good()){
//First IF statement facilitates movement to a new line.
if (finder % 13 == 12){
getline(file, value, '\n');
row++;
finder++;
}
//skips first line
if (row == -1){
getline(file, value, '\n');
row++;
}
else{
getline(file, value, ',');

if (finder % 13 == 7){

xCoord[row] = std::stoi(value);
}
if (finder % 13 == 8){
yCoord[row] = std::stoi(value);
}
if (finder % 13 == 9){
zCoord[row] = std::stoi(value);
}
finder++;
}
 
}
file.close();
/*===========================================================================================*/

/*============================================================================================
In the following block of code we determine the minimum and maximum values in the dataset and we 
assign individual actors to the expression values. We then use the colorVal function to assign 
the color value to each point. 
==============================================================================================*/
std::vector<vtkSmartPointer<vtkActor> > actors;
float minimum=*std::min_element(geneValue,geneValue+946);
float maximum=*std::max_element(geneValue,geneValue+946);
cout<<"min is:"<<minimum<<"\n";
cout<<"max is:"<<maximum<<"\n";
float numProbes=(float)(probeInfo.size());
for(i = 0; i < 946; i++)
  {
  vtkSmartPointer<vtkSphereSource> sphereSource = 
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetCenter(xCoord[i], yCoord[i], zCoord[i]);
  sphereSource->SetRadius(1);

  vtkSmartPointer<vtkPolyDataMapper> mapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(sphereSource->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = 
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->SetScale(.91);

//float geneValueFloat=geneValue[946];
  float *anotherArray=colorVal(geneValue[i],minimum, maximum,numProbes);
  float redColor=anotherArray[0];
  float blueColor=anotherArray[2];
  actor->GetProperty()->SetColor(redColor,0,blueColor);

  actors.push_back(actor);
  }
/*============================================================================================*/

/*=============================================================================================
In the following code block, we set a reference point so the user can easily explore the data 
===============================================================================================*/

vtkSmartPointer<vtkCubeSource> cubeSource = 
  vtkSmartPointer<vtkCubeSource>::New();
cubeSource->SetCenter(0, 0,0);
  cubeSource->SetXLength(5);
  cubeSource->SetYLength(5);
  cubeSource->SetZLength(5);
vtkSmartPointer<vtkPolyDataMapper> mapper3 = 
  vtkSmartPointer<vtkPolyDataMapper>::New();
mapper3->SetInputConnection(cubeSource->GetOutputPort());

vtkSmartPointer<vtkActor> actor3 = 
  vtkSmartPointer<vtkActor>::New();
actor3->SetMapper(mapper3);
actor3->GetProperty()->SetColor(5,0,5);

actors.push_back(actor3);

/*============================================================================================*/
//Manual Data entry for region selection
int startRegion = 1;
int endRegion = 83;
if (startRegion == 0 && endRegion == 0){
cout << "Enter starting Region: ";
cin >> startRegion;
cout << endl << endl;
cout << "Enter ending Region: ";
cin >> endRegion;
cout << endl;
}

/*==============================================================================================
In the following code block, we render the final images. Note that there are certain code mini-
blocks which are commented out. These are for the deceased brain image and were used for the 
approximation of a scalar for the IU brain image. If the user intends to make the approximation
more precise, they can uncomment the code and do so. 
=============================================================================================*/

// Load data
vtkSmartPointer<vtkNIFTIImageReader> reader =
vtkSmartPointer<vtkNIFTIImageReader>::New();
reader->SetFileName(fileName);

/******************************************/

vtkSmartPointer<vtkNIFTIImageReader> reader2 =
vtkSmartPointer<vtkNIFTIImageReader>::New();
reader2->SetFileName(fileName2);

// Create a 3D model using discrete marching cubes
vtkSmartPointer<vtkDiscreteMarchingCubes> mc =
vtkSmartPointer<vtkDiscreteMarchingCubes>::New();
mc->SetInputConnection(reader->GetOutputPort());
mc->ComputeScalarsOn();
mc->ComputeGradientsOn();

/*===========================================*/

//vtkSmartPointer<vtkDiscreteMarchingCubes> mc2 =
//vtkSmartPointer<vtkDiscreteMarchingCubes>::New();
//mc2->SetInputConnection(reader2->GetOutputPort());


int x = 0; //x will number off the contours

//For loop builds the contours based off of the regions of interest
/*
for (int i = startRegion; i < endRegion; i++){
  mc->SetValue(x, i);
  x=x+1;
}
*/
mc->Update();
//mc2->Update();

// To remain largest region
vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter =
  vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
confilter->SetInputConnection(mc->GetOutputPort());
confilter->SetExtractionModeToAllRegions();
//extract every region
//confilter->SetExtractionModeToLargestRegion(); //only model the largest Region
confilter->Update();

/*======================================================*/
//vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter2 =
//  vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
//confilter2->SetInputConnection(mc2->GetOutputPort());
//confilter2->SetExtractionModeToAllRegions();


// Create a mapper
vtkSmartPointer<vtkPolyDataMapper> mapper =
  vtkSmartPointer<vtkPolyDataMapper>::New();
mapper->SetInputConnection(confilter->GetOutputPort());
mapper->ScalarVisibilityOff(); 


/*******************************************/
//vtkSmartPointer<vtkPolyDataMapper> mapper2 =
  //vtkSmartPointer<vtkPolyDataMapper>::New();
//mapper2->SetInputConnection(confilter2->GetOutputPort());

//mapper2->ScalarVisibilityOff(); 
//mapper2->SetScalarRange(startRegion, endRegion);
//mapper2->SetScalarModeToUseCellData();
//mapper2->SetColorModeToMapScalars();
//mapper2->SetLookupTable(hueLut);

//Transform 
vtkSmartPointer<vtkTransform> transform =
  vtkSmartPointer<vtkTransform>::New();
  transform->Translate(-45,15,210);

// Visualize
vtkSmartPointer<vtkActor> actor =
  vtkSmartPointer<vtkActor>::New();
actor->GetProperty()->SetColor(1,9,1);
actor->GetProperty()->SetOpacity(.03);
//actor->SetOrigin(40,210,-70);
actor->RotateX(-90);
//actor->RotateY(30);
//actor->RotateZ(5);
//actor->SetScale(1);
//actor->Translate(5,0,0);
actor->SetUserTransform(transform);
//actor->SetOrigin(actor->GetPosition());
actor->SetMapper(mapper);
//actor->SetScale(1.1);

/****************************************/

//vtkSmartPointer<vtkActor> actor2 =
//  vtkSmartPointer<vtkActor>::New();
//actor2->GetProperty()->SetColor(1,3,2);
//actor2->GetProperty()->SetOpacity(0.5);
//actor2->SetMapper(mapper2);

//actor2->RotateY(20);
//actor2->RotateX(74);
//actor2->RotateZ(27);
//actor2->SetPosition(40,185,-60);
//actor2->SetScale(.91);

vtkSmartPointer<vtkRenderer> renderer =
  vtkSmartPointer<vtkRenderer>::New();
renderer->AddActor(actor);
//renderer->AddActor(actor2);


for(i = 0; i < actors.size(); i++)
  {
  renderer->AddActor(actors[i]);
  }
  

/**************************************/
vtkSmartPointer<vtkRenderer> renderer2 =
  vtkSmartPointer<vtkRenderer>::New();

vtkSmartPointer<vtkRenderWindow> renwin =
  vtkSmartPointer<vtkRenderWindow>::New();
renwin->AddRenderer(renderer);

vtkSmartPointer<vtkRenderWindowInteractor> iren =
  vtkSmartPointer<vtkRenderWindowInteractor>::New();
iren->SetRenderWindow(renwin);
iren->Initialize();
iren->Start();

return EXIT_SUCCESS;
}



 
	

