#include "vtkAegis.h"
#include "simpleLogger.h"


vtkAegis::vtkAegis(std::string particleTrace)
{
  this->unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  if (particleTrace == "yes") {
    this->drawParticleTracks = true;
  }
}

void vtkAegis::init_Ptrack_root()
{
  // Initalise root
    this->multiBlockRoot = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    this->multiBlockBranch = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    this->multiBlockRoot->SetBlock(0, this->multiBlockBranch); // set block 
    this->multiBlockRoot->GetMetaData(static_cast<int>(0)) // name block
                  ->Set(vtkCompositeDataSet::NAME(), "Particle Tracks");
    LOG_INFO << "Initialising particle_tracks root ";
}

void vtkAegis::init_Ptrack_branch(std::string branchName)
{
  if (this->particleTracks.find(branchName) == this->particleTracks.end())
  {
    this->particleTracks[branchName] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    int staticCast = this->multiBlockCounters.size();
    multiBlockBranch->SetBlock(staticCast, this->particleTracks[branchName]); // set block 
    multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                   ->Set(vtkCompositeDataSet::NAME(), branchName);
    this->multiBlockCounters[branchName] = 0;
    std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchName << std::endl;
  }
}

vtkNew<vtkPolyData> vtkAegis::new_track(std::string branchName, vtkPoints* vtkpoints, double heatflux)
{
  vtkNew<vtkPolyLine> vtkpolyline;
  int nVTKPts = this->vtkpoints->GetNumberOfPoints();
  vtkpolyline->GetPointIds()->SetNumberOfIds(nVTKPts);
  for (unsigned int i=0; i<nVTKPts; i++)
  {
    vtkpolyline->GetPointIds()->SetId(i, i);
  }
  
  vtkNew<vtkPolyData> vtkPolydata;
  vtkNew<vtkCellArray> vtkcellarray;
  vtkcellarray->InsertNextCell(vtkpolyline);

  vtkSmartPointer<vtkPolyData> vtkpolydata;
  vtkPolydata->SetPoints(this->vtkpoints);
  vtkPolydata->SetLines(vtkcellarray);

  vtkNew<vtkDoubleArray> vtkHeatflux;
  vtkHeatflux->SetNumberOfComponents(1);
  vtkHeatflux->SetName("Heat_Flux");
  vtkHeatflux->InsertNextTuple1(heatflux);
  vtkPolydata->GetFieldData()->AddArray(vtkHeatflux);
  this->multiBlockCounters[branchName] +=1;

  return vtkPolydata;

}



void vtkAegis::new_vtkArray(std::string arrName, int nComponents)
{
  vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New(); 
  tempArray->SetNumberOfComponents(nComponents);
  tempArray->SetName(arrName.data());
  this->arrays.insert(std::make_pair(arrName.data(), tempArray));
  LOG_INFO << "Initialised new vtkDoubleArray '" << arrName  << "' with nComponents = " << nComponents;
}

void vtkAegis::add_vtkArrays(std::string vtk_input_file) // read stl and add arrays
{
  // Read in STL file 
  vtkNew<vtkSTLReader> vtkstlReader; // STL reader 
  vtkstlReader->SetFileName(vtk_input_file.data());
  vtkstlReader->Update();
  
  LOG_INFO << "Initialising vtkUnstructuredGrid... ";


  // Transform PolyData to vtkUnstructuredGrid datatype using append filter
  vtkNew<vtkAppendFilter> appendFilter;
  vtkPolyData* vtkTargetPD = vtkstlReader->GetOutput(); 
  appendFilter->AddInputData(vtkTargetPD);
  appendFilter->Update();
  this->unstructuredGrid->ShallowCopy(appendFilter->GetOutput());

  for (const auto &arr: this->arrays)
  {
    this->unstructuredGrid->GetCellData()->AddArray(arr.second);
    //LOG_INFO << "Added array '" << arr.second << "' to vtkUnstructuredGrid";
  }
}

void vtkAegis::write_unstructuredGrid(std::string vtk_input_file, std::string fileName)
{
  this->add_vtkArrays(vtk_input_file);
  vtkNew<vtkUnstructuredGridWriter> vtkUstrWriter;
  vtkUstrWriter->SetFileName(fileName.data());
  vtkUstrWriter->SetInputData(this->unstructuredGrid);
  vtkUstrWriter->Write();
}

void vtkAegis::write_particle_track(std::string branchName, double heatflux){
  
  if (!this->drawParticleTracks) {return;} // early return if drawing particle tracks disabled
  else{
    this->init_Ptrack_branch(branchName);
    vtkNew<vtkPolyData> polydataTrack;
    polydataTrack = this->new_track(branchName, this->vtkpoints, heatflux);
    this->particleTracks[branchName]->SetBlock(multiBlockCounters[branchName], polydataTrack);
  }
}

void vtkAegis::write_multiBlockData(std::string fileName){
  vtkNew<vtkXMLMultiBlockDataWriter> vtkMBWriter;
  vtkMBWriter->SetFileName("particle_tracks.vtm");
  vtkMBWriter->SetInputData(this->multiBlockRoot);
  vtkMBWriter->Write();
}
