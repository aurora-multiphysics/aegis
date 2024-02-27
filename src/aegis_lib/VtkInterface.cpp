#include "VtkInterface.h"
#include "SimpleLogger.h"



VtkInterface::VtkInterface(const std::shared_ptr<InputJSON> &inputs)
{
  json vtkNamelist;
  unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (inputs->data.contains("vtk_params"))
  {
    vtkNamelist = inputs->data["vtk_params"];
    drawParticleTracks = vtkNamelist["draw_particle_tracks"];
    vtk_input_file = vtkNamelist["VTK"];
  }

}

void VtkInterface::init_Ptrack_root()
{
  // Initalise root
    multiBlockRoot = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    multiBlockBranch = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    multiBlockRoot->SetBlock(0, multiBlockBranch); // set block 
    multiBlockRoot->GetMetaData(static_cast<int>(0)) // name block
                  ->Set(vtkCompositeDataSet::NAME(), "Particle Tracks");
}

void VtkInterface::init_Ptrack_branch(std::string branchName)
{
  if (particleTracks.find(branchName) == particleTracks.end())
  {
    particleTracks[branchName] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    int staticCast = multiBlockCounters.size();
    multiBlockBranch->SetBlock(staticCast, particleTracks[branchName]); // set block 
    multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                   ->Set(vtkCompositeDataSet::NAME(), branchName);
    multiBlockCounters[branchName] = 0;
  }
}

void VtkInterface::init(){
  
  if (drawParticleTracks){
    init_Ptrack_root();
  }

  new_vtkArray("Q", 1);
  // new_vtkArray("B.n_direction", 1);
  // new_vtkArray("Normal", 3);
  // new_vtkArray("B_field", 3);
  // new_vtkArray("Psi_Start", 1);
  // new_vtkArray("B.n", 1);

  
  if (drawParticleTracks && rank == 0){
    std::cout << "vtkMultiBlockDataSet initialised for particle tracks" << std::endl; 
  }
}

vtkNew<vtkPolyData> VtkInterface::new_track(std::string branchName, vtkSmartPointer<vtkPoints> points, double heatflux)
{
  vtkNew<vtkPolyLine> vtkpolyline;
  int nVTKPts = particleTrackPoints->GetNumberOfPoints();
  vtkpolyline->GetPointIds()->SetNumberOfIds(nVTKPts);
  for (unsigned int i=0; i<nVTKPts; i++)
  {
    vtkpolyline->GetPointIds()->SetId(i, i);
  }
  
  vtkNew<vtkPolyData> vtkPolydata;
  vtkNew<vtkCellArray> vtkcellarray;
  vtkcellarray->InsertNextCell(vtkpolyline);

  vtkSmartPointer<vtkPolyData> vtkpolydata;
  vtkPolydata->SetPoints(particleTrackPoints);
  vtkPolydata->SetLines(vtkcellarray);

  vtkNew<vtkDoubleArray> vtkHeatflux;
  vtkHeatflux->SetNumberOfComponents(1);
  vtkHeatflux->SetName("Heat_Flux");
  vtkHeatflux->InsertNextTuple1(heatflux);
  vtkPolydata->GetFieldData()->AddArray(vtkHeatflux);
  multiBlockCounters[branchName] +=1;

  return vtkPolydata;

}



void VtkInterface::new_vtkArray(std::string arrName, int nComponents)
{
  vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New(); 
  tempArray->SetNumberOfComponents(nComponents);
  tempArray->SetName(arrName.data());
  arrays.insert(std::make_pair(arrName.data(), tempArray));
  std::stringstream newVtkArrayOut;
  newVtkArrayOut << "Initialised new vtkDoubleArray '" << arrName  << "' with nComponents = " << nComponents;
  log_string(LogLevel::INFO, newVtkArrayOut.str());
}

void VtkInterface::add_vtkArrays() // read stl and add arrays
{
  // Read in STL file 
  vtkNew<vtkSTLReader> vtkstlReader; // STL reader 
  vtkstlReader->SetFileName(vtk_input_file.data());
  vtkstlReader->Update();
  
  log_string(LogLevel::INFO,"Initialising vtkUnstructuredGrid...");


  // Transform PolyData to vtkUnstructuredGrid datatype using append filter
  vtkNew<vtkAppendFilter> appendFilter;
  vtkPolyData* vtkTargetPD = vtkstlReader->GetOutput(); 
  appendFilter->AddInputData(vtkTargetPD);
  appendFilter->Update();
  unstructuredGrid->ShallowCopy(appendFilter->GetOutput());

  for (const auto &arr: arrays)
  {
    unstructuredGrid->GetCellData()->AddArray(arr.second);
  }
}


void VtkInterface::write_unstructuredGrid(std::string fileName)
{
  add_vtkArrays();
  vtkNew<vtkUnstructuredGridWriter> vtkUstrWriter;
  vtkUstrWriter->SetFileName(fileName.data());
  vtkUstrWriter->SetInputData(unstructuredGrid);
  vtkUstrWriter->Write();
}

void VtkInterface::write_particle_track(std::string branchName, double heatflux){
  
  if (!drawParticleTracks) {return;} // early return if drawing particle tracks disabled
  else{
    init_Ptrack_branch(branchName);
    vtkNew<vtkPolyData> polydataTrack;
    polydataTrack = new_track(branchName, particleTrackPoints, heatflux);
    particleTracks[branchName]->SetBlock(multiBlockCounters[branchName], polydataTrack);
  }
}

void VtkInterface::write_multiBlockData(std::string fileName){
  if (!drawParticleTracks) {return;} // early return if drawing particle tracks disabled
  else
  {
    vtkNew<vtkXMLMultiBlockDataWriter> vtkMBWriter;
    vtkMBWriter->SetFileName("particle_tracks.vtm");
    vtkMBWriter->SetInputData(multiBlockRoot);
    vtkMBWriter->Write();
  }
}

void VtkInterface::insert_next_uStrGrid(std::string arrayName, std::vector<double> valuesToAdd){
  if (valuesToAdd.size() > 3) 
  { 
    return;
  }
  arrays[arrayName]->InsertNextTuple3(valuesToAdd[0], valuesToAdd[1], valuesToAdd[2]);
}

void VtkInterface::insert_next_uStrGrid(std::string arrayName, double valueToAdd){
  arrays[arrayName]->InsertNextTuple1(valueToAdd);
}

void VtkInterface::insert_zero_uStrGrid()
{
  for (const auto &i:arrays) // loop through all arrays and insert next value as 0.0
  {
    // check number of components in array before inserting 0
    if (i.second->GetSize() == 3) {i.second->InsertNextTuple3(0.0, 0.0, 0.0);}
    else if (i.second->GetSize() == 1) {i.second->InsertNextTuple1(0.0);}
  }
}

void VtkInterface::init_new_vtkPoints(){
  particleTrackPoints = vtkSmartPointer<vtkPoints>::New();
}

void VtkInterface::insert_next_point_in_track(std::vector<double> pointsToAdd){
  particleTrackPoints->InsertNextPoint(pointsToAdd[0], pointsToAdd[1], pointsToAdd[2]);
}

void VtkInterface::mpi_write_uStrGrid(std::string vtk_input_file, std::vector<double> heatfluxVector){



}