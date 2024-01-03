#include "vtkAegis.h"
#include "simpleLogger.h"


vtkAegis::vtkAegis()
{
  vtkstlReader = vtkSmartPointer<vtkSTLReader>::New();
  vtkTargetUstr = vtkSmartPointer<vtkUnstructuredGrid>::New();
  multiBlockRoot = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  multiBlockBranch = vtkSmartPointer<vtkMultiBlockDataSet>::New();

}

void vtkAegis::init_Ptrack_root()
{
  // Initalise root

    multiBlockRoot->SetBlock(0, multiBlockBranch); // set block 
    multiBlockRoot->GetMetaData(static_cast<int>(0)) // name block
                  ->Set(vtkCompositeDataSet::NAME(), "Particle Tracks");
    LOG_INFO << "Initialising particle_tracks root ";
}

void vtkAegis::init_Ptrack_branch(const char* branchName)
{

  if (particleTracks.find(branchName) == particleTracks.end())
  {
    int staticCast = multiBlockCounters.size();
    multiBlockBranch->SetBlock(staticCast, particleTracks[branchName]); // set block 
    multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                   ->Set(vtkCompositeDataSet::NAME(), branchName);
    multiBlockCounters[branchName] = 0;
    std::cout << "vtkMultiBlock Particle_track Branch Initialised - " << branchName << std::endl;
    particleTracks[branchName] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  }
}

vtkNew<vtkPolyData> vtkAegis::new_track(const char* branchName, vtkPoints* vtkpoints, double heatflux)
{
  vtkNew<vtkPolyLine> vtkpolyline;
  int nVTKPts = vtkpoints->GetNumberOfPoints();
  vtkpolyline->GetPointIds()->SetNumberOfIds(nVTKPts);
  for (unsigned int i=0; i<nVTKPts; i++)
  {
    vtkpolyline->GetPointIds()->SetId(i, i);
  }
  
  vtkNew<vtkPolyData> vtkPolydata;
  vtkNew<vtkCellArray> vtkcellarray;
  vtkcellarray->InsertNextCell(vtkpolyline);

  vtkSmartPointer<vtkPolyData> vtkpolydata;
  vtkPolydata->SetPoints(vtkpoints);
  vtkPolydata->SetLines(vtkcellarray);

  vtkNew<vtkDoubleArray> vtkHeatflux;
  vtkHeatflux->SetNumberOfComponents(1);
  vtkHeatflux->SetName("Heat_Flux");
  vtkHeatflux->InsertNextTuple1(heatflux);
  vtkPolydata->GetFieldData()->AddArray(vtkHeatflux);
  multiBlockCounters[branchName] +=1;

  return vtkPolydata;

}

void vtkAegis::add_track(const char* branchName, vtkPoints* vtkpoints, double heatflux)
{
  init_Ptrack_branch(branchName);
  vtkSmartPointer<vtkPolyData> polydataTrack = vtkSmartPointer<vtkPolyData>::New();
  polydataTrack = new_track(branchName, vtkpoints, heatflux);
  particleTracks[branchName]->SetBlock(multiBlockCounters[branchName], polydataTrack);
  arrays["Q"]->InsertNextTuple1(heatflux);
}

void vtkAegis::new_vtkArray(std::string arrName, int nComponents)
{
  vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New(); 
  tempArray->SetNumberOfComponents(nComponents);
  tempArray->SetName(arrName.data());
  arrays.insert(std::make_pair(arrName.data(), tempArray));
  LOG_INFO << "Initialised new vtkDoubleArray '" << arrName  << "' with nComponents = " << nComponents;

}



