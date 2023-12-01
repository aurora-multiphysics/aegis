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

void vtkAegis::init_Ptrack_branch(const char* branchName, std::map<std::string, vtkNew<vtkMultiBlockDataSet>> vtkParticleTracks)
{

  if (vtkParticleTracks.find(branchName) == vtkParticleTracks.end())
  {
    int staticCast = multiBlockCounters.size();
    std::cout << "STATICAST INT = " << staticCast << std::endl;
    multiBlockBranch->SetBlock(0, vtkParticleTracks[branchName]); // set block 
    multiBlockBranch->GetMetaData(static_cast<int>(staticCast)) // name block
                   ->Set(vtkCompositeDataSet::NAME(), branchName);
    multiBlockCounters[branchName] = 0;
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