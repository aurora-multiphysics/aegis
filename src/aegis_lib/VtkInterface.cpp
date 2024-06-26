#include "VtkInterface.h"
#include "SimpleLogger.h"

VtkInterface::VtkInterface(const std::shared_ptr<JsonHandler> & inputs)
{
  nlohmann::json vtkNamelist;

  if (inputs->data().contains("vtk_params"))
  {
    vtkNamelist = inputs->data()["vtk_params"];
    drawParticleTracks = vtkNamelist["draw_particle_tracks"];
  }

  if (drawParticleTracks)
  {
    init_Ptrack_root();

    // log_string(LogLevel::WARNING, "vtkMultiBlockDataSet initialised for particle tracks");
  }
}

void
VtkInterface::init_Ptrack_root()
{
  // Initalise root
  multiBlockRoot = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  multiBlockBranch = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  multiBlockRoot->SetBlock(0, multiBlockBranch); // set block
  multiBlockRoot
      ->GetMetaData(static_cast<int>(0)) // name block
      ->Set(vtkCompositeDataSet::NAME(), "Particle Tracks");
}

void
VtkInterface::init_Ptrack_branch(std::string branchName)
{
  if (particleTracks.find(branchName) == particleTracks.end())
  {
    particleTracks[branchName] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    int staticCast = multiBlockCounters.size();
    multiBlockBranch->SetBlock(staticCast, particleTracks[branchName]); // set block
    multiBlockBranch
        ->GetMetaData(static_cast<int>(staticCast)) // name block
        ->Set(vtkCompositeDataSet::NAME(), branchName);
    multiBlockCounters[branchName] = 0;
  }
}

void
VtkInterface::init()
{
}

vtkNew<vtkPolyData>
VtkInterface::new_track(std::string branchName, vtkSmartPointer<vtkPoints> points, double heatflux)
{
  vtkNew<vtkPolyLine> vtkpolyline;
  int nVTKPts = particleTrackPoints->GetNumberOfPoints();
  vtkpolyline->GetPointIds()->SetNumberOfIds(nVTKPts);
  for (unsigned int i = 0; i < nVTKPts; i++)
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
  multiBlockCounters[branchName] += 1;

  return vtkPolydata;
}

void
VtkInterface::write_particle_track(std::string branchName, double heatflux)
{

  if (drawParticleTracks)
  {
    init_Ptrack_branch(branchName);
    vtkNew<vtkPolyData> polydataTrack;
    polydataTrack = new_track(branchName, particleTrackPoints, heatflux);
    particleTracks[branchName]->SetBlock(multiBlockCounters[branchName], polydataTrack);
  }
}

void
VtkInterface::write_multiBlockData(std::string fileName)
{
  if (!drawParticleTracks)
  {
    return;
  } // early return if drawing particle tracks disabled
  else
  {
    vtkNew<vtkXMLMultiBlockDataWriter> vtkMBWriter;
    vtkMBWriter->SetFileName("particle_tracks.vtm");
    vtkMBWriter->SetInputData(multiBlockRoot);
    vtkMBWriter->Write();
  }
}

void
VtkInterface::init_new_vtkPoints()
{
  particleTrackPoints = vtkSmartPointer<vtkPoints>::New();
}

void
VtkInterface::insert_next_point_in_track(std::vector<double> pointsToAdd)
{
  if (drawParticleTracks == true)
  {
    particleTrackPoints->InsertNextPoint(pointsToAdd[0], pointsToAdd[1], pointsToAdd[2]);
  }
}
