#ifndef VtkInterface__
#define VtkInterface__

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyLine.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkCompositeDataSet.h>
#include <vtkInformation.h>
#include <vtkSTLReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkAppendFilter.h>
#include "Inputs.h"
#include "AegisBase.h"

class VtkInterface : public AegisBase
{
  public:
  VtkInterface(const std::shared_ptr<InputJSON> &JSONsettings);
  void init_Ptrack_root();
  void init_Ptrack_branch(std::string branchName);
  void init(); // initialise unstructured grid and particle tracks multiblock
  vtkNew<vtkPolyData> new_track(std::string branchName, vtkSmartPointer<vtkPoints> points, double heatflux);
  void write_particle_track(std::string branchName, double heatflux);
  void write_multiBlockData(std::string fileName);
  void init_new_vtkPoints();
  void insert_next_point_in_track(std::vector<double> pointsToAdd);  


  private:
  int vtkPointCounter;
  std::map<std::string, vtkSmartPointer<vtkMultiBlockDataSet>> particleTracks;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockRoot;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockBranch;
  std::map<std::string, int> multiBlockCounters; // map of various counters for each branch
  vtkSmartPointer<vtkPoints> particleTrackPoints;
  std::vector<std::array<double,3>> points;
  bool drawParticleTracks = false;

};














#endif