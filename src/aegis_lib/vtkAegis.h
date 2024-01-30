#ifndef vtkAegis__
#define vtkAegis__

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


class vtkAegis
{
  public:
  std::unordered_map<std::string, vtkSmartPointer<vtkDoubleArray>> arrays; ; // map of vtkDoubleArrays

  int vtkPointCounter;
  vtkSmartPointer<vtkPoints> vtkpoints;

  vtkAegis(std::string particleTrace);
  void init_Ptrack_root();
  void init_Ptrack_branch(std::string branchName);
  vtkNew<vtkPolyData> new_track(std::string branchName, vtkPoints* vtkpoints, double heatflux);
  void new_vtkArray(std::string arrName, int nComponents);
  void add_vtkArrays(std::string vtk_input_file);
  void write_unstructuredGrid(std::string vtk_input_file, std::string fileName);
  void write_particle_track(std::string branchName, double heatflux);
  void write_multiBlockData(std::string fileName);

  private:
  std::map<std::string, vtkSmartPointer<vtkMultiBlockDataSet>> particleTracks;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockRoot;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockBranch;
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  std::map<std::string, int> multiBlockCounters; // map of various counters for each branch

  bool drawParticleTracks = false;

};














#endif