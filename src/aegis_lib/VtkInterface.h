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

#include <settings.hpp>

class VtkInterface
{
  public:
  VtkInterface(const std::shared_ptr<InputJSON> &JSONsettings);
  void init_Ptrack_root();
  void init_Ptrack_branch(std::string branchName);
  void init(); // initialise unstructured grid and particle tracks multiblock
  vtkNew<vtkPolyData> new_track(std::string branchName, vtkSmartPointer<vtkPoints> points, double heatflux);
  void new_vtkArray(std::string arrName, int nComponents);
  void add_vtkArrays();
  void write_unstructuredGrid(std::string fileName);
  void write_particle_track(std::string branchName, double heatflux);
  void write_multiBlockData(std::string fileName);
  void insert_next_uStrGrid(std::string arrayName, std::vector<double> valuesToAdd);
  void insert_next_uStrGrid(std::string arrayName, double valueToAdd);
  void init_new_vtkPoints();
  void insert_next_point_in_track(std::vector<double> pointsToAdd);  
  void mpi_write_uStrGrid(std::string vtk_input_file, std::vector<double> heatfluxVector);


  private:
  std::unordered_map<std::string, vtkSmartPointer<vtkDoubleArray>> arrays; ; // map of vtkDoubleArrays
  int vtkPointCounter;
  std::map<std::string, vtkSmartPointer<vtkMultiBlockDataSet>> particleTracks;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockRoot;
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockBranch;
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  std::map<std::string, int> multiBlockCounters; // map of various counters for each branch
  vtkSmartPointer<vtkPoints> particleTrackPoints;
  std::string vtk_input_file;
  bool drawParticleTracks = false;

};














#endif