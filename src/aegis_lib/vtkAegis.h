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
  vtkSmartPointer<vtkSTLReader> vtkstlReader; // STL reader 
  vtkSmartPointer<vtkUnstructuredGrid> vtkTargetUstr; // Unstructured grib for heat flux map
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockRoot; // multiblock data strucure root
  vtkSmartPointer<vtkMultiBlockDataSet> multiBlockBranch; // branches for multiblock data structure
  std::map<std::string, int> multiBlockCounters; // map of various counters for each branch
  std::unordered_map<std::string, vtkSmartPointer<vtkDoubleArray>> arrays; ; // map of vtkDoubleArrays
  std::map<std::string, vtkSmartPointer<vtkMultiBlockDataSet>> particleTracks;

  vtkAegis();
  void init_Ptrack_root();
  void init_Ptrack_branch(const char* branchName);
  vtkNew<vtkPolyData> new_track(const char* branchName, vtkPoints* vtkpoints, double heatflux);
  void add_track(const char* branchName, vtkPoints* vtkpoints, double heatflux);
  void new_vtkArray(std::string arrName, int nComponents);
};














#endif