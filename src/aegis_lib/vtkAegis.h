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
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  std::map<std::string, int> multiBlockCounters; // map of various counters for each branch
  std::unordered_map<std::string, vtkSmartPointer<vtkDoubleArray>> arrays; ; // map of vtkDoubleArrays
  std::map<std::string, vtkSmartPointer<vtkMultiBlockDataSet>> particleTracks;

  vtkAegis();
  void init_Ptrack_root(vtkSmartPointer<vtkMultiBlockDataSet> &multiBlockRoot, vtkSmartPointer<vtkMultiBlockDataSet> &multiBlockBranch);
  void init_Ptrack_branch(const char* branchName, vtkSmartPointer<vtkMultiBlockDataSet> &multiBlockBranch, vtkSmartPointer<vtkMultiBlockDataSet> &track);
  vtkNew<vtkPolyData> new_track(const char* branchName, vtkPoints* vtkpoints, double heatflux);
  void new_vtkArray(std::string arrName, int nComponents);
  void add_vtkArrays(const char* vtk_input_file);
  void write_unstructuredGrid(const char* vtk_input_file, const char* fileName);
};














#endif