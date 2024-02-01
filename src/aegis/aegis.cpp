#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>  // for setprecision
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>  // for min/max values
#include <set>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <time.h>
#include <any>

#include <vtkCellArray.h>
#include <vtkNew.h>
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

#include "DagMC.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include <moab/OrientedBoxTreeTool.hpp>
#include "settings.hpp"
#include "simpleLogger.h"
#include "equData.h"
#include "source.h"
#include "integrator.h"
#include "coordtfm.h"
#include "alglib/interpolation.h"
#include "vtkAegis.h"
#include "particle.h"
#include "aegisClass.h"


int main() {
  clock_t start = clock();
  
  AegisClass aegis;
  aegis.Execute();

  clock_t end = clock();
  double elapsed = double(end - start)/CLOCKS_PER_SEC;

  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Elapsed Aegis run time = " << elapsed << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;

  return 0;
}


