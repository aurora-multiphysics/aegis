#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "simpleLogger.h"
#include "equData.h"


#include "alglib/interpolation.h"

// Read eqdsk file
void equData::read_eqdsk(std::string filename)
{
    eqdsk_file.open(filename);
    std::stringstream header_ss;
    std::string temp;
    std::vector<int> header_ints;

    LOG_WARNING << "eqdsk file to be read - " << filename;

    // Extract header information
    std::getline(eqdsk_file, header);
    header_ss << header;
    int number_found;
    
    while (header_ss >> temp)
    {
      if (std::stringstream(temp) >> number_found)
      {
        header_ints.push_back(number_found);
      }
    }

    // Store nw and nh from header information (last two numbers in header)
    nw = header_ints[header_ints.size()-2];
    nh = header_ints[header_ints.size()-1];
    LOG_WARNING << "Number of grid points in R (nw) = " << nw;
    LOG_WARNING << "Number of grid points in Z (nh) = " << nh;

    // Read first four lines of data
    eqdsk_file >> rdim >> zdim >> rcentr >> rgrid >> zmid;
    eqdsk_file >> rmaxis >> zmaxis >> psimag1 >> psibdry1 >> bcentr;
    eqdsk_file >> cpasma >> psimag2 >> xdum >> rmaxis >> xdum;
    eqdsk_file >> zmaxis >> xdum >> psibdry2 >> xdum >> xdum;

    LOG_WARNING << "Geometry parameters from EFIT:";
    LOG_WARNING << "Domain size in R rdim " << rdim;
    LOG_WARNING << "Domain size in Z zdim " << zdim;
    LOG_WARNING << "R at centre " << rcentr;
    LOG_WARNING << "Domain start in R rgrid " << rgrid;
    LOG_WARNING << "Domain centre in Z zmid " << zmid;
    LOG_WARNING << "Plasma parameters from EFIT:";
    LOG_WARNING << "B at rcentre " << bcentr;
    LOG_WARNING << "Flux at centre ssimag1 " << psimag1;
    LOG_WARNING << "Flux at boundary ssibry1 "<< psibdry1;
    LOG_WARNING << "Plasma centre in R rmaxis " << rmaxis;
    LOG_WARNING << "Plasma centre in Z zmaxis " << zmaxis;
    LOG_WARNING << "Plasma current " << cpasma;

    // Read 1D array data

    if (nw > 0)
    {
      fpol = read_array(nw, "fpol");
      pres = read_array(nw, "pres");
      ffprime = read_array(nw, "ffprime");
      pprime = read_array(nw, "pprime");
    }
    else
    {
      LOG_FATAL << "Error reading 1D data from eqdsk";
    }

    // Read Psi(R,Z) data
    if (nh > 0)
    {
      psi = read_2darray(nw, nh, "Psi(R,Z)");
    }
    else
    {
      LOG_FATAL << "Error reading 2D data from eqdsk";
    }

    // Read the rest of the data
    qpsi = read_array(nw, "qpsi");
    eqdsk_file >> nbdry >> nlim;
    LOG_WARNING << "Number of boundary points " << nbdry;
    LOG_WARNING << "Number of limiter points " << nlim;

    if (nbdry > 0)
    {
      rbdry.resize(nbdry);
      zbdry.resize(nbdry);
      LOG_WARNING << "Reading rbdry and zbdry...";
      for(int i=0; i<nbdry; i++) // Read in n elements into vector from file
      {
        eqdsk_file >> rbdry[i] >> zbdry[i];
      }
      LOG_WARNING << "Number of rbdry/zbdry values read " << nbdry;
    }
    else
    {
      LOG_FATAL << "Error reading boundary data from eqdsk";
    }

    if (nlim > 0)
    {
      LOG_WARNING << "Reading rlim and zlim...";
      rlim.resize(nlim);
      zlim.resize(nlim);
      for (int i=0; i<nlim; i++)
      {
        eqdsk_file >> rlim[i] >> zlim[i];
      } 
      LOG_WARNING << "Number of rlim/zlim values read " << nlim;
    }
    else
    {
      LOG_FATAL << "Error reading limiter data from eqdsk";
    }
}

// Read 1D array from eqdsk
std::vector<double> equData::read_array(int n, std::string varName)
{
  std::vector<double> work1d(n); // working vector of doubles of size n
  LOG_WARNING << "Reading " << varName << " values...";
  for(int i=0; i<nw; i++) // Read in n elements into vector from file
  {
    eqdsk_file >> work1d[i];
  }
  LOG_WARNING << "Number of " << varName << " values read = " << n;
  return work1d;
}

// Read 2D array from eqdsk
std::vector<std::vector<double>> equData::read_2darray(int nx, int ny, std::string varName)
{
  std::vector<std::vector<double>> work2d(nw, std::vector<double>(nh));
  LOG_WARNING << "Reading " << varName << " values...";
  for (int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      eqdsk_file >> work2d[i][j];
      work2d[i][j] = -work2d[i][j];
    }
  }
  LOG_WARNING << "Number of " << varName << " values read = " << nx*ny;
  return work2d;
}

// Write out eqdsk data back out in eqdsk format
void equData::write_eqdsk_out()
{
  std::ofstream eqdsk_out("eqdsk.out");
  eqdsk_out << std::setprecision(9) << std::scientific;
  eqdsk_out << std::setiosflags(std::ios::uppercase);
  double element;
  int counter=0;

  // Write out header information
  eqdsk_out << header << std::endl;

  // Write out initial four lines of floats
  double parameters[20] = {rdim, zdim, rcentr, rgrid, zmid, rmaxis, zmaxis, psimag1, psibdry1, 
                bcentr, cpasma, psimag2, xdum, rmaxis, xdum, zmaxis, xdum, psibdry2, xdum, xdum};

  for (int i=0; i<20; i++)
  {
    element = parameters[i];
    counter = eqdsk_line_out(eqdsk_out, element, counter);
  }

  // Write out data
  eqdsk_write_array(eqdsk_out, fpol, counter); // write fpol array
  eqdsk_write_array(eqdsk_out, pres, counter); // write pres array 
  eqdsk_write_array(eqdsk_out, ffprime, counter); // write ffprime array
  eqdsk_write_array(eqdsk_out, pprime, counter); // write pprime array 

  // Write out psi(R,Z)
  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      element = psi[i][j];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
    }
  }
  if (counter < 5)
  {
    counter = 0;
    eqdsk_out << std::endl;
  }

  eqdsk_write_array(eqdsk_out, qpsi, counter); // write qpsi array
  eqdsk_out << std::endl << "  " << nbdry << "   " << nlim << std::endl; // write nbdry and nlim

  // write rbdry and zbdry
  for (int i=0; i<nbdry; i++)
  { 
      element = rbdry[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = zbdry[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = rbdry[i];
  }
  if (counter < 5)
  {
    counter = 0;
    eqdsk_out << std::endl;
  }

  // write rlim and zlim arrays
  for (int i=0; i<nlim; i++)
  { 
      element = rlim[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = zlim[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
  }

  if (counter < 5)
  {
    counter = 0;
    eqdsk_out << std::endl;
  }
}

// Write singular line out in EQDSK format
int equData::eqdsk_line_out(std::ofstream &file, double element, int counter)
    {
      if (counter == 5)
      {
        file << std::endl;
        counter = 0;
      }
      if (element >= 0)
      {
        file << " " << element;
      }
      else
      {
        file << element;
      }
      counter +=1;
      return counter;
    }

// Write out eqdsk array
void equData::eqdsk_write_array(std::ofstream &file, std::vector<double> array, int counter)
{
  double element;
  for(int i=0; i<nw; i++)
  {
    element = array[i];
    counter = eqdsk_line_out(file, element, counter);
  }
}

// Initialise the 1D arrays and spline functions
void equData::init_interp_splines()
{
  rmin = rgrid;
  zmin = zmid - zdim/2;
  rmax = rgrid+rdim;
  zmax = zmid + zdim/2;

  nr = nw-1;
  nz = nh-1;

  dr = (rmax-rmin)/nr;
  dz = (zmax-zmin)/nz;

  double r_pts[nw];
  double z_pts[nh];
  r_pts[0] = rmin;
  z_pts[0] = zmin;

  for (int i=0; i<nw; i++)
  {
    r_pts[i+1] = r_pts[i]+dr;
  }

  for (int i=0; i<nh; i++)
  {
    z_pts[i+1] = z_pts[i]+dz;
  }

  double psi_pts[nw*nh];
  int count = 0;
  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      psi_pts[count] = psi[i][j];
      count += 1;
    }
  }

  // set 1d arrays for R grid, Z grid and Psi grid
  r_grid.setcontent(nw, r_pts);
  z_grid.setcontent(nh, z_pts);
  psi_grid.setcontent(count, psi_pts);

  // Construct the spline interpolant to be used later 
  alglib::spline2dbuildbilinearv(r_grid, nw, z_grid, nh, psi_grid, 1, psiSpline);

  
} 

void equData::gnuplot_out()
{
  std::ofstream psiRZ_out;
  psiRZ_out.open("psi_RZ.gnu"); 

  for (int j=0; j<nh; j++)
  {
    for (int i=0; i<nw; i++)
    {
      psiRZ_out << i << " " << j << " " <<  r_grid[i] << " " << z_grid[j] << " " << 
                  spline2dcalc(psiSpline, r_grid[i], z_grid[j]) << std::endl;
    }
    psiRZ_out << std::endl;
  }

}