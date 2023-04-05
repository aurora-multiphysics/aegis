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
    LOG_WARNING << "Flux at centre ssimag1 " << psimag1; // psiaxis in SMARDDA
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

    psiqbdry = psibdry1;
    psiaxis = psimag1;
    set_rsig();
}

// Read 1D array from eqdsk (PRIVATE)
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

// Read 2D array from eqdsk (PRIVATE)
std::vector<std::vector<double>> equData::read_2darray(int nx, int ny, std::string varName)
{
  std::vector<std::vector<double>> work2d(nw, std::vector<double>(nh));
  LOG_WARNING << "Reading " << varName << " values...";
  for (int i=0; i<nx; i++)
  {
    for(int j=0; j<ny; j++)
    {
      eqdsk_file >> work2d[i][j];
      // work2d[i][j] = -work2d[i][j]; // reverse sign of psi if needed
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

// Write singular line out in EQDSK format (PRIVATE) 
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

// Write out eqdsk arrays (PRIVATE)
void equData::eqdsk_write_array(std::ofstream &file, std::vector<double> array, int counter)
{
  double element;
  for(int i=0; i<nw; i++)
  {
    element = array[i];
    counter = eqdsk_line_out(file, element, counter);
  }
}

// Initialise the 1D arrays and 2d spline functions
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
  
  // Alglib::spline2ddiff function can return a value of spline at (R,Z) as well as derivatives. 
  // I.e no need to have a separate spline for each derivative dPsidR and dPsidZ 

} 

// Write out psi(R,Z) data for gnuplotting
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

// Set sign to determine if psi decreasing or increasing away from centre
void equData::set_rsig() 
{
  if (psiqbdry-psiaxis < 0)
  {
    rsig = -1.0;
  }
  else
  {
    rsig = 1.0;
  }

  LOG_INFO << "Value of rsig (psiqbdry-psiaxis) = " << rsig;
}

// Find central psi extrema
void equData::centre()
{
  int igr; // R position index of global extremum
  int igz; // Z position index of global extremum
  int soutr = 10; // maximum number of outer searches
  int sinr = 10; // maximum number of inner searches
  int isr; // R search (increment) direction 
  int isz; // Z search (increment) direction
  int ir; // current R position as index
  int iz; // current Z position as index
  double zpp; // current value of psi
  double zpg; // value of psi at global extremum


  // using central index as guess (See beq_centre for other cases for initial guess)
  igr = (rmax-rmin)/(2*dr);
  igz = (zmax-zmin)/(2*dz); 
  
  //set initially positive search directions for increasing psi
  isr = 1;
  isz = 1;

  for (int j=1; j<=soutr; j++)
  {
    ir = igr;
    iz = igz;
    zpg = spline2dcalc(psiSpline, rmin+(ir-1)*dr, zmin+(iz-1));

    // search in r
    // change direction if necessary
    ir = igr+isr;
    zpp = spline2dcalc(psiSpline, rmin+(ir-1)*dr, zmin+(iz-1)*dz);
    if ((zpp-zpg)*rsig < 0)
    {
      igr = ir;
      zpg = zpp;
    }
    else
    {
      isr = -1;
      ir = igr;
    }

    for (int i=1; i<=sinr; i++)
    {
      ir = ir+isr;
      zpp = spline2dcalc(psiSpline, rmin+(ir-1)*dr, zmin+(iz-1)*dz);
      if ((zpp-zpg)*rsig <0)
      {
        igr = ir;
        zpg = zpp;
      }
      else
      {
        break;
      }
    }

    // search in Z
    // change direction if necessary

    iz = igz+isz;
    zpp = spline2dcalc(psiSpline, rmin+(ir-1)*dr, zmin+(iz-1)*dz);
    if ((zpp-zpg)*rsig < 0)
    {
      igz = iz;
      zpg = zpp;
    }
    else
    {
      isz = -1;
      iz = igz;
    }

    for (int i=1; i<=sinr; i++)
    {
      iz = iz+isz;
      zpp = spline2dcalc(psiSpline, rmin+(ir-1)*dr, zmin+(iz-1)*dz);
      if ((zpp-zpg)*rsig < 0)
      {
        igz = iz;
        zpg = zpp;
      }
      else
      {
        break;
      }
    }

    // check global minimum
    zpp = spline2dcalc(psiSpline, rmin+(igr-2)*dr, zmin+(igz-1)*dz);
    if ((zpp-zpg)*rsig < 0 ) {continue;}
    zpp = spline2dcalc(psiSpline, rmin+igr*dr, zmin+(igz-1)*dz);
    if ((zpp-zpg)*rsig < 0 ) {continue;}
    zpp = spline2dcalc(psiSpline, rmin+(igr-1)*dr, zmin+(igz-2)*dz);
    if ((zpp-zpg)*rsig < 0) {continue;}
    zpp = spline2dcalc(psiSpline, rmin+(igr-1)*dr, zmin+igz*dz);
    if ((zpp-zpg)*rsig < 0) {continue;}
    break;
  }
  rcen = rmin+(igr-1)*dr;
  zcen = zmin+(igz-1)*dz;

  LOG_WARNING << "Rcen value calculated from equData.centre() = " << rcen;
  LOG_WARNING << "Zcen value calculated from equData.centre() = " << zcen;

}

// Create 2d spline structures for R(psi,theta) and Z(psi,theta)
void equData::rz_splines()
{
  int iext; // maximum allowed number of knots for spline in r
  int iknot; // actual number of knots for splinie in r
  int intheta; // actual number of angles defining R,Z(psi,theta)
  int intv; // interval in which spline inverse found 
  double zsig; // sign of dpsi/dr value (same as rsig)
  int isig; // zero if psi decreases outward, else unity
  double zdsrmin; // floor to Delta r_i
  double cpsi; // constant for estimating Delta r_i
  double theta; // theta_j
  double tmin; // minimum theta for R,Z(psi,theta)
  double dpdr; // dpsi/dr
  double dpdz; // dpsi/dz
  double sr; // r
  double costheta; // cos(theta_j)
  double sintheta; // sin(theta_j)
  double re; // R_i
  double ze; // Z_i
  double dpdsr; // {dspi/dr}_{i}
  double dsr; // Delta r_i
  double dpdsrl; // {dspi/dr}_{i-1}
  double psi1; // psi
  double psi2; // psi
  
  

}