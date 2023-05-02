#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <vector>


#include "simpleLogger.h"
#include "equData.h"
#include "alglib/interpolation.h"
#include "coordtfm.h"


// Return eqdsk data structure
eqdskData equData::get_eqdsk_struct()
{
  return eqdsk; 
}

// Read eqdsk file
void equData::read_eqdsk(std::string filename)
{
    eqdsk_file.open(filename);
    std::stringstream header_ss;
    std::string temp;
    std::vector<int> header_ints;

    LOG_WARNING << "eqdsk file to be read - " << filename;

    // Extract header information
    std::getline(eqdsk_file, eqdsk.header);
    header_ss << eqdsk.header;
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
    eqdsk_file >> eqdsk.rdim >> eqdsk.zdim >> eqdsk.rcentr >> eqdsk.rgrid >> eqdsk.zmid;
    eqdsk_file >> eqdsk.rmaxis >> eqdsk.zmaxis >> eqdsk.psimag1 >> eqdsk.psibdry1 >> eqdsk.bcentr;
    eqdsk_file >> eqdsk.cpasma >> eqdsk.psimag2 >> eqdsk.xdum >> eqdsk.rmaxis >> eqdsk.xdum;
    eqdsk_file >> eqdsk.zmaxis >> eqdsk.xdum >> eqdsk.psibdry2 >> eqdsk.xdum >> eqdsk.xdum;

    LOG_WARNING << "Geometry parameters from EFIT:";
    LOG_WARNING << "Domain size in R eqdsk.rdim " << eqdsk.rdim;
    LOG_WARNING << "Domain size in Z eqdsk.zdim " << eqdsk.zdim;
    LOG_WARNING << "R at centre " << eqdsk.rcentr;
    LOG_WARNING << "Domain start in R eqdsk.rgrid " << eqdsk.rgrid;
    LOG_WARNING << "Domain centre in Z eqdsk.zmid " << eqdsk.zmid;
    LOG_WARNING << "Plasma parameters from EFIT:";
    LOG_WARNING << "B at rcentre " << eqdsk.bcentr;
    LOG_WARNING << "Flux at centre ssimag1 " << eqdsk.psimag1; // psiaxis in SMARDDA
    LOG_WARNING << "Flux at boundary ssibry1 "<< eqdsk.psibdry1;
    LOG_WARNING << "Plasma centre in R eqdsk.rmaxis " << eqdsk.rmaxis;
    LOG_WARNING << "Plasma centre in Z eqdsk.zmaxis " << eqdsk.zmaxis;
    LOG_WARNING << "Plasma current " << eqdsk.cpasma;

    // Read 1D array data

    if (nw > 0)
    {
      eqdsk.fpol = read_array(nw, "eqdsk.fpol");
      eqdsk.pres = read_array(nw, "eqdsk.pres");
      eqdsk.ffprime = read_array(nw, "eqdsk.ffprime");
      eqdsk.pprime = read_array(nw, "eqdsk.pprime");
    }
    else
    {
      LOG_FATAL << "Error reading 1D data from eqdsk";
    }

    // Read Psi(R,Z) data
    if (nh > 0)
    {
      eqdsk.psi = read_2darray(nw, nh, "Psi(R,Z)");
    }
    else
    {
      LOG_FATAL << "Error reading 2D data from eqdsk";
    }

    // Read the rest of the data
    eqdsk.qpsi = read_array(nw, "eqdsk.qpsi");
    eqdsk_file >> eqdsk.nbdry >> eqdsk.nlim;
    LOG_WARNING << "Number of boundary points " << eqdsk.nbdry;
    LOG_WARNING << "Number of limiter points " << eqdsk.nlim;

    if (eqdsk.nbdry > 0)
    {
      eqdsk.rbdry.resize(eqdsk.nbdry);
      eqdsk.zbdry.resize(eqdsk.nbdry);
      LOG_WARNING << "Reading eqdsk.rbdry and eqdsk.zbdry...";
      for(int i=0; i<eqdsk.nbdry; i++) // Read in n elements into vector from file
      {
        eqdsk_file >> eqdsk.rbdry[i] >> eqdsk.zbdry[i];
      }
      LOG_WARNING << "Number of eqdsk.rbdry/eqdsk.zbdry values read " << eqdsk.nbdry;
    }
    else
    {
      LOG_FATAL << "Error reading boundary data from eqdsk";
    }

    if (eqdsk.nlim > 0)
    {
      LOG_WARNING << "Reading eqdsk.rlim and eqdsk.zlim...";
      eqdsk.rlim.resize(eqdsk.nlim);
      eqdsk.zlim.resize(eqdsk.nlim);
      for (int i=0; i<eqdsk.nlim; i++)
      {
        eqdsk_file >> eqdsk.rlim[i] >> eqdsk.zlim[i];
      } 
      LOG_WARNING << "Number of eqdsk.rlim/eqdsk.zlim values read " << eqdsk.nlim;
    }
    else
    {
      LOG_FATAL << "Error reading limiter data from eqdsk";
    }

    // set equData attributes
    rmin = eqdsk.rgrid;
    zmin = eqdsk.zmid - eqdsk.zdim/2;
    rmax = eqdsk.rgrid+eqdsk.rdim;
    zmax = eqdsk.zmid + eqdsk.zdim/2;
    nr = nw-1; // why does smardda do this?
    nz = nh-1;
    dr = (rmax-rmin)/nr;
    dz = (zmax-zmin)/nz;
    psiqbdry = eqdsk.psibdry1;
    psiaxis = eqdsk.psimag1;
    psinorm = fabs(psiqbdry-psiaxis)/2;
    set_rsig();
    dpsi = fabs(rsig*(psiqbdry-psiaxis)/nw);

    LOG_WARNING << "DPSI =  " << dpsi;
    LOG_WARNING << "PSIQBDRY = " << psiqbdry;
    LOG_WARNING << "PSIAXIS = " << psiaxis;
    LOG_WARNING << "RSIG = " << rsig;
    LOG_WARNING << "RMIN = " << rmin;
    LOG_WARNING << "RMAX = " << rmax;
    LOG_WARNING << "ZMIN = " << zmin;
    LOG_WARNING << "ZMAX = " << zmax;
    LOG_WARNING << "dR = " << dr;
    LOG_WARNING << "dZ = " << dz;
    LOG_WARNING << "PSINORM = " << psinorm;



    // scale for psibig (TODO)
    // if psibig 
    // scale by 2pi
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
  eqdsk_out << eqdsk.header << std::endl;

  // Write out initial four lines of floats
  double parameters[20] = {eqdsk.rdim, eqdsk.zdim, eqdsk.rcentr, eqdsk.rgrid, eqdsk.zmid, eqdsk.rmaxis, eqdsk.zmaxis, eqdsk.psimag1, eqdsk.psibdry1, 
                eqdsk.bcentr, eqdsk.cpasma, eqdsk.psimag2, eqdsk.xdum, eqdsk.rmaxis, eqdsk.xdum, eqdsk.zmaxis, eqdsk.xdum, eqdsk.psibdry2, eqdsk.xdum, eqdsk.xdum};

  for (int i=0; i<20; i++)
  {
    element = parameters[i];
    counter = eqdsk_line_out(eqdsk_out, element, counter);
  }

  // Write out data
  eqdsk_write_array(eqdsk_out, eqdsk.fpol, counter); // write eqdsk.fpol array
  eqdsk_write_array(eqdsk_out, eqdsk.pres, counter); // write eqdsk.pres array 
  eqdsk_write_array(eqdsk_out, eqdsk.ffprime, counter); // write eqdsk.ffprime array
  eqdsk_write_array(eqdsk_out, eqdsk.pprime, counter); // write eqdsk.pprime array 

  // Write out psi(R,Z)
  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      element = eqdsk.psi[i][j];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
    }
  }
  if (counter < 5)
  {
    counter = 0;
    eqdsk_out << std::endl;
  }

  eqdsk_write_array(eqdsk_out, eqdsk.qpsi, counter); // write eqdsk.qpsi array
  eqdsk_out << std::endl << "  " << eqdsk.nbdry << "   " << eqdsk.nlim << std::endl; // write eqdsk.nbdry and eqdsk.nlim

  // write eqdsk.rbdry and eqdsk.zbdry
  for (int i=0; i<eqdsk.nbdry; i++)
  { 
      element = eqdsk.rbdry[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = eqdsk.zbdry[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = eqdsk.rbdry[i];
  }
  if (counter < 5)
  {
    counter = 0;
    eqdsk_out << std::endl;
  }

  // write eqdsk.rlim and eqdsk.zlim arrays
  for (int i=0; i<eqdsk.nlim; i++)
  { 
      element = eqdsk.rlim[i];
      counter = eqdsk_line_out(eqdsk_out, element, counter);
      element = eqdsk.zlim[i];
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

  double r_pts[nw];
  double z_pts[nh];
  r_pts[0] = rmin;
  z_pts[0] = zmin;
  
  
// loop over Z creating spline knots
  for (int i=0; i<nw; i++)
  {
    r_pts[i+1] = r_pts[i]+dr;
  }
// loop over Z creating spline knots
  for (int i=0; i<nh; i++)
  {
    z_pts[i+1] = z_pts[i]+dz;
  }

// loop over (R,Z) creating for spline knots
  double psi_pts[nw*nh];
  int count = 0;
  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      psi_pts[count] = eqdsk.psi[i][j];
      count += 1;
    }
  }

  double f_pts[nw];
  double psi_1dpts[nw];
  std::copy(eqdsk.fpol.begin(), eqdsk.fpol.end(), f_pts);


  // loop over R to create 1d knots of psi 
  psi_1dpts[0] = psiaxis;
  dpsi = rsig*dpsi; // set correct sign of dpsi depending on if increase/decrease outwards
  for (int i=0; i<nh; i++)
  {
    psi_1dpts[i+1] = psi_1dpts[i]+dpsi;
  }

  

  // set 1d arrays for R grid, Z grid and Psi grids
  r_grid.setcontent(nw, r_pts);
  z_grid.setcontent(nh, z_pts);
  psi_grid.setcontent(count, psi_pts);
  psi_1dgrid.setcontent(nw, psi_1dpts);
  f_grid.setcontent(nw, f_pts);


  // Construct the spline interpolant for flux function psi(R,Z) 
  alglib::spline2dbuildbilinearv(r_grid, nw, z_grid, nh, psi_grid, 1, psiSpline);
  
  // Construct the spline interpolant for toroidal component flux function I(psi) aka f(psi) 
  alglib::spline1dbuildlinear(psi_1dgrid ,f_grid, fSpline);

  // Alglib::spline2ddiff function can return a value of spline at (R,Z) as well as derivatives. 
  // I.e no need to have a separate spline for each derivative dPsidR and dPsidZ 



  // create 1d spline for f(psi) aka eqdsk.fpol
  //alglib::spline1dbuildlinear(f_grid,)

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
// (+1 -> Increase outwards, -1 -> Decrease outwards)
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

void equData::r_extrema()
{
  int nrsrsamp; // number of samples in r
  double zpsi; // psi
  double zdpdr; // dpsi/dr
  double zdpdz; // dpsi/dz
  double ztheta; // theta_j
  double zsrr; // estimate for maximum |R-R_c| in domain
  double zszz; // estimate for maximum |Z-Z_c| in domain
  double zsrmin; // estimate for starting r
  double zsrmax; // estimate for maximum r in domain
  double zcostheta; // cos(theta_j)
  double zsintheta; // sin(theta_j)
  double re; // R_i
  double ze; // Z_i
  double zdpdsr; // {dspi/dr}_{i}
  double zsr; // r_i
  double zdsr; // Delta r_i
  double zdpdsrl; // (dspi/dr)_{i-1}
  double zrpmin; // largest r giving psi < psi_{min}
  double zrpmax; // smallest r giving psi > psi_{max} 
  int isr; // flag that value psi < psi_{min} has been found
  int idsplet; // flag that (dpsi/dr)_{i-1} set
  int im; // number of angles in largest interval of acceptable theta
  int il; // marks lower bound of current interval of acceptable theta
  int ilm; // marks lower bound of largest interval of acceptable theta
  std::vector<double> work1(ntheta+1); // 1d work array

  std::fill(work1.begin(), work1.end(), 0); // fill working array with zeroes


  zsrr = std::max(fabs(rmax-rcen), fabs(rmin-rcen));
  zszz = std::max(fabs(zmax-zcen), fabs(zmin-zcen));
  zsrmax = sqrt(pow(zsrr,2) + pow(zszz, 2));
  nrsrsamp = nr+nz;
  zdsr = zsrmax/nrsrsamp; 
  zsrmin = zsrmax/10;
  dtheta = (thetamax-thetamin)/ntheta;

  // loop over angle
  for (int j=1; j<=ntheta+1; j++)
  {
    ztheta = thetamin + j*dtheta;
    zcostheta = cos(ztheta);
    zsintheta = sin(ztheta);

    // loop over distance from centre. Start a little away from origin
    zsr = zsrmin;
    isr = 0;
    idsplet = 0;

    for (int i=1; i<=nrsrsamp; i++)
    {
      re = rcen + zsr*zcostheta;
      if (re>=rmax || re<=rmin)
      {
        LOG_WARNING << "Rejecting direction-extremum not found. Rejecting in R";
        work1[j] = 2;
        break;
      }
      ze = zcen + zsr*zsintheta;
      if (ze>=zmax || ze<=zmin)
      {
        LOG_WARNING << "Rejecting direction-extremum not found. Rejecting in Z";
        work1[j] = 2;
        break;
      }
      zpsi = alglib::spline2dcalc(psiSpline, re, ze);
    }
  }

  if (rsig>0) // i.e if psiaxis < psiqbdry meaning psi increases outwards
  {
    //if (zpsi<psimin)
  }

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
  double ztheta; // theta_j
  double tmin; // minimum theta for R,Z(psi,theta)
  double zdpdr; // dpsi/dr
  double zdpdz; // dpsi/dz
  double zsr; // r
  double zcostheta; // cos(theta_j)
  double zsintheta; // sin(theta_j)
  double re; // R_i
  double ze; // Z_i
  double zdpdsr; // {dspi/dr}_{i}
  double zdsr; // Delta r_i
  double zdpdsrl; // {dspi/dr}_{i-1}
  double zpsi; // psi
  double zpsi_i; // psi_i
  

  // loop over angle

  int counter = 0;
  for (int j=0; j<ntheta; j++)
  { 
    counter += 1; 
    ztheta = thetamin + j*dtheta;
    zcostheta = cos(ztheta);
    zsintheta = sin(ztheta);


  }

}

// Caculate B field vector (in toroidal polars) at given position
std::vector<double> equData::b_field(std::vector<double> position,
                                     std::string startingFrom)
{
  std::vector<double> bVector(3); // B in toroidal polars
  double zr; // local R from position vector supplied
  double zz; // local Z from position vector supplied
  double zpsi; // local psi returned from spline calc
  double zdpdr; // local dpsi/dr from spline calc
  double zdpdz; // local dpsi/dz from spline calc
  double zf; // local f(psi) = RB_T (not sure what this means, copied comment from SMARDDA)
  double null; // second derivative of psi not needed

  // if position vector is already in polars skip coord transform and calculate B
  if (startingFrom == "polar")
  {
    zr = position[0];
    zz = position[1];
  }
  // otherwise transform from cartesian -> polar before calculating B
  else
  {
    std::vector<double> polarPosition;
    polarPosition = coordTfm::cart_to_polar(position, "forwards");
    zr = polarPosition[0];
    zz = polarPosition[1];
  }

  // evaluate psi and psi derivs at given (R,Z) coords
  alglib::spline2ddiff(psiSpline, zr, zz, zpsi, zdpdr, zdpdz, null); 

  // evaluate I aka f at psi (I aka f is the flux function) 
  zf = alglib::spline1dcalc(fSpline, zpsi);

  // calculate B in cylindrical toroidal polars
  bVector[0] = -zdpdz/zr; // B_R
  bVector[1] = zdpdr/zr; // B_Z
  bVector[2] = zf/zr; // B_T - toroidal component of field directed along phi

  // return the calculated B vector 
  return bVector;
}

std::vector<double> equData::b_field_cart(std::vector<double> polarBVector, double phi)
{
  std::vector<double> cartBVector(3);
  double zbx; // cartesian Bx component
  double zby; // cartesian By component
  double zbz; // Bz component
  double zbr; // polar Br component
  double zbt; // polar toroidal Bt component

  zbr = polarBVector[0];
  zbz = polarBVector[1];
  zbt = -polarBVector[2];

  zbx = zbr*cos(-phi) - zbt*sin(-phi);
  zby = zbr*sin(-phi) + zbt*cos(-phi);

  cartBVector[0] = zbx;
  cartBVector[1] = zby;
  cartBVector[2] = zbz; 

  return cartBVector;
}

void equData::write_bfield(bool plotRZ, bool plotXYZ)
{
  std::vector<double> polarPos(3); // cartesian position P(x,y,z)
  std::vector<double> polarB(3); // polar toroidal magnetic field B(Br, Bz, Bphi)
  
  if (plotRZ) // write out polar toroidal P(r,z) and B(Br,Bz)
  {
    std::ofstream BField_out_rz;
    BField_out_rz.open("BField_rz.txt");
   // BField_out_rz << "rPos" << " " << "zPos" << " " << "Br" << " " << "Bz" << std::endl;
   // BField_out_rz << std::setprecision(6) << std::fixed;

    for (int j=0; j<nh; j++)
    {
      for (int i=0; i<nw; i++)
      {
        polarPos[0] = r_grid[i];
        polarPos[1] = z_grid[j];
        polarB = b_field(polarPos, "polar");
        BField_out_rz << polarPos[0] << " " << polarPos[1] << " " 
                      << polarB[0] << " " << polarB[1] << std::endl;
      }
      BField_out_rz << std::endl;
    }
  }

  if (plotXYZ) // write out cartesian P(x,y,z) and cartesian B(Bx,By,Bz)
  {
    std::ofstream BField_out_xyz;
    BField_out_xyz.open("BField_xyz.txt");
   // BField_out_xyz << "xPos" << " " <<  "yPos" << " " << "zPos" << " " 
    //               << "Bx" << " " << "By" << " " << "Bz" << std::endl;
    BField_out_xyz << std::setprecision(6) << std::fixed;

    std::vector<double> cartPos(3); // polar toroidal position P(r,z,phi)
    std::vector<double> cartB(3); // cartesian magnetic field B(Bx, By, Bz)
    std::string tfmDir = "backwards"; // set transform direction to polar -> cart
    int phiSamples = 100;
    double dphi = 2*M_PI/phiSamples;

    for (int k=0; k<phiSamples; k++) // loop through phi
    {
      for (int j=0; j<nh; j++) // loop through z 
      {
        for (int i=0; i<nw; i++) // loop through r
        {
          polarPos[0] = r_grid[i]; // pull r position from eqdsk
          polarPos[1] = z_grid[j]; // pull z position from eqdsk
          // test for negative Rs
          if (polarPos[0] < 0 ) 
          {
            LOG_ERROR << "Negative R Values when attempting plot magnetic field. Fixup eqdsk required";
          }

          polarB = b_field(polarPos, "polar"); // calculate B(R,Z,phi)
          cartB = b_field_cart(polarB, polarPos[2]); // transform to B(x,y,z)
          cartPos = coordTfm::cart_to_polar(polarPos, tfmDir); // transform position to cartesian

          // write out magnetic field data for plotting
          if (cartPos[0] > 0 )
          {
            BField_out_xyz << " ";
          }
          BField_out_xyz << cartPos[0] << " " << cartPos[1] << " " << cartPos[2] << " "
                         << cartB[0] << " " << cartB[1] << " " << cartB[2] << " " 
                         << polarPos[0] << " " << polarPos[1] << " " << polarPos[2] << std::endl;
        }
      }
      polarPos[2] += dphi;
      BField_out_xyz << std::endl;
    }
  }
}
