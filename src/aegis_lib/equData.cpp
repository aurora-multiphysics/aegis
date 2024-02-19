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
#include "coordtfm.h"


void equData::setup(const std::shared_ptr<InputJSON> &inputs){

  json equilNamelist;
  if (inputs->data.contains("equil_params"))
  { 
    equilNamelist = inputs->data["equil_params"];
    eqdskFilepath = equilNamelist["eqdsk"];
    cenopt = equilNamelist["cenopt"];
    powerSOL = equilNamelist["P_sol"];
    lambdaQ = equilNamelist["lambda_q"];
    psiref = equilNamelist["psiref"];
    rOutrBdry = equilNamelist["r_outrbdry"];
    rmove = equilNamelist["rmove"];
    drawEquRZ = equilNamelist["draw_equil_rz"];
    drawEquXYZ = equilNamelist["draw_equil_xyz"];
    debug = equilNamelist["print_debug_info"];
  }

  std::cout << "eqdskFilePath = " << eqdskFilepath << std::endl;

  read_eqdsk(eqdskFilepath);
}
void equData::psiref_override()
{
  psibdry = psiref;
}

// Return eqdsk data structure
eqdskData equData::get_eqdsk_struct()
{
  return eqdsk;
}

// Read eqdsk file
void equData::read_eqdsk(std::string filename)
{
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "--------------------READING EQDSK---------------------" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    LOG_TRACE << "-----equData.read_eqdsk()-----";
    eqdsk_file.open(filename);
    std::stringstream header_ss;
    std::string temp;
    std::vector<int> header_ints;

    LOG_WARNING << "eqdsk file to be read - " << filename;
    if (!eqdsk_file.is_open())
    {
      std::cout << "Error! Could not open eqdsk file " << filename << std::endl;
    }
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


    // Read first four lines of data
    eqdsk_file >> eqdsk.rdim >> eqdsk.zdim >> eqdsk.rcentr >> eqdsk.rgrid >> eqdsk.zmid;
    eqdsk_file >> eqdsk.rqcen >> eqdsk.zqcen >> eqdsk.psimag1 >> eqdsk.psibdry1 >> eqdsk.bcentr;
    eqdsk_file >> eqdsk.cpasma >> eqdsk.psimag2 >> eqdsk.xdum >> eqdsk.rqcen >> eqdsk.xdum;
    eqdsk_file >> eqdsk.zqcen >> eqdsk.xdum >> eqdsk.psibdry2 >> eqdsk.xdum >> eqdsk.xdum;

    if (debug)
    { 
      LOG_WARNING << "Number of grid points in R (nw) = " << nw;
      LOG_WARNING << "Number of grid points in Z (nh) = " << nh;
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
      LOG_WARNING << "Plasma centre in R eqdsk.rmaxis " << eqdsk.rqcen;
      LOG_WARNING << "Plasma centre in Z eqdsk.zmaxis " << eqdsk.zqcen;
      LOG_WARNING << "Plasma current " << eqdsk.cpasma;
    }
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
      LOG_INFO << "Reading eqdsk.rbdry and eqdsk.zbdry...";
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
      LOG_INFO << "Reading eqdsk.rlim and eqdsk.zlim...";
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

    //initialise rest of equData attributes
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
   
    ivac = eqdsk.fpol[0];
    psifac = 1.0;

    // Fix for ITER eqdsks

    OVERRIDE_ITER = true;
    
    if (OVERRIDE_ITER == true)
    {
     LOG_INFO << "Override for ITER case applied. Sign of psi flipped";
     if (eqdsk.psimag1 > eqdsk.psibdry1)
      {
        psiaxis = -eqdsk.psimag1*psifac; // flip sign
        psiqbdry = -eqdsk.psibdry1*psifac; // flip sign
        psibdry = psiqbdry;
        for (auto& row:eqdsk.psi)
        {
          for (auto& i:row) 
            {
              i = -i;
              //std::cout << i << std::endl;
            }
        }
      }
      else 
      {
        for (auto i:eqdsk.fpol) {i = -i;}
      }
    } 
    
    set_rsig();
    dpsi = std::abs(rsig*(psiqbdry-psiaxis)/nw);
    psinorm = std::fabs(psiqbdry-psiaxis)/2;

    if (debug)
    { 
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
      LOG_WARNING << "IVAC = " << ivac;
    }


    // scale for psibig (TODO)
    // if psibig
    // scale by 2pi


}

// Read 1D array from eqdsk (PRIVATE)
std::vector<double> equData::read_array(int n, std::string varName)
{
  LOG_TRACE << "-----equData.read_array()-----";

  std::vector<double> work1d(n); // working vector of doubles of size n
  LOG_INFO << "Reading " << varName << " values...";
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
  LOG_TRACE << "-----equData.read_2darray()-----";

  std::vector<std::vector<double>> work2d(nw, std::vector<double>(nh));
  LOG_INFO << "Reading " << varName << " values...";
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
  LOG_TRACE << "-----equData.write_eqdsk_out()-----";
  std::ofstream eqdsk_out("eqdsk.out");
  eqdsk_out << std::setprecision(9) << std::scientific;
  eqdsk_out << std::setiosflags(std::ios::uppercase);
  double element;
  int counter=0;

  // Write out header information
  eqdsk_out << eqdsk.header << std::endl;

  // Write out initial four lines of floats
  double parameters[20] = {eqdsk.rdim, eqdsk.zdim, eqdsk.rcentr, eqdsk.rgrid, eqdsk.zmid, eqdsk.rqcen, eqdsk.zqcen, eqdsk.psimag1, eqdsk.psibdry1,
                eqdsk.bcentr, eqdsk.cpasma, eqdsk.psimag2, eqdsk.xdum, eqdsk.rqcen, eqdsk.xdum, eqdsk.zqcen, eqdsk.xdum, eqdsk.psibdry2, eqdsk.xdum, eqdsk.xdum};

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
      LOG_TRACE << "-----equData.eqdsk_line_out()-----";
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
  LOG_TRACE << "-----equData.eqdsk_write_array()-----";
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
  LOG_TRACE << "-----equData.init_interp_splines()-----";
  std::vector<double> r_pts(nw);
  std::vector<double> z_pts(nh);
  r_pts[0] = rmin;
  z_pts[0] = zmin;

  //set_rsig();


// loop over Z creating spline knots
  for (int i=1; i<nw; i++)
  {
    r_pts[i] = r_pts[i-1]+dr;
  }
// loop over Z creating spline knots
  for (int i=1; i<nh; i++)
  {
    z_pts[i] = z_pts[i-1]+dz;
  }

// loop over (R,Z) creating for spline knots
  std::vector<double> psi_pts(nw*nh);
  int count = 0;
  for (int i=0; i<nw; i++) // flatten Psi(R,Z) array into Psi(index)
  {
    for(int j=0; j<nh; j++)
    {
      psi_pts[count] = eqdsk.psi[i][j];
      count += 1;
    }
  }

  std::vector<double> psi_1dpts;
  psi_1dpts.reserve(nw);

  // loop over R to create 1d knots of psi
  psi_1dpts.push_back(psiaxis);

  dpsi = rsig*dpsi; // set correct sign of dpsi depending on if increase/decrease outwards
  for (int i=1; i<nw; i++)
  {
    psi_1dpts.push_back(psi_1dpts[i-1]+dpsi);
  }

  // set 1d arrays for R grid, Z grid and Psi grids
  r_grid.setcontent(nw, r_pts.data());
  z_grid.setcontent(nh, z_pts.data());
  psi_grid.setcontent(count, psi_pts.data());
  psi_1dgrid.setcontent(nw, psi_1dpts.data());
  f_grid.setcontent(nw, eqdsk.fpol.data());

//  for (int i=0 ; i<nh; i++)
//  {
//    psi1d_out << psi_1dgrid[i] << std::endl;
//  }


  // Construct the spline interpolant for flux function psi(R,Z)
  alglib::spline2dbuildbicubicv(r_grid, nw, z_grid, nh, psi_grid, 1, psiSpline);

  // Construct the spline interpolant for toroidal component flux function I(psi) aka f(psi)
  alglib::spline1dbuildcubic(psi_1dgrid ,f_grid, fSpline);

  // Alglib::spline2ddiff function can return a value of spline at (R,Z) as well as derivatives.
  // I.e no need to have a separate spline for each derivative dPsidR and dPsidZ

  double R,Z,PSI;
  R = 6;
  Z = 0.41;
  PSI = 0;
  PSI = alglib::spline2dcalc(psiSpline, R, Z);

  // create 1d spline for f(psi) aka eqdsk.fpol
  //alglib::spline1dbuildlinear(f_grid,)

}

// Write out psi(R,Z) data for gnuplotting
void equData::gnuplot_out()
{
  LOG_TRACE << "-----equData.gnuplot_out()-----";
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
  LOG_TRACE << "-----equData.set_rsig()-----";
  if (psiqbdry-psiaxis < 0)
  {
    rsig = -1.0;
  }
  else if (psiqbdry-psiaxis > 0)
  {
    rsig = +1.0;
  }
  
  if (debug) {
    LOG_WARNING << "Value of rsig (psiqbdry-psiaxis) = " << rsig;
  }
}

// Find central psi extrema
void equData::centre(int cenopt)
{
  LOG_TRACE << "-----equData.centre-----";
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


  switch (cenopt)
  {
  case 1:
    rcen = eqdsk.rqcen;
    zcen = eqdsk.zqcen;
    LOG_INFO << "(Rcen,Zcen) values taken from eqdsk (Rmaxis,Zmaxis)";
    break; // end case 1

  case 2:
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
    

    LOG_INFO << "New (Rcen,Zcen) values calculated";
    break; // end case 2
  }
  psicen = alglib::spline2dcalc(psiSpline, rcen, zcen);
  LOG_WARNING << "Rcen value calculated from equData.centre() = " << rcen;
  LOG_WARNING << "Zcen value calculated from equData.centre() = " << zcen;
  LOG_WARNING << "Psicen value calculated from equData.centre() = " << psicen;
  

}

// ***TODO***
void equData::r_extrema()
{
  LOG_TRACE << "-----equData.r_extrema()-----";
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
  int idplset; // flag that (dpsi/dr)_{i-1} set
  int im; // number of angles in largest interval of acceptable theta
  int il; // marks lower bound of current interval of acceptable theta
  int ilm; // marks lower bound of largest interval of acceptable theta
  std::vector<double> work1(ntheta+1); // 1d work array

  std::fill(work1.begin(), work1.end(), 0); // fill working array with zeroes


  zsrr = std::max(std::abs(rmax-rcen), std::abs(rmin-rcen));
  zszz = std::max(std::abs(zmax-zcen), std::abs(zmin-zcen));
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
    idplset = 0;

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


// Create 2d spline structures for R(psi,theta) and Z(psi,theta) ***TODO***
void equData::rz_splines()
{
  LOG_TRACE << "-----equData.rz_splines()-----";
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

  double ntmax = 33;
  double ntmin = 1;

  intheta = ntmax - ntmin+1;

  // loop over angle

  int counter = 0;
  for (int j=0; j<ntheta; j++)
  {
    counter += 1;
    ztheta = thetamin + j*dtheta;
    zcostheta = cos(ztheta);
    zsintheta = sin(ztheta);
    std::cout << "THETAMIN = " << thetamin << std::endl;


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
  double zf; // local f(psi) = RB_T toroidal component of B
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
  
  if (zr < rmin || zr > rmax || zz < zmin || zz > zmax )
  {
    return std::vector<double>();
  }
  else
  {
    // evaluate psi and psi derivs at given (R,Z) coords
    alglib::spline2ddiff(psiSpline, zr, zz, zpsi, zdpdr, zdpdz, null);

    // evaluate I aka f at psi (I aka f is the flux function)
    zf = alglib::spline1dcalc(fSpline, zpsi); 

    // calculate B in cylindrical toroidal polars
    bVector[0] = -zdpdz/zr; // B_R
    bVector[1] = zdpdr/zr; // B_Z
    bVector[2] = zf/zr; // B_T - toroidal component of field directed along phi
  }

  // return the calculated B vector
  return bVector;

}

std::vector<double> equData::b_field_cart(std::vector<double> polarBVector, double phi, int normalise)
{
  std::vector<double> cartBVector(3);
  double zbx; // cartesian Bx component
  double zby; // cartesian By component
  double zbz; // Bz component
  double zbr; // polar Br component
  double zbt; // polar toroidal Bt component

  zbr = polarBVector[0];
  zbz = polarBVector[1];
  zbt = polarBVector[2];

  zbx = zbr*cos(phi) - zbt*sin(phi);
  zby = -(zbr*sin(phi) + zbt*cos(phi));

  cartBVector[0] = zbx;
  cartBVector[1] = zby;
  cartBVector[2] = zbz;

  if (normalise == 1)
  {
    double norm;
    norm = sqrt(pow(cartBVector[0],2) + pow(cartBVector[1],2) + pow(cartBVector[2],2));
    cartBVector[0] = cartBVector[0]/norm;
    cartBVector[1] = cartBVector[1]/norm;
    cartBVector[2] = cartBVector[2]/norm;
  }

  return cartBVector;
}



void equData::write_bfield(bool plotRZ, bool plotXYZ)
{
  LOG_TRACE << "-----equData.b_field_write()-----";
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
    std::ofstream aegisB;
    aegisB.open("aegis_b.txt");
    BField_out_xyz.open("BField_xyz.txt");
   // BField_out_xyz << "xPos" << " " <<  "yPos" << " " << "zPos" << " "
    //               << "Bx" << " " << "By" << " " << "Bz" << std::endl;
    BField_out_xyz << std::setprecision(8) << std::fixed;
    aegisB << std::setprecision(8) << std::fixed;

    std::vector<double> cartPos(3); // polar toroidal position P(r,z,phi)
    std::vector<double> cartB(3); // cartesian magnetic field B(Bx, By, Bz)
    int phiSamples = 12;
    double dphi = 2*M_PI/phiSamples;

    BField_out_xyz << "x" << " " << "y" << " " << "z" << " "
		   << "Bx" << " " << "By" << " " << "Bz" << " "
		   << "r" << " " << "z" << " " << "phi" << std::endl;

    for (int k=0; k<=phiSamples; k++) // loop through phi
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
          cartB = b_field_cart(polarB, polarPos[2], 0); // transform to B(x,y,z)
          cartPos = coordTfm::cart_to_polar(polarPos, "backwards"); // transform position to cartesian

          // write out magnetic field data for plotting
          if (cartPos[0] > 0 )
          {
            BField_out_xyz << " ";
          }
          BField_out_xyz << cartPos[0] << " " << cartPos[1] << " " << cartPos[2] << " "
                         << cartB[0] << " " << cartB[1] << " " << cartB[2] << " "
                         << polarPos[0] << " " << polarPos[1] << " " << polarPos[2] << std::endl;

          aegisB << polarB[0] << " " << polarB[1] << " " << polarB[2] << std::endl;
        //polarPos[2] += dphi; // Incorrect place but much faster delaunay (leads to inconsistencies)

        }
      }
      polarPos[2] += dphi; // Correct place for this to be. But significantly slower delaunay 3d in paraview
    //  BField_out_xyz << std::endl;
    }
  }

}

// get toroidal ripple term in magnetic field from external file or otherwise. By default use MAST term
// std::vector<double> equData::b_ripple(std::vector<double> pos, std::vector<double> bField)
// {
//   LOG_TRACE << "-----equData.b_ripple()-----";
//   std::vector<double> bRipple(3); // toroidal ripple term
//   double zr; // local R
//   double zz; // local Z
//   double zphi; // local phi
//   double zh; // local H
//   double zdhdr; // frac{partial H}{partial R}
//   double zdhdz; // frac{partial H}{partial Z}

//   // eventually can add in the ability to read in files but for now this will just define the MAST ripple





// }




void equData::boundary_rb()
{
  // mostly copied from r_extrema()

  LOG_TRACE << "-----equData.boundary_rb()-----";
  std::string functionName = ".boundary_rb(): ";
  int nrsrsamp; // number of samples in r

  double zddpdz; // second derivative of psi not needed from spline2ddiff

  double zpsi; // psi
  double zpsiinr; // estimate for inner limit of psi
  double zpsioutr; // estimate for outer limit of psi
  double zdpdr; // dpsi/dr
  double zdpdz; // dpsi/dz
  double ztheta; // theta_j
  double zsrr; // estimate for maximum |R-R_c| in domain
  double zszz; // estimate for maximum |Z-Z_c| in domain
  double zsrmin; // estimate for starting r
  double zsrmax; // estimate for maximum r in domain
  double zsrinr; // smallest r in range
  double zsroutr; // largest r in range
  double zsrsta; // estimate for starting r
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
  int idplset; // flag that (dpsi/dr)_{i-1} set
  int im; // number of angles in largest interval of acceptable theta
  int il; // marks lower bound of current interval of acceptable theta
  int ilm; // marks lower bound of largest interval of acceptable theta



  //ztheta = M_PI; // either pick inboard
  ztheta = 0.0; // or outboard

  zpsiinr = (psiqbdry+psiaxis)/2;
  zsrr = std::max(std::abs(rmax-rcen), std::abs(rmin-rcen));
  zszz = std::max(std::abs(zmax-zcen), std::abs(zmin-zcen));
  zsrmax = sqrt(pow(zsrr,2) + pow(zszz, 2));
  nrsrsamp = nr+nz;
  zdsr = zsrmax/nrsrsamp;
  zsrsta = zsrmax/10;
  zpsioutr = zpsiinr;


  zcostheta = cos(ztheta);
  zsintheta = sin(ztheta);

  // loop over distance from centre. Start a little away from origin
  zsr = zsrsta;
  zsrinr = zsr;
  isr = 0;
  idplset = 0;

  for (int i=1; i<=nrsrsamp; i++)
  {
    re = rcen + zsr*zcostheta;
    if (re>=rmax || re<=rmin)
    {
      // extremal R reached
      break;
    }
    ze = zcen + zsr*zsintheta;
    if (ze>=zmax || ze<=zmin)
    {
      // extremal Z reached
      break;
    }
    zpsi = alglib::spline2dcalc(psiSpline, re, ze);


    if (rsig>0) // psiaxis < psiqbdry (psi increasing outwards)
    {
      if ((zpsi-zpsiinr)*rsig<0)
      {
        // record while psi < psi-inner
        isr = i;
        zsrinr = zsr;
      }
      else
      {
        if (isr==0)
        {
          // sr sta is too large, reduce
          zsr = 4*zsr/5;
          continue;
        }
      }
    }
    else // psiaxis > psiqbdry (psi decreasing outwards)
    {
      if ((zpsi-zpsiinr)*rsig<0)
      {
        // record while psi > psi-inner
        isr = i;
        zsrinr = zsr;
      }
      else
      {
        if (isr==0)
        {
          // sr sta is too large, reduce
          zsr = 4*zsr/5;
          continue;
        }
      }
    }

    // check monotone
    alglib::spline2ddiff(psiSpline, re, ze, zpsi, zdpdr, zdpdz, zddpdz);
    zdpdsr = zdpdr*zcostheta + zdpdz*zsintheta;
    if (idplset==1)
    {
      if (zdpdsrl*zdpdsr<=0)
      {
        // extremum in psi reached
        // try to use monotone range
        break;
      }
    }
    zsr = zsr+zdsr;
    zdpdsrl = zdpdsr;
    idplset = 1;
    zsroutr = zsr;
    zpsioutr = zpsi;
  }


  if ((zpsiinr-zpsioutr)*rsig >=0)
  {
    LOG_FATAL << className << functionName << "No suitable psi range exists";
  }

  zsrmin = zsrinr;
  zsrmax = zsroutr;

  // end of r_extrema() equivalent code


  int iext=4*npsi; // maximum allowed number of knots for spline in r
  int iknot; // actual number of knots for spline in r
  double cpsi; // constant for estimating delta r_i
  double zdsrmin; // floor to delta r_i
  std::vector <double> wvextn(iext); // 1D work array extended size for nodes
  std::vector <double> wvext(iext); // 1D work array extended size for samples
  std::vector <double> wvextd(iext); // 1D work array extended size for derivatives


  cpsi = std::abs(zpsioutr-zpsiinr)/(2*npsi);
  zdsrmin = (zsrmax-zsrmin)/(4*npsi);

  zsr = zsrmin;

  // loop over r to define 1-D splines as a function of r
  for (int i=1; i<=iext; i++)
  {
    re = rcen + zsr*zcostheta;
    ze = zcen + zsr*zsintheta;
    alglib::spline2ddiff(psiSpline, re, ze, zpsi, zdpdr, zdpdz, zddpdz);
    zdpdsr = zdpdr*zcostheta + zdpdz*zsintheta;

    if (i>1 && zdpdsrl*zdpdsr<=0)
    {
      // should only be small non-monotone region, try to ignore
      LOG_ERROR << className << functionName << "Lack of monotincity";
      // fix up
      zdpdsr = zdpdsrl;
    }
    wvextn[i] = zsr;
    wvext[i] = zpsi;
    wvextd[i] = zdpdsr;
    iknot = i;
    zdsr = std::abs(cpsi/zdpdsr);
    zsr = zsr + std::max(std::min(zdsr, 3*zdsrmin), zdsrmin);
    if (i>1 && (zsr-zsrmax>0))
    {
      break;
    }
    zdpdsrl = zdpdsr;
  }

  // LOG_FATAL << className << functionName << "Too many nodes";

  int intv; // interval in which spline inverse found

  for (int j=1; j<=iknot-1; j++)
  {
    intv = j;
    if ((wvext[j]-psiqbdry)*(wvext[j+1]-psiqbdry) <= 0) {break;}
  }
  zsr = wvextn[intv];


  zpsi = psiqbdry;
  LOG_WARNING << "psiqbdry from EQDSK = " << zpsi;
  re = rcen + zsr*zcostheta;
  ze = zcen + zsr*zsintheta;
  alglib::spline2ddiff(psiSpline, re, ze, zpsi, zdpdr, zdpdz, zddpdz);
  psibdry = zpsi;

  LOG_WARNING << "psibdry calculated from splines = " << zpsi;



  // store R_m and B_pm

  rbdry = re;
  zbdry = ze;
  bpbdry = (1/re)*sqrt(std::max(0.0, (pow(zdpdr,2) + pow(zdpdz, 2))));

  double zbr; // radial Bfield component
  double zbz; // vertical Bfield component
  double zbt; // toroidal Bfield component
  double zf; // toroidal component of B
  double cylj; // cylindrical current

  zf = alglib::spline1dcalc(fSpline, zpsi);
  btotbdry = sqrt(std::max(0.0, pow(bpbdry, 2) + pow((zf/re), 2)));

  LOG_WARNING << className << functionName << "Reference Boundary values";
  LOG_WARNING << "psiqbdry (psi_m) " << psiqbdry;
  LOG_WARNING << "rbdry (R_m) " << rbdry;
  LOG_WARNING << "zbdry (Z_m) " << ze;
  LOG_WARNING << "bpbdry (B_pm) " << bpbdry;
  LOG_WARNING << "btotbdry (B_m) " << btotbdry;

  cylj = std::abs(zsr*bpbdry/2.0e-7);
  std::cout << "ZSR = " << zsr << " BPBDRY = " << bpbdry << std::endl;
  LOG_WARNING << "Estimated cylindrical current " << cylj;

  zbr = -(1/re)*zdpdz;
  zbz = (1/re)*zdpdr;
  zbt = zf/re;

  LOG_WARNING << "Radial BField component (Br) " << zbr;
  LOG_WARNING << "Vertical BField component (Bz) " << zbz;
  LOG_WARNING << "Toroidal BField component (Bt) " << zbt;
}

double equData::omp_power_dep(double psi, double bn, std::string formula)
{
  double heatFlux;
  double exponential, fpfac;
  double zpsim = psibdry;
  double zbpm = bpbdry;
  double zrm = rbdry;
  double zrbfac; // power normalisation factor
  double zrblfac; // power factor scaled by decay length
  double rblfac;
  double power_split;
  double psign = -1.0;

  //zrm = 4.0422466206709089;
  //zbpm = 0.98301390593792659;

  zrbfac = 1/(2*M_PI*zrm*zbpm);


  zrblfac = zrbfac/lambdaQ;
  power_split = 0.5;

  // std::cout << "RBDRY = " << rbdry << std::endl;
  // std::cout << "BPBDRY = " << bpbdry << std::endl;
  // std::cout << "PSIBDRY = " << psibdry << std::endl;






  rblfac = 2*M_PI*zrblfac*(-1.0)*psign;



  if (formula == "exp")
  {
    fpfac = power_split*powerSOL*zrblfac;
    heatFlux = fpfac*bn*exp(rblfac*psi); // added bn into this formula because it was missing
  }
  //std::cout << "FPFAC = " << fpfac << std::endl;
  //std::cout << "RBLFAC = " << rblfac << std::endl;
    //std::cout << psi << std::endl;
  return std::fabs(heatFlux);
}


void equData::psi_limiter(std::vector<std::vector<double>> vertices)
{
  std::ofstream fluxVert("flux_vertex.txt");
  std::ofstream cartVert("cart_vertex.txt");

  double zpsimin=1.E+8;
  double zpsimax=-1.E+8;
  double zthetamin=1.E+8;
  double zthetamax=-1.E+8;
  double zrmin, zrmax;
  double zzmin, zzmax;
  double zr, zz, zpsi;
  double null;
  for (auto &i:vertices)
  {
    std::vector<double> cartPos = i;
    std::vector<double> polarPos = coordTfm::cart_to_polar(i, "forwards");
    std::vector<double> fluxPos = coordTfm::polar_to_flux(polarPos, "forwards", *this);
    zpsi = fluxPos[0];
    double ztheta = fluxPos[1];
    zr = polarPos[0];
    zz = polarPos[1];

    cartVert << cartPos[0] << " " << cartPos[1] << " " << cartPos[2] << std::endl;
    fluxVert << fluxPos[0] << " " << fluxPos[1] << " " << fluxPos[2] << std::endl;

    if (zpsi < zpsimin)
    {
      zpsimin = zpsi;
      zrmin = zr;
      zzmin = zz;
    }
    if (zpsi > zpsimax)
    {
      zpsimax = zpsi;
      zrmax = zr;
      zzmax = zz;
    }

    zthetamin = std::min(zthetamin, ztheta);
    zthetamax = std::max(zthetamax, ztheta);




  }

  double zdpdr;
  double zdpdz;
  alglib::spline2ddiff(psiSpline, zr, zz, zpsi, zdpdr, zdpdz, null);

  double zbpbdry, zbtotbdry;
  zbpbdry=(1/zr)*sqrt( std::max(0.0,(pow(zdpdr,2) + pow(zdpdz,2))));
  double zf = alglib::spline1dcalc(fSpline, zpsi);
  zbtotbdry = sqrt(std::max(0.0, (pow(zdpdr, 2) + pow(zdpdz, 2))));

  rbdry = zr;
  bpbdry = zbpbdry;
  btotbdry = zbtotbdry;

  std::cout << "PSI_LIMITER() RBDRY = " << rbdry << std::endl;
  std::cout << "PSI_LIMITER() BPBDRY = " << bpbdry << std::endl;
  std::cout << "PSI_LIMITER() BTOTBDRY = " << btotbdry << std::endl;


}

void equData::move()
{
  double zgfac; // geometrical factor
  std::cout << "BField Data has been moved and scaled by:" << std::endl;
  std::cout << "RMOVE = " << rmove << std::endl;
  std::cout << "ZMOVE = " << zmove << std::endl;
  std::cout << "FSCALE = " << fscale << std::endl;

  // move (R,Z) values
  rmin = rmin + rmove;
  rmax = rmax + rmove;
  zmin = zmin + zmove;
  zmax = zmax + zmove;
  eqdsk.rqcen = eqdsk.rqcen + rmove;
  eqdsk.zqcen = eqdsk.zqcen + zmove;

  zgfac = eqdsk.rqcen/(eqdsk.rqcen-rmove);
  std::cout << "MOVE FUNCTION INCOMPLETE. NEED TO FIX" << std::endl;
  // adjust \f$ f \f$ approximately (should use R instead of R_c,
  // but then f becomes a function of \f$ \theta \f$)

  for (int i=0; i<nw; i++) {i*zgfac;}
  // scale f
  for (int i=0; i<nw; i++) {i*fscale;}

}

std::array<double, 3> equData::get_midplane_params()
{
  std::array<double, 3> midplaneParams;
  midplaneParams[0] = rbdry;
  midplaneParams[1] = rOutrBdry;
  midplaneParams[2] = zcen;

  return midplaneParams;
}


// void equData::set_omp()
// {
//   // Ensure equData::centre() has been called first so R,Z values of centre known
//   if (rcen == 0 && zcen == 0)
//   {
//     LOG_ERROR << "Values RCEN and ZCEN are not set. Ensure these are set before attempting to
//                   define outer-midplane";
//   }

//   else
//   {
//     std::vector
//   }

// }
