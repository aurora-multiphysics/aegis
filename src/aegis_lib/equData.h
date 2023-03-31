#ifndef equData__
#define equData__

#include <stdio.h>
#include <vector>
#include <fstream>
#include "alglib/interpolation.h"



class equData{
  private:
  double rsig; // sign of dpsi/dr value

  public:

  // File streams
  std::ifstream eqdsk_file;
  std::ofstream eqdsk_out;

  // Header string
  std::string header; // string containing header information

  // Integer data
  int nw; // Number of horizontal R points in grid
  int nh; // Number of vertical Z points in grid
  int nbdry; // Number of boundary points
  int nlim; // Number of limiter points

  int nr; // nw-1 (used for finite difference)
  int nz; // nh-1 (used for finite difference)


  // eqdsk Floating point data
  double rdim; // EFIT Horizontal dimension in metre of computational box
  double zdim; // EFIT Vertical dimension in metre of computational box
  double rcentr; // EFIT Radial dimension at centre of plasma
  double rgrid; // EFIT Minimum R in metre of rectangular computational box
  double zmid; // EFIT centre of computational box in metres
  double rmaxis; // EFIT plasma centre (magnetic axis) in metre (R at central psi extremum)
  double zmaxis; // EFIT plasma centre (magnetic axis) in metre (Z at central psi extremum)
  double psimag1; // EFIT poloidal flux at the magnetic axis (Wb/rad)
  double psibdry1; // EFIT poloidal flux at the boundary (Wb/rad)
  double bcentr; // EFIT vacuum toroidal field at r=rcentr IGNORED BY SMARDDA
  double cpasma; // EFIT computed plasma current in A IGNORED BY SMARDDA
  double psimag2; // EFIT poloidal flux at the magnetic axis (Wb/rad) IGNORED BY SMARDDA
  double psibdry2; // EFIT poloidal flux at the boundary (Wb/rad) IGNORED BY SMARDDA
  double xdum; // empty dummy variables in data

  // attributes from SMARDDA
  double psiaxis;
  double psiqbdry;
  double rcen; // rcen calculated from centre()
  double zcen; // zcen calculated from centre()


  double rmin; // min R value in equillibrium data
  double zmin; // min Z value in equillibrium data
  double rmax; // max R value in equillibrium data
  double zmax; // max Z value in equillibrium data

  double dr; // step size in R (rmax-rmin)
  double dz; // step size in Z (zmax-zmin)


  // 1D array data
  std::vector<double> fpol; // fpol size(nw)
  std::vector<double> pres; // pres size(nw)
  std::vector<double> ffprime; // ffprime(nw)
  std::vector<double> pprime; // pprime(nw)

  alglib::real_1d_array r_grid, z_grid, psi_grid;

  std::vector<double> qpsi; // qpsi(nw)
  std::vector<double> rbdry; // Boundary points in R
  std::vector<double> zbdry; // Boundary points in Z
  std::vector<double> rlim; // Limiter points in R  
  std::vector<double> zlim; // Limiter points in Z  


  // 2D array data
  std::vector<std::vector<double>> psi; // psi(R,Z) size(nw*nh)

  // Splines
  alglib::spline2dinterpolant psiSpline; // 2d spline interpolant for Psi(R,Z)

  // Methods
  void read_eqdsk(std::string filename);
  std::vector<double> read_array(int n, std::string varName);
  std::vector<std::vector<double>> read_2darray(int nx, int ny, std::string varName);
  void write_eqdsk_out();
  int eqdsk_line_out(std::ofstream &file, double element, int counter);
  void eqdsk_write_array(std::ofstream &file, std::vector<double> array, int counter);
  void centre();
  void set_rsig();
  void init_interp_splines();
  void gnuplot_out();
};

#endif

                                          

