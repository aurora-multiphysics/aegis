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



// Read eqdsk


void equData::read_eqdsk(std::string filename)
{
    std::ifstream eqdsk_file(filename);
    std::string line;
    std::stringstream line_ss;
    std::string temp;
    std::vector<int> header_ints;

    LOG_WARNING << "eqdsk file to be read - " << filename;

    // Extract header information
    std::getline(eqdsk_file, line);
    line_ss << line;
    int number_found;

    while (!line_ss.eof()) // Search for solitary numbers in header information
    {
      line_ss >> temp;
      if (std::stringstream(temp) >> number_found)
      {
        header_ints.push_back(number_found);
      }
    }
    // Store nw and nh from header information (last two numbers in header)
    nw = header_ints[header_ints.size()-2];
    nh = header_ints[header_ints.size()-1];
    LOG_WARNING << "nw found from eqdsk file to be = " << nw;
    LOG_WARNING << "read found from eqdsk file to be = " << nh;

    // Read first four lines of data
    eqdsk_file >> xdim >> zdim >> rcentr >> rgrid >> zmid;
    eqdsk_file >> rmaxis >> zmaxis >> psimag1 >> psibdry1 >> bcentr; 
    eqdsk_file >> cpasma >> psimag2 >> xdum >> rmaxis >> xdum; 
    eqdsk_file >> zmaxis >> xdum >> psibdry2 >> xdum >> xdum; 

    LOG_WARNING << "Geometry parameters from EFIT:";
    LOG_WARNING << "Domain size in R xdim " << xdim;
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

    if (nw>0)
    {
      LOG_WARNING << "Reading fpol values...";
      fpol = read_array(eqdsk_file);
      LOG_WARNING << "Reading pres values";
      pres = read_array(eqdsk_file);
      LOG_WARNING << "Reading ffprime values";
      ffprime = read_array(eqdsk_file);
      LOG_WARNING << "Reading pprime values";
      pprime = read_array(eqdsk_file);      
    }
    else
    {
      LOG_FATAL << "Unable to read 1D data from eqdsk";
    }

    // Read Psi(R,Z) data
    
    LOG_WARNING << "Reading Psi(R,Z) values:";
    psi = read_2darray(eqdsk_file);
    
}



// Read 1D array from eqdsk
std::vector<double> equData::read_array(std::ifstream &eqdsk_file)
{
  std::vector<double> work1d(nw); // working vector of doubles of size n 
  for(int i=0; i<nw; i++) // Read in n elements into vector from file
  {
    eqdsk_file >> work1d[i];
    LOG_WARNING << work1d[i];
  }
  return work1d;
}

// Read 2D array from eqdsk
std::vector<std::vector<double>> equData::read_2darray(std::ifstream &eqdsk_file)
{
  std::vector<std::vector<double>> work2d(nw, std::vector<double>(nh));

  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      eqdsk_file >> work2d[i][j];
    }
  }
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

  // Write out initial four lines of floats
  eqdsk_out << " " << xdim << " " << zdim << " " << rcentr << " " << rgrid << " " << zmid << std::endl;
  eqdsk_out << " " << rmaxis << " " << zmaxis << " " << psimag1 << " " << psibdry1 << " " << bcentr << std::endl; 
  eqdsk_out << " " << cpasma << " " << psimag2 << " " << xdum << " " << rmaxis << " " << xdum << std::endl; 
  eqdsk_out << " " << zmaxis << " " << xdum << " " << psibdry2 << " " << xdum << " " << xdum << std::endl; 
  
  // Write out fpol
  for (int i=0; i<nw; i++)
  {
    element = fpol[i];
    counter = write_line_out(eqdsk_out, element, counter);
  }

  // Write out pres
  for (int i=0; i<nw; i++)
  {
    element = pres[i];
    counter = write_line_out(eqdsk_out, element, counter);
  }

  // Write out ffprime
  for (int i=0; i<nw; i++)
  {
    element = ffprime[i];
    counter = write_line_out(eqdsk_out, element, counter);
  }

  // Write out pprime
  for (int i=0; i<nw; i++)
  {
    element = pprime[i];
    counter = write_line_out(eqdsk_out, element, counter);
  }

  // Write out psi(R,Z) 
  for (int i=0; i<nw; i++)
  {
    for(int j=0; j<nh; j++)
    {
      element = psi[i][j];
      counter = write_line_out(eqdsk_out, element, counter);
    }
  } 
}

// Write singular line out in EQDSK format
int equData::write_line_out(std::ofstream &file, double element, int counter)
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
