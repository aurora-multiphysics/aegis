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


