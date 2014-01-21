#include "Types.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"
#include <Debug.hh>
#include <random>
#include <iostream>
#include "FileReader.hh"
#include "SORSolver.hh"

int main()
{
   using namespace std;
   FileReader params;
   FluidSimulator::registerModule(params);
   CHECK(params.readFile("./compfg.par"));
   FluidSimulator fsim(params);

   auto &u = fsim.grid().u();
   auto &v = fsim.grid().v();
   auto &f = fsim.grid().f();
   auto dx = fsim.grid().dx();
   auto dy = fsim.grid().dy();
   auto dt = params.getRealParameter("dt");
   auto Re = params.getRealParameter("Re");
   auto imax = params.getIntParameter("imax");
   auto jmax = params.getIntParameter("jmax");
   
   // u = x²y²
   u.fill([dx,dy](int i, int j, int k) -> real {return i * i * dx * dx * (j-0.5) * (j-0.5) * dy * dy;});
   // v = x²y²
   v.fill([dx,dy](int i, int j, int k) -> real {return (i-0.5) * (i-0.5) * dx * dx * j * j * dy * dy;});

   fsim.computeFG();
   
   Array<real> test(imax+1, jmax+2);
   // test = x²y² + dt * (2/Re (x² + y²) - 4y⁴x³ - 4x⁴y³ + gx)
   test.fill([dx,dt,dy,Re](int i, int j, int k) -> real { return i * i * dx * dx * (j-0.5) * (j-0.5) * dy * dy + dt * (
             2.0/Re * ((j-0.5) * (j-0.5) * dy * dy + i * i * dx * dx) - 
             4.0 * (j-0.5) * (j-0.5) * (j-0.5) * (j-0.5) * dy * dy * dy * dy * i * i * i * dx * dx * dx - 
             4.0 * i * i * i * i * dx * dx * dx * dx * (j-0.5) * (j-0.5) * (j-0.5) * dy * dy * dy); } );
   
   real maxerror = 0.0;
   for (int j = 1; j <= jmax; j++) {
      for (int i = 1; i < imax; i++) {
        maxerror = max(maxerror, fabs(test(i,j) - f(i,j)));
      }
   }
   
   cout << "Checking f" << endl;
   CHECK(maxerror < 0.0001)
   cout << "CHECK passed" << endl;;
   return 0;
}
