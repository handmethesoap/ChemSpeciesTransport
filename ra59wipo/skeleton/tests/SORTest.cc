#include "Types.hh"
#include "SORSolver.hh"
#include "StaggeredGrid.hh"
#include <Debug.hh>
#include <random>
#include <iostream>

void initGridSetup1( StaggeredGrid & grid )
{
   // Setup 1:
   //    - grid.p   : init with random values
   //    - grid.rhs : init with zero
   std::default_random_engine generator;
   std::uniform_real_distribution<double> distribution(0.0,10.0);
   
   auto & p = grid.p();
   auto & rhs = grid.rhs();

   for (int j = 0; j < p.getSize(1); j++) {
     for (int i = 0; i < p.getSize(0); i++) {
       p(i,j) = distribution(generator);
       rhs(i,j) = 0.0;
     }
   }
}

void initGridSetup2( StaggeredGrid & grid )
{
   // Setup 2:
   //    - grid.p   : init with random values
   //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
   std::default_random_engine generator;
   std::uniform_real_distribution<real> distribution(0.0,1.0);

   auto & p = grid.p();
   auto & rhs = grid.rhs();
   real dx = grid.dx(); 
   for (int j = 0; j < p.getSize(1); j++) {
     for (int i = 0; i < p.getSize(0); i++) {
       p(i,j) = distribution(generator);
       rhs(i,j) = sin(2 * (i+0.5) * dx * M_PI);
     }
   }
}



int main()
{
   int itermax = 1000;
   int checkfrequency = 5;
   int normalizationfrequency = 5;
   real eps = 1e-4;
   real omg = 1.8;
   int imax = 30;
   int jmax = 30;
   real xlength = 1.0;
   real ylength = 1.0;
      
   StaggeredGrid grid(imax,jmax,xlength/imax,ylength/jmax);
   SORSolver solver(itermax,checkfrequency,normalizationfrequency,eps,omg);

   PRG_MSG("RHS zero")
   initGridSetup1( grid );
   CHECK(solver.solve( grid ));
   
   PRG_MSG("RHS sin")
   initGridSetup2( grid );
   CHECK(solver.solve( grid ));
   
   std::cout << "SOR-Test passed." << std::endl;

   return 0;
}
