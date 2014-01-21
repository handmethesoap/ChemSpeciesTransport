#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include <FileReader.hh>

class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver ( const int itermax, const int checkfrequency, const int normalizationfrequency, const real eps, const real omg );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );


   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );
   
   // register parameters in configuration file
   static void registerModule( FileReader & configuration );

private:
   int itermax_;
   int checkfrequency_;
   int normalizationfrequency_;
   real eps_;
   real omg_;
  
};

#endif //SOR_SOLVER_HH
