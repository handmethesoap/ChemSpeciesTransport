
//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================
#include "SORSolver.hh"
//#include <bits/mathcalls.h>
#include <Debug.hh>
#include <cmath>
#include <iostream>

SORSolver::SORSolver(const int itermax, const int checkfrequency, const int normalizationfrequency, const real eps, const real omg) : itermax_(itermax), checkfrequency_(checkfrequency), normalizationfrequency_(normalizationfrequency), eps_(eps), omg_(omg) {}

SORSolver::SORSolver(const FileReader& configuration)
{
    itermax_ = configuration.getIntParameter("itermax");
    CHECK_MSG((itermax_ >= 0), "itermax is negative. Check your configuration.")
    eps_ = configuration.getRealParameter("eps");
    CHECK_MSG((eps_ >= 0), "eps (residual bound) is negative. Check your configuration.")
    omg_ = configuration.getRealParameter("omg");
    CHECK_MSG((omg_ > 0) && (omg_ < 2), "omg (SOR paramater) is not between (0,2). Check your configuration.")  
    checkfrequency_ = configuration.getIntParameter("checkfrequency");
    CHECK_MSG((checkfrequency_ >= 0), "checkfrequency is negative. Check your configuration.")
    normalizationfrequency_ = configuration.getIntParameter("normalizationfrequency");
    CHECK_MSG((normalizationfrequency_ >= 0), "normalizationfrequency is negative. Check your configuration.")
}

void SORSolver::registerModule(FileReader& configuration)
{
   configuration.registerIntParameter("checkfrequency", 0);
   configuration.registerIntParameter("normalizationfrequency", 0);
}


bool SORSolver::solve(StaggeredGrid& grid)
{
    Array<real> & p = grid.p();
    Array<real> & rhs = grid.rhs();
    real dx = grid.dx();
    real dy = grid.dy();    
    int imax = grid.imax();
    int jmax = grid.jmax();

    CHECK_MSG((imax > 0) && (jmax >0), "Domain too small, imax and jmax zero or negative.");

    real residual = 0.0;

    for (int iter = 1; iter <= itermax_; iter++) {
        
	    // Copy the values of the inner points to the boundary points
        for (int i = 1; i <= imax; i++) {
            p(i,0) = p(i,1);
            p(i,jmax+1) = p(i,jmax);
        }
        for (int j = 1; j <= jmax; j++) {
            p(0,j) = p(1,j);
            p(imax+1,j) = p(imax,j);
        }
        
        // Perform normalization
        if ((normalizationfrequency_ != 0) && (iter % normalizationfrequency_ == 0)) {

            real average = 0.0;
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (!grid.isFluid(i,j)) continue;
                    average += p(i,j);
                }
            }
            average /= grid.getNumFluid();
            for (int j = 1; j <= jmax; j++) {
                for (int i = 1; i <= imax; i++) {
                    if (!grid.isFluid(i,j)) continue;
                    p(i,j) -= average;
                }
            }
        }
	    
        // Copy the values of the inner points to the boundary points
        for (int i = 1; i <= imax; i++) {
            p(i,0) = p(i,1);
            p(i,jmax+1) = p(i,jmax);
        }
        for (int j = 1; j <= jmax; j++) {
            p(0,j) = p(1,j);
            p(imax+1,j) = p(imax,j);
        }
        
        // Calculate residual
        if ((checkfrequency_ != 0) && (iter % checkfrequency_ == 0)) {
          residual = 0.0;
          for(int j = 1; j <= jmax; j++) {
            for(int i = 1; i <= imax; i++) {
                if (!grid.isFluid(i,j)) continue;
                real temp = (grid.p(i,j,EAST)+grid.p(i,j,WEST))/(dx*dx) - (2/(dx*dx) + 2/(dy*dy)) * p(i,j) + (grid.p(i,j,NORTH)+grid.p(i,j,SOUTH))/(dy*dy) - rhs(i,j);
                residual += temp * temp;
            }
          }
          residual = std::sqrt(residual/grid.getNumFluid());
    
          std::stringstream is;
          is  << "Iteration: " << iter << " Residual: " << residual;// << std::endl;
          PRG_MSG(is.str())
          if (residual < eps_) break;
        }

        // Perform SOR on the inner points
       for (int j = 1; j <= jmax; j++) {
            for (int i = (j+1)%2+1; i <= imax; i+=2) {
                if (!grid.isFluid(i,j)) continue;
                p(i,j) = (1-omg_)*p(i,j) + omg_/(2.0/(dx*dx)+2.0/(dy*dy)) *
                         ( (grid.p(i,j,EAST)+grid.p(i,j,WEST))/(dx*dx) + (grid.p(i,j,NORTH)+grid.p(i,j,SOUTH))/(dy*dy) - rhs(i,j) );
            }
        }

        for (int j = 1; j <= jmax; j++) {
            for (int i = j%2+1; i <= imax; i+=2) {
                if (!grid.isFluid(i,j)) continue;
                p(i,j) = (1-omg_)*p(i,j) + omg_/(2.0/(dx*dx)+2.0/(dy*dy)) *
                         ( (grid.p(i,j,EAST)+grid.p(i,j,WEST))/(dx*dx) + (grid.p(i,j,NORTH)+grid.p(i,j,SOUTH))/(dy*dy) - rhs(i,j) );
            }
        }
    }

    return  (residual < eps_);
}
