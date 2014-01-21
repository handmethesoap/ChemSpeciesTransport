#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "Debug.hh"
#include <string>
#include <functional>
#include "VTKWriter.hh"

class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid() {return sg_;}
      const StaggeredGrid & grid() const {return sg_;}
      
      void computeFG();
      void determineNextDT ( real const & limit );
      void composeRHS();
      void updateVelocities();
      void refreshBoundaries();
      void assignInitialValues();
      void simulateGeneral(std::function<bool(real,unsigned int)> criterion);

      static void registerModule(FileReader & params);
  private:

  protected:
      real dt_;
      real Re_;
      real gamma_;
      real gx_;
      real gy_;
      real safetyfactor_;

      std::string boundary_condition_S_;
      std::string boundary_condition_N_;
      std::string boundary_condition_W_;
      std::string boundary_condition_E_;

      real boundary_velocity_S_;
      real boundary_velocity_N_;
      real boundary_velocity_W_;
      real boundary_velocity_E_;
      
      real U_INIT_;
      real V_INIT_;
      real P_INIT_;
      
      StaggeredGrid sg_;
      SORSolver sor_;
      
      VTKWriter vtkWriter_;
      int outputinterval_;
};



#endif
