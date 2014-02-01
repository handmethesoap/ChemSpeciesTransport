#include "Types.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"
#include <Debug.hh>
#include <random>
#include <iostream>
#include "FileReader.hh"
#include "SORSolver.hh"

int main(){
  
  using namespace std;
  FileReader params;
  FluidSimulator::registerModule(params);
  CHECK(params.readFile("concentration2.par"));
  FluidSimulator fsim(params);
  
  int timesteps = params.getIntParameter("timesteps");
  CHECK_MSG(timesteps > 0, "timesteps non-positive. Check your config.")
  fsim.simulateTimeStepCount(timesteps);
  
  return 1;
  
}