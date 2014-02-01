#include "Types.hh"
#include "StaggeredGrid.hh"
#include "FluidSimulator.hh"
#include <Debug.hh>
#include <random>
#include <iostream>
#include "FileReader.hh"
#include "SORSolver.hh"

int main(){
    
  FileReader params;
  FluidSimulator::registerModule(params);
  CHECK(params.readFile("./concentration.par"));
  FluidSimulator fsim(params);
  
  for(int i = 0; i < fsim.grid().numSpecies(); ++i){
    fsim.grid().c(i).print();
  }
  
}