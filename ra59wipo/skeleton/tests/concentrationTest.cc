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
    
  fsim.updateCBoundaries();
  std::cout <<  "<<<<<<<<<<<<<<<<< Testing Dirichlet Boundary Conditions <<<<<<<<<<<<<<<<<" << std::endl << std::endl; 
  for(int i = 0; i < fsim.grid().numSpecies(); ++i){
    fsim.grid().c(i).print();
  }
  
  std::cout <<  "<<<<<<<<<<<<<<<<< Testing Neumann Boundary Conditions <<<<<<<<<<<<<<<<<" << std::endl << std::endl; 
  
  fsim.grid().c(3) = fsim.grid().c(2);
  fsim.grid().c(3).print();
  fsim.updateCBoundaries();
  fsim.grid().c(3).print();
  
  return 1;
}