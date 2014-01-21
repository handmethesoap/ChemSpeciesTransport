#include "Array.hh"
#include "FluidSimulator.hh"
#include "FileReader.hh"

#include <iostream>
#include <algorithm>

int main( int argc, char** argv )
{
    using namespace std;
    FileReader params;
    FluidSimulator::registerModule(params);
    CHECK(params.readFile("backstep.par"));
    FluidSimulator fsim(params);
    
    int timesteps = params.getIntParameter("timesteps");
    CHECK_MSG(timesteps > 0, "timesteps non-positive. Check your config.")
    fsim.simulateTimeStepCount(timesteps);
    
    return 0;
}
