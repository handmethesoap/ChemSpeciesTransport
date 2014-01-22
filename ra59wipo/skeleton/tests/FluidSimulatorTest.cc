#include <iostream>

#include "FluidSimulator.hh"
#include "FileReader.hh"
#include "Debug.hh"

int main(int argc, char** argv) {

    if(argc != 2) {
        std::cerr << "Usage: fluidsimulator_test <inputfile>\n";
        exit(EXIT_FAILURE);
    }

    FileReader reader;

    std::cout << "MOOOOOOOOOOO\n";

    reader.registerIntParameter("imax");
    reader.registerIntParameter("jmax");
    reader.registerIntParameter("itermax");
    reader.registerIntParameter("timesteps");
    reader.registerIntParameter("checkfrequency");
    reader.registerIntParameter("normalizationfrequency");
    reader.registerIntParameter("outputinterval");

    reader.registerIntParameter("numSpecies",2);
    
    reader.registerRealParameter("gx");
    reader.registerRealParameter("gy");
    reader.registerRealParameter("Re");
    reader.registerRealParameter("xlength");
    reader.registerRealParameter("ylength");
    reader.registerRealParameter("dt");
    reader.registerRealParameter("eps");
    reader.registerRealParameter("omg");
    reader.registerRealParameter("gamma");
    reader.registerRealParameter("safetyfactor");
    reader.registerRealParameter("U_INIT");
    reader.registerRealParameter("V_INIT");
    reader.registerRealParameter("P_INIT");
    reader.registerRealParameter("duration", 2.0);

    reader.registerStringParameter("boundary_condition_N", "NOSLIP");
    reader.registerRealParameter("boundary_velocity_N", 0.0);
    reader.registerStringParameter("boundary_condition_E", "NOSLIP");
    reader.registerRealParameter("boundary_velocity_E", 0.0);
    reader.registerStringParameter("boundary_condition_S", "NOSLIP");
    reader.registerRealParameter("boundary_velocity_S", 0.0);
    reader.registerStringParameter("boundary_condition_W", "NOSLIP");
    reader.registerRealParameter("boundary_velocity_W", 0.0);
    reader.registerStringParameter("name");
    
    reader.readFile(argv[1]);

    CHECK_MSG(reader.getIntParameter("imax")>0,"imax must be larger than 0");
    CHECK_MSG(reader.getIntParameter("jmax")>0,"jmax must be larger than 0");
    CHECK_MSG(reader.getIntParameter("itermax")>1,"itermax must be larger than 1");
    CHECK_MSG(fabs(reader.getRealParameter("omg")-1.8) <= 0.1,"omg must be between 1.7 and 1.9");
    CHECK_MSG(reader.getRealParameter("eps")>0.0,"eps must be larger than 0.0");
    CHECK_MSG(reader.getRealParameter("xlength")>0.0, "xlength must be larger than 0.0");
    CHECK_MSG(reader.getRealParameter("ylength")>0.0, "ylength must be larger than 0.0");
    CHECK_MSG(reader.getRealParameter("dt")>0.0, "dt must be larger than 0.0");
  
    FluidSimulator simulator(reader);

    for(int i = 0; i < simulator.grid().numSpecies(); ++i){
        simulator.grid().c(i).fill(i);
    }

    for(int i = 0; i < simulator.grid().numSpecies(); ++i){
        simulator.grid().c(i).print();
    }

    std::cout << "lambda1: " << simulator.grid().lambda(0) << "\nlambda2: " << simulator.grid().lambda(1) << "\n";

    //simulator.simulateTimeStepCount(reader.getIntParameter("timesteps"));

    //simulator.simulate(reader.getRealParameter("duration"));

}
