
include_directories( ../src )
add_executable( array_test       ../src/Debug.cc ArrayTest.cc      )
add_executable( filereader_test ../src/FileReader.cc ../src/Debug.cc FileReaderTest.cc )
add_executable( sor_test        ../src/Debug.cc ../src/SORSolver.cc  ../src/FileReader.cc ../src/StaggeredGrid.cc SORTest.cc ../src/GrayScaleImage.cc ../src/lodepng.cc )
add_executable( comp_fg_test ../src/VTKWriter.cc CompFGTest.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/FluidSimulator.cc ../src/FileReader.cc ../src/SORSolver.cc ../src/GrayScaleImage.cc ../src/lodepng.cc )
add_executable( concentration_test ../src/VTKWriter.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/FluidSimulator.cc ../src/FileReader.cc ../src/SORSolver.cc ../src/GrayScaleImage.cc ../src/lodepng.cc  concentrationTest.cc)
add_executable( concentration_visualisation_test ../src/VTKWriter.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/FluidSimulator.cc ../src/FileReader.cc ../src/SORSolver.cc ../src/GrayScaleImage.cc ../src/lodepng.cc  concentrationVisualisationTest.cc)

add_executable( fluidsimulator_test ../src/VTKWriter.cc FluidSimulatorTest.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/FluidSimulator.cc ../src/FileReader.cc ../src/SORSolver.cc ../src/GrayScaleImage.cc ../src/lodepng.cc )

file(COPY ../poisson.par DESTINATION .)
file(COPY compfg.par DESTINATION .)
file(COPY FileReaderTestInput.txt DESTINATION .)
