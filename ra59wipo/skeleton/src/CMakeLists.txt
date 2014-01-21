
include_directories( . )

add_executable( nusif VTKWriter.cc FluidSimulator.cc StaggeredGrid.cc SORSolver.cc FileReader.cc Debug.cc GrayScaleImage.cc lodepng.cc main.cc )
file(COPY ../dcavity.par DESTINATION .)
file(COPY ../basin.par DESTINATION .)
file(COPY ../backstep.par DESTINATION .)
