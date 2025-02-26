CC          = clang++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -O3
#OPT   =  -DDEBUG

#-----------------------------------------

TARGETS = Array3DTest StaggeredGridTest FluidSimulator ParticleViewer

OBJECTS = jsoncpp.o Particle.o PressureSolver.o SimulationParameters.o StaggeredGrid.o

#-----------------------------------------

CCOPTS = $(OPT) -g -std=c++17 -I/opt/homebrew/include/eigen3 -target arm64-apple-macos11
LDOPTS = -std=c++17 -g -target arm64-apple-macos11 -framework OpenGL -framework GLUT -L/opt/homebrew/lib

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

Array3DTest: Array3DTest.o
	$(CC) $(LDOPTS) Array3DTest.o -o Array3DTest

StaggeredGridTest: $(OBJECTS) StaggeredGridTest.o
	$(CC) $(LDOPTS) $(OBJECTS) StaggeredGridTest.o -o StaggeredGridTest

FluidSimulator: $(OBJECTS) FluidSimulator.o
	$(CC) $(OBJECTS) $(LDOPTS) FluidSimulator.o -o FluidSimulator


ParticleViewer: ParticleViewer.o
	$(CC) $(LDOPTS) ParticleViewer.o -o ParticleViewer 
#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------
