CC          = g++ 

#-----------------------------------------
#Optimization ----------------------------
OPT   =  -O3

#-----------------------------------------

TARGETS = FluidSimulator

OBJECTS = jsoncpp.o Particle.o PressureSolver.o SimulationParameters.o StaggeredGrid.o

#-----------------------------------------

CCOPTS = $(OPT) -I/usr/include/eigen3/
LDOPTS = 

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

FluidSimulator: $(OBJECTS) FluidSimulator.o
	$(CC) $(OBJECTS) $(LDOPTS) FluidSimulator.o -o FluidSimulator

#-----------------------------------------
#-----------------------------------------

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

#-----------------------------------------
#-----------------------------------------















