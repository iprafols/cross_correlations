# lines starting with '#' are comments
# usage:    make "target"
# "target" specifies the task to be performed. Default is all


#
# VARIABLES:
#


# this variable sets the compiler to use
CC = g++

# this variable sets the library headers files' path
INCLUDE = -I/usr/local/include/CCfits -I/opt/ccfits/include/CCfits
#INCLUDE = -I/opt/ccfits/include/CCfits

# this variable sets the library name and path
LIBRARY = -lCCfits -lcfitsio -L/usr/local/lib/
#LIBRARY = -lCCfits -lcfitsio -L/usr/local/lib/ -L/opt/ccfits/lib/

# this variables sets the name of the executable
EXECUTABLE_CORRELATION = programs/correlation.run
EXECUTABLE_PLATE_NEIGHBOURS = programs/plates.run

# this variable contains the list of sources
SOURCES = $(wildcard *.cpp)
SOURCES_CORRELATION = main_correlation.cpp astro_object.cpp astro_object_dataset.cpp correlation_plate.cpp correlation_results.cpp dataset.cpp function_compute_plate_neighbours.cpp input.cpp interpolation_map.cpp lya_pixel.cpp lya_spectra_dataset.cpp lya_spectrum.cpp plate.cpp plate_neighbours.cpp plots_object.cpp sphere_point.cpp
SOURCES_PLATE_NEIGHBOURS = main_plates.cpp global_variables.cpp plate.cpp plate_neighbours.cpp sphere_point.cpp
SOURCES_PLOT = $(wildcard ./output/plots/*.py)

# this variable contains the list of object files
OBJECTS = $(patsubst %.cpp,build/%.o,$(SOURCES))
OBJECTS_CORRELATION = $(patsubst %.cpp,build/%.o,$(SOURCES_CORRELATION))
OBJECTS_PLATE_NEIGHBOURS = $(patsubst %.cpp,build/%.o,$(SOURCES_PLATE_NEIGHBOURS))
OBJECTS_PLOT = $(patsubst %.py,%.pyc,$(SOURCES_PLOT))

# this variables sets the options passed to the compiler for compilation only
CFLAGS = -c $(INCLUDE) $(LIBRARY)
# this variable sets the options passed to the compiler for linking
LFLAGS = $(INCLUDE) $(LIBRARY) -fopenmp -Wall


#
# TARGETS:
#

all: correlation

correlation: $(EXECUTABLE_CORRELATION)

plates: $(OBJECTS_PLATE_NEIGHBOURS) ../plate_neighbours.dat

plots: $(OBJECTS_PLOT)



#
# RUN TARGETS
#
../plate_neighbours.dat: $(EXECUTABLE_PLATE_NEIGHBOURS)
	./$(EXECUTABLE_PLATE_NEIGHBOURS)



#
# COMPILATION TARGETS
#
build/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

%.pyc: %.py
	python $<

$(EXECUTABLE_CORRELATION): $(OBJECTS_CORRELATION)
	$(CC) $(LFLAGS) $(OBJECTS_CORRELATION) -o $(EXECUTABLE_CORRELATION)
	#./$(EXECUTABLE_CORRELATION)

$(EXECUTABLE_PLATE_NEIGHBOURS): $(OBJECTS_PLATE_NEIGHBOURS)
	$(CC) $(LFLAGS) $(OBJECTS_PLATE_NEIGHBOURS) -o $(EXECUTABLE_PLATE_NEIGHBOURS)



#
# CLEAN TARGETS
#
clean: clean_objects clean_programs

clean_objects:
	-rm build/*.o

clean_programs:
	-rm programs/*.run

clean_plots:
	-rm ../plots/*.pyc


