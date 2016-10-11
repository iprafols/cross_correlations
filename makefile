# lines starting with '#' are comments
# usage:    make "target"
# "target" specifies the task to be performed. Default is all


#
# VARIABLES:
#


# this variable sets the compiler to use
CC = g++

# this variable sets the library headers files' path
INCLUDE = -I/usr/local/include/CCfits -I/opt/ccfits/include/CCfits -I/./lib
#INCLUDE = -I/opt/ccfits/include/CCfits

# this variable sets the library name and path
LIBRARY = -lCCfits -lcfitsio -L/usr/local/lib/
#LIBRARY = -lCCfits -lcfitsio -L/usr/local/lib/ -L/opt/ccfits/lib/

# this variables sets the name of the executable
EXECUTABLE_CORRELATION = programs/correlation.run
EXECUTABLE_LYA1D = programs/compute_lya1d.run
EXECUTABLE_PROJECTION_CORRECTION = programs/compute_projection_correction.run
EXECUTABLE_TEST_LYA_DELTAS = programs/test_lya_deltas.run

# this variable contains the list of sources
SOURCES = $(wildcard *.cpp)
SOURCES_CORRELATION = main_correlation.cpp astro_object.cpp astro_object_dataset.cpp civ_spectra_dataset.cpp correlation_plate.cpp correlation_results.cpp covariance_matrix.cpp covariance_plate.cpp dataset.cpp distortion_matrix.cpp distortion_plate.cpp dla_dataset.cpp function_compute_plate_neighbours.cpp input.cpp interpolation_map.cpp lya_auto_interpolation_map.cpp lya_mean_projected_deltas_interpolation_map.cpp lya_pixel.cpp lya_spectra_dataset.cpp lya_spectrum.cpp pair.cpp pair_dataset.cpp plate.cpp plate_neighbours.cpp plots_object.cpp quasar_dataset.cpp spectra_dataset.cpp sphere_point.cpp z_dist_interpolation_map.cpp
SOURCES_LYA1D = main_compute_lya_1d.cpp astro_object.cpp astro_object_dataset.cpp civ_spectra_dataset.cpp function_compute_plate_neighbours.cpp input.cpp interpolation_map.cpp lya_auto_interpolation_map.cpp lya_pixel.cpp lya_spectra_dataset.cpp lya_spectrum.cpp plate.cpp plate_neighbours.cpp plots_object.cpp spectra_dataset.cpp sphere_point.cpp
SOURCES_PROJECTION_CORRECTION = main_compute_projection_correction.cpp astro_object.cpp astro_object_dataset.cpp civ_spectra_dataset.cpp function_compute_plate_neighbours.cpp input.cpp interpolation_map.cpp lya_auto_interpolation_map.cpp lya_pixel.cpp lya_spectra_dataset.cpp lya_spectrum.cpp plate.cpp plate_neighbours.cpp plots_object.cpp spectra_dataset.cpp sphere_point.cpp
SOURCES_PLOT = $(wildcard ./output/plots/*.py)
SOURCES_TEST_LYA_DELTAS = test_lya_deltas.cpp astro_object.cpp astro_object_dataset.cpp civ_spectra_dataset.cpp function_compute_plate_neighbours.cpp input.cpp interpolation_map.cpp lya_auto_interpolation_map.cpp lya_pixel.cpp lya_spectra_dataset.cpp lya_spectrum.cpp plate.cpp plate_neighbours.cpp plots_object.cpp spectra_dataset.cpp sphere_point.cpp

# this variable contains the list of object files
OBJECTS = $(patsubst %.cpp,build/%.o,$(SOURCES))
OBJECTS_CORRELATION = $(patsubst %.cpp,build/%.o,$(SOURCES_CORRELATION))
OBJECTS_LYA1D = $(patsubst %.cpp,build/%.o,$(SOURCES_LYA1D))
OBJECTS_PROJECTION_CORRECTION = $(patsubst %.cpp,build/%.o,$(SOURCES_PROJECTION_CORRECTION))
OBJECTS_PLATE_NEIGHBOURS = $(patsubst %.cpp,build/%.o,$(SOURCES_PLATE_NEIGHBOURS))
OBJECTS_PLOT = $(patsubst %.py,%.pyc,$(SOURCES_PLOT))
OBJECTS_TEST_LYA_DELTAS = $(patsubst %.cpp,build/%.o,$(SOURCES_TEST_LYA_DELTAS))

# this variables sets the options passed to the compiler for compilation only
CFLAGS = -c $(INCLUDE) $(LIBRARY) -fopenmp
# this variable sets the options passed to the compiler for linking
#LFLAGS = $(INCLUDE) $(LIBRARY) -Wall
LFLAGS = $(INCLUDE) $(LIBRARY) -fopenmp -Wall


#
# TARGETS:
#

all: correlation lya_1d projection_correction

correlation: $(EXECUTABLE_CORRELATION)

lya_1d: $(EXECUTABLE_LYA1D)

projection_correction: $(EXECUTABLE_PROJECTION_CORRECTION)

test_lya_deltas: $(EXECUTABLE_TEST_LYA_DELTAS)

plots: $(OBJECTS_PLOT)



#
# COMPILATION TARGETS
#
build/%.o: %.cpp
	$(CC) $(CFLAGS) src/$< -o $@

%.pyc: %.py
	python $<

$(EXECUTABLE_CORRELATION): $(OBJECTS_CORRELATION)
	$(CC) $(LFLAGS) $(OBJECTS_CORRELATION) -o $(EXECUTABLE_CORRELATION)

$(EXECUTABLE_LYA1D): $(OBJECTS_LYA1D)
	$(CC) $(LFLAGS) $(OBJECTS_LYA1D) -o $(EXECUTABLE_LYA1D)

$(EXECUTABLE_PROJECTION_CORRECTION): $(OBJECTS_PROJECTION_CORRECTION)
	$(CC) $(LFLAGS) $(OBJECTS_PROJECTION_CORRECTION) -o $(EXECUTABLE_PROJECTION_CORRECTION)


$(EXECUTABLE_TEST_LYA_DELTAS): $(OBJECTS_TEST_LYA_DELTAS)
	$(CC) $(LFLAGS) $(OBJECTS_TEST_LYA_DELTAS) -o $(EXECUTABLE_TEST_LYA_DELTAS)


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


