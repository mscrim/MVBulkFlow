CC=g++

FFLAGS = -03
LFLAGS = -lgsl -lblas
# May be needed instead of the above: LFLAGS = -lgsl -lgslcblas

CPPFILES = usefulfuncs.cpp cosmology.cpp MV.cpp MLE.cpp MV_bulkflow.cpp
OBJECTS = $(CPPFILES:.cpp=.o)

MV_bulkflow: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o MV_bulkflow

%.o: %.cpp
	$(CC) $(FFLAGS) -c $(CPPFILES)

clean:
	rm -rdf *.o MV_bulkflow