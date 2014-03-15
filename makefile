# Basic makefile build on top of snap
include makefile.config

TARGETS=netloc_mpi

all: $(TARGETS)

netloc_mpi: netloc_mpi.cpp master_mpi.o worker_mpi.o utils.o hypothesis.o epfl_solver.o ec2.h $(CSNAP)/Snap.o
	$(MPICXX) $(CXXFLAGS) $^ -isystem $(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%_mpi.o: %_mpi.cpp %_mpi.h
	$(MPICXX) $(CXXFLAGS) $< -c -isystem $(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp %.h
	$(CC) $(CXXFLAGS) $< -c -isystem $(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -isystem $(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

generators/synthetic_generator: hypothesis.o

clean:
	rm -f *.o $(TARGETS)
