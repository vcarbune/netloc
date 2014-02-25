# Basic makefile build on top of snap
include makefile.config

TARGETS=netloc_mpi

all: $(TARGETS)

netloc_mpi: netloc_mpi.cpp master_mpi.o worker_mpi.o utils.o hypothesis.o ec2.h $(CSNAP)/Snap.o
	$(MPICXX) $(CXXFLAGS) $^ -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

netloc: netloc.cpp utils.o hypothesis.o ec2.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) $^ -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%_mpi.o: %_mpi.cpp %_mpi.h
	$(MPICXX) $(CXXFLAGS) $< -c -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp %.h
	$(CC) $(CXXFLAGS) $< -c -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o *.dat $(TARGETS)
