# Basic makefile build on top of snap
include makefile.config

TARGETS=generator simple netloc

all: $(TARGETS)

netloc: hypothesis.o ec2.h netloc.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o *.dat $(TARGETS)
