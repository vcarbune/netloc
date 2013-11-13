# Basic makefile build on top of snap
include makefile.config

TARGETS=generator simple netloc

all: $(TARGETS)

netloc: netloc.cpp utils.o hypothesis.o ec2.h $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) $^ -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp %.h
	$(CC) $(CXXFLAGS) $< -c -I$(CGLIB) $(LDFLAGS) $(LIBS) -o $@

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o *.dat $(TARGETS)
