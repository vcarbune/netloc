# Basic makefile build on top of snap
TARGETS=simple

include makefile.config

all: $(TARGETS)

simple: simple.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o *.dat $(TARGETS)
