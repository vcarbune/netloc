# Basic makefile build on top of snap
include makefile.config

TARGETS=generator simple

all: $(TARGETS)

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -o $@ $^ -I$(CGLIB) $(LDFLAGS) $(LIBS)

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o *.dat $(TARGETS)
