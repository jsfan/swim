#FC=gfortran
#FCFLAGS=-Wall -fdefault-real-8
FC=g95
FCFLAGS=-Wall -r8

LDFLAGS=-lnetcdf -lnetcdff -lsigwatch -L/usr/local/lib -L/usr/lib
INC=-I/usr/local/include -I/usr/include

PROGRAMS=swim swim-debug

OBJS=grid.o init.o output.o localconfig.o sighandler.o helpers.o integration.o input.o
OBJS_DEBUG=grid.o init.o output_debug.o localconfig.o sighandler.o helpers.o integration.o input.o 

.PHONY: debug default all

default: $(OBJS) swim

debug: $(OBJS_DEBUG) swim-debug

all: $(PROGRAMS)

swim-debug: $(OBJS_DEBUG) swim.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(INC)

swim: $(OBJS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) $(INC)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $< $(INC)

clean:
	rm -f *.o *.mod $(PROGRAMS)
