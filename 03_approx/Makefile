
OBJS1 = numtype.o cheby.o test_ch.o# object files

PROG1 = cheb # code name

F90 = gfortran

F90FLAGS = -O3 -funroll-loops # -fexternal-blas # optimization 

#LIBS = -llapack -lblas #-framework Accelerate # library 

LDFLAGS = $(LIBS)

all: $(PROG1) 

$(PROG1): $(OBJS1)
	$(F90) $(LDFLAGS) -o $@ $(OBJS1) -llapack -lblas

clean:
	rm -f $(PROG1) *.{o,mod} fort.*

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<



