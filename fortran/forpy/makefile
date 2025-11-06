FC=gfortran
FFLAGS=-Wall -g -O2
PY_VERSION=3

EXE=ED

ifeq ($(FC),ifort)
  FFLAGS=-warn all -g -O2
endif

PYTHON=python3

ifeq ($(PY_VERSION),2)
  PYTHON=python2
  FFLAGS+=-DPYTHON2
endif

PYTHON_CONFIG=$(PYTHON)-config

ifeq ($(PY_DEBUG), 1)
  PYTHON_CONFIG=$(PYTHON)-dbg-config
  FFLAGS+= -DPy_DEBUG
endif

LDFLAGS=`$(PYTHON_CONFIG) --ldflags --embed`

.PHONY: all clean

all: $(EXE)

$(EXE): ED.o forpy_mod.o
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE)

ED.o: forpy_mod.o