SOFA_DIR?=../../sofa_f/f77/src

SOFA_SOURCES=$(addprefix $(SOFA_DIR)/,$(filter-out t_sofa_f.for,$(notdir $(wildcard $(SOFA_DIR)/*.for))))

.PHONY: all

all: epv00_data.rs test_data.rs

epv00_data.rs: dump_epv00
	./dump_epv00 >epv00_data.rs

dump_epv00: dump_epv00.f90 epv00.for
	gfortran -Wall dump_epv00.f90 epv00.for -o dump_epv00

test_data.rs: test_sofa
	./test_sofa >test_data.rs

test_sofa: test_sofa.f90 $(SOFA_SOURCES)
	gfortran -Wall test_sofa.f90 $(SOFA_SOURCES) -o test_sofa

test_sun: test_sun.f90 $(SOFA_SOURCES)
	gfortran -Wall test_sun.f90 $(SOFA_SOURCES) -o test_sun
