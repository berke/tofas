SOFA_DIR?=../../sofa_f/f77/src

SOFA_SOURCES=$(addprefix $(SOFA_DIR)/,$(filter-out t_sofa_f.for,$(notdir $(wildcard $(SOFA_DIR)/*.for))))
#SOFA_SOURCES=$(filter-out t_sofa_f.for,$(notdir $(wildcard $(SOFA_DIR)/*.for)))

$(info Hello $(SOFA_SOURCES) in $(SOFA_DIR))

.PHONY: all

all: epv00_data.rs test_sofa

epv00_data.rs: dump_epv00
	./dump_epv00 >epv00_data.rs

dump_epv00: dump_epv00.f90 epv00.for
	gfortran -Wall dump_epv00.f90 epv00.for -o dump_epv00

test_sofa: test_sofa.f90 $(SOFA_SOURCES)
	gfortran -Wall test_sofa.f90 $(SOFA_SOURCES) -o test_sofa

