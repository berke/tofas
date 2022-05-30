epv00_data.rs: dump_epv00
	./dump_epv00 >epv00_data.rs

dump_epv00: dump_epv00.f90 epv00.for
	gfortran -Wall dump_epv00.f90 epv00.for -o dump_epv00

