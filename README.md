# tofas

This is a translation and an adaptation of selected subroutines from
the IAU SOFA library, based on the Fortran version.

Using pure Rust code has benefits such as better optimization, no
compilation or linking headaches and the safety guarantees from the
compiler.

I also provide a nicer interface than the Fortran (or C) routines.
For example the different kinds of time (TAI, TT, UT1, etc.)
are wrapped in separate types.

The subroutines are sufficient for LEO Earth observation evaluation
purposes, in particular for computing the Sun angles.

The dependencies have been kept to a minimum.

Test vectors generated using the 2017-04-20 Fortran 77 version
of SOFA are included, as well as the F90 code to generate them.

## Author

Berke DURAK <bd@exhrd.fr>
