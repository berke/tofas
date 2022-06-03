# tofas

This is a **partial** translation and an adaptation of selected
subroutines from the IAU SOFA library, based on the Fortran version.

Using pure Rust code has benefits such as better optimization, no
compilation or linking headaches and the safety guarantees from the
compiler.

I also provide a nicer interface than the Fortran (or C) routines.
For example the different kinds of time (TAI, TT, UT1, etc.)
are wrapped in separate types.

The subroutines are sufficient for LEO Earth observation evaluation
purposes, in particular for computing the Sun angles.

Please note that the `cio_locator` function has not been translated
yet and returns 0.

The dependencies have been kept to a minimum.

## Author

Berke DURAK <bd@exhrd.fr>
