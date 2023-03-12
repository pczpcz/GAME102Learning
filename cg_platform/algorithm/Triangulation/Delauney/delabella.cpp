/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

//#define DELABELLA_AUTOTEST

// in case of troubles, allows to see if any assert pops up.
// define it globally (like with -DDELABELLA_AUTOTEST)

// this is to work around bug in the DekkersProduct()
// note Splitter constant is also wrong but not used outside DekkersProduct()

// benching hack, fixme!

#include <stdint.h>

uint64_t sorting_bench = 0;



