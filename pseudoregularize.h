// Header file for shift field refinement
/* Copyright 2020 Kevin Cowtan & University of York all rights reserved */


#ifndef PSEUDOREGULARIZE
#define PSEUDOREGULARIZE

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>

class PseudoRegularize{

 public:
  static bool regularize( clipper::MiniMol& mol1, const clipper::MiniMol& mol2 );
  static bool regularize( clipper::MPolymer& mp1, const clipper::MPolymer& mp2 );
};

#endif

