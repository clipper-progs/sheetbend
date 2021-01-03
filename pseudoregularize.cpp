// Header file for shift field refinement
/* Copyright 2020 Kevin Cowtan & University of York all rights reserved */


#include "pseudoregularize.h"


bool PseudoRegularize::regularize( clipper::MiniMol& mol1, const clipper::MiniMol& mol2 )
{
  // make a list of work atom groups
  double crad = 2.4;
  std::vector<std::vector<clipper::MAtomIndex> > groups;
  for ( int c = 0; c < mol2.size(); c++ ) {
    std::vector<clipper::MAtomIndex> group;
    for ( int r0 = 0; r0 < mol2[c].size(); r0++ ) {
      // add this monomer
      group.push_back( clipper::MAtomIndex(c,r0,0) );
      // check if next monomer is connected to this
      bool iscon = false;
      int r1 = r0 + 1;
      if ( r1 < mol2[c].size() ) {
        for ( int a0 = 0; a0 < mol2[c][r0].size(); a0++ )
          for ( int a1 = 0; a1 < mol2[c][r1].size(); a1++ )
            if ( ( mol2[c][r0][a0].coord_orth() -
                   mol2[c][r1][a1].coord_orth() ).lengthsq() < crad*crad )
              iscon = true;
      }
      // if not, add this group
      if ( !iscon ) {
        groups.push_back( group );
        group.clear();
      }
    }
  }

  // now loop over groups and regularize in turn
  for ( int i = 0; i < groups.size(); i++ ) {
    if ( groups[i].size() > 0 ) {  // should be redundant?
      clipper::MPolymer mp1, mp2;
      for ( int j = 0; j < groups[i].size(); j++ ) {
        const clipper::MAtomIndex& m = groups[i][j];
        mp1.insert( mol1[m.polymer()][m.monomer()] );
        mp2.insert( mol2[m.polymer()][m.monomer()] );
      }
      regularize( mp1, mp2 );
      for ( int j = 0; j < groups[i].size(); j++ ) {
        const clipper::MAtomIndex& m = groups[i][j];
        mol1[m.polymer()][m.monomer()] = mp1[j];
      }
    }
  }

  return true;
}


bool PseudoRegularize::regularize( clipper::MPolymer& mp1, const clipper::MPolymer& mp2 )
{
  // find key atom coords for mp2, which are used to determine per atom weights
  std::vector<clipper::Coord_orth> keys;
  for ( int r = 0; r < mp2.size(); r++ ) {
	  int a = mp2[r].lookup( " CA ", clipper::MM::ANY );
	  if ( a < 0 ) a = mp2[r].lookup( " C1*", clipper::MM::ANY );
    if ( a < 0 ) a = 0;
    keys.push_back( mp2[r][a].coord_orth() );
  }

  // make table of distances by atom
  std::vector<std::vector<double> > w0, w1, w2;
  for ( int r = 0; r < mp2.size(); r++ ) {
    double d0(1.0e6), d2(1.0e6);
    if ( r-1 >= 0         ) d0 = sqrt((keys[r]-keys[r-1]).lengthsq());
    if ( r+1 < mp2.size() ) d2 = sqrt((keys[r]-keys[r+1]).lengthsq());
    std::vector<double> w0tmp, w1tmp, w2tmp;
    for ( int a = 0; a < mp2[r].size(); a++ ) {
      double r0(1.0e6), r2(1.0e6);   //double r0(1.0e6), r1(1.0e6), r2(1.0e6);
      const clipper::Coord_orth co(mp2[r][a].coord_orth());
      if ( r-1 >= 0         ) r0 = sqrt((co-keys[r-1]).lengthsq());
      //r1 = sqrt((co-keys[r  ]).lengthsq());
      if ( r+1 < mp2.size() ) r2 = sqrt((co-keys[r+1]).lengthsq());
      const double w00 = std::min( 1.0 - r0/d0, 0.5 );
      const double w02 = std::min( 1.0 - r2/d2, 0.5 );
      const double w01 = 1.0 - (w00+w02);
      w0tmp.push_back(w00);
      w1tmp.push_back(w01);
      w2tmp.push_back(w02);
    }
    w0.push_back( w0tmp );
    w1.push_back( w1tmp );
    w2.push_back( w2tmp );
  }

  // now superpose a list of pentamer fragments
  clipper::MPolymer f0(mp1), f1(mp1), f2(mp1);
  const int dr = 2;
  for ( int r = 0; r < mp2.size(); r++ ) {
    std::vector<clipper::Coord_orth> co1, co2;
    for ( int r1 = r-dr; r1 <= r+dr; r1++ ) {
      int r2 = std::max( std::min( r1, mp2.size()-1 ), 0 );
      for ( int a = 0; a < mp1[r2].size(); a++ ) {
        co1.push_back( mp1[r2][a].coord_orth() );
        co2.push_back( mp2[r2][a].coord_orth() );
      }
    }
    clipper::RTop_orth rt( co2, co1 );

    // rebuild
    int r2 = std::max( r-1, 0            );
    int r0 = std::min( r+1, mp2.size()-1 );
    for ( int a = 0; a < mp2[r2].size(); a++ )
      f2[r2][a].set_coord_orth( rt * mp2[r2][a].coord_orth() );
    for ( int a = 0; a < mp2[r ].size(); a++ )
      f1[r ][a].set_coord_orth( rt * mp2[r ][a].coord_orth() );
    for ( int a = 0; a < mp2[r0].size(); a++ )
      f0[r0][a].set_coord_orth( rt * mp2[r0][a].coord_orth() );
  }

  // make weighted combination
  for ( int r = 0; r < mp1.size(); r++ ) {
    for ( int a = 0; a < mp1[r].size(); a++ ) {
      const clipper::Coord_orth co( w0[r][a]*f0[r][a].coord_orth() +
                                    w1[r][a]*f1[r][a].coord_orth() +
                                    w2[r][a]*f2[r][a].coord_orth() );
      mp1[r][a].set_coord_orth( co );
    }
  }

  return true;
}
