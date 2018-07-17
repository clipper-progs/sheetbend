// Header file for shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */


#include "shiftfield.h"


bool Shift_field_refine::shift_field_coord( const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
                                            clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
                                            float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid), cffty(grid), cfftz(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    std::complex<float> i(0.0,1.0);
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const std::complex<float> cdata = i * cfftx.cplx_data(c);
          cfftx.cplx_data(c) = float(clipper::Util::twopi()*hkl.h()) * cdata;
          cffty.cplx_data(c) = float(clipper::Util::twopi()*hkl.k()) * cdata;
          cfftz.cplx_data(c) = float(clipper::Util::twopi()*hkl.l()) * cdata;
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cffty.fft_h_to_x( 1.0 / cmap.cell().volume() );
    cfftz.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        	x2map[iw] = cffty.real_data( iw.coord() );
        	x3map[iw] = cfftz.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    mapout.open_write( "grad2.map" );
    mapout.export_xmap( x2map );
    mapout.close_write();
    mapout.open_write( "grad3.map" );
    mapout.export_xmap( x3map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> mmap( spgr, cell, grid );
  clipper::Xmap<float> y0map( spgr, cell, grid );
  clipper::Xmap<float> y1map( spgr, cell, grid );  // const clipper::Xmap<float> y1map = x1map;
  clipper::Xmap<float> y2map( spgr, cell, grid );  // const clipper::Xmap<float> y2map = x2map;
  clipper::Xmap<float> y3map( spgr, cell, grid );  // const clipper::Xmap<float> y3map = x3map;
  clipper::Xmap<float> x00map( spgr, cell, grid );
  clipper::Xmap<float> x01map( spgr, cell, grid );
  clipper::Xmap<float> x02map( spgr, cell, grid );
  clipper::Xmap<float> x03map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );
  clipper::Xmap<float> x12map( spgr, cell, grid );
  clipper::Xmap<float> x13map( spgr, cell, grid );
  clipper::Xmap<float> x22map( spgr, cell, grid );
  clipper::Xmap<float> x23map( spgr, cell, grid );
  clipper::Xmap<float> x33map( spgr, cell, grid );

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) mmap[ix]   = mask[ix];

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = 1.0      *1.0      ;
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = 1.0      *x1map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x02map[ix] = 1.0      *x2map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x03map[ix] = 1.0      *x3map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = x1map[ix]*x2map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = x1map[ix]*x3map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = x2map[ix]*x2map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = x2map[ix]*x3map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = x3map[ix]*x3map[ix];

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = 1.0      *ymap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = x2map[ix]*ymap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = x3map[ix]*ymap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // mask maps
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = mmap[ix]*y0map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = mmap[ix]*y1map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y2map[ix]  = mmap[ix]*y2map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y3map[ix]  = mmap[ix]*y3map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = mmap[ix]*x00map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = mmap[ix]*x01map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x02map[ix] = mmap[ix]*x02map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x03map[ix] = mmap[ix]*x03map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = mmap[ix]*x11map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x12map[ix] = mmap[ix]*x12map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x13map[ix] = mmap[ix]*x13map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x22map[ix] = mmap[ix]*x22map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x23map[ix] = mmap[ix]*x23map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x33map[ix] = mmap[ix]*x33map[ix];

  fltr( mmap, mmap );
  fltr( y0map, y0map );
  fltr( y1map, y1map );
  fltr( y2map, y2map );
  fltr( y3map, y3map );
  fltr( x00map, x00map );
  fltr( x01map, x01map );
  fltr( x02map, x02map );
  fltr( x03map, x03map );
  fltr( x11map, x11map );
  fltr( x12map, x12map );
  fltr( x13map, x13map );
  fltr( x22map, x22map );
  fltr( x23map, x23map );
  fltr( x33map, x33map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {
    std::vector<float> v(4);
    v[0] = y0map[ix];
    v[1] = y1map[ix];
    v[2] = y2map[ix];
    v[3] = y3map[ix];
    clipper::Matrix<float> m(4,4);
    m(0,0) = x00map[ix];
    m(1,1) = x11map[ix];
    m(2,2) = x22map[ix];
    m(3,3) = x33map[ix];
    m(0,1) = m(1,0) = x01map[ix];
    m(0,2) = m(2,0) = x02map[ix];
    m(0,3) = m(3,0) = x03map[ix];
    m(1,2) = m(2,1) = x12map[ix];
    m(1,3) = m(3,1) = x13map[ix];
    m(2,3) = m(3,2) = x23map[ix];
    v = m.solve(v);
    x1map[ix] = v[1];
    x2map[ix] = v[2];
    x3map[ix] = v[3];
  }

  return true;
}



bool Shift_field_refine::shift_field_u_iso( const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
                                            clipper::Xmap<float>& x1map,
                                            float rad, int filter )
{
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  const clipper::Spacegroup&    spgr = clipper::Spacegroup::p1();
  const clipper::Cell&          cell = cmap.cell();
  const clipper::Grid_sampling& grid = cmap.grid_sampling();
  const clipper::Xmap<float>&   ymap = dmap;

  { // being map preparation
    clipper::FFTmap_p1 cfftx(grid);
    clipper::Xmap<float>::Map_reference_coord i0( cmap, clipper::Coord_grid(0,0,0) );
    clipper::Xmap<float>::Map_reference_coord iu, iv, iw;

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	cfftx.real_data( iw.coord() ) = cmap[iw];
        }

    // calculate map coefficients
    cfftx.fft_x_to_h( cmap.cell().volume() );

    // calculate gradient map coefficients
    const clipper::Grid_sampling& g = cfftx.grid_real();
    const clipper::Grid&         gr = cfftx.grid_reci();
    clipper::Coord_grid ch( g.nu()/2, g.nv()/2, g.nw()/2 );
    clipper::Coord_grid c;
    for ( c.u() = 0; c.u() < gr.nu(); c.u()++ )
      for ( c.v() = 0; c.v() < gr.nv(); c.v()++ )
        for ( c.w() = 0; c.w() < gr.nw(); c.w()++ ) {
          const clipper::HKL hkl = clipper::HKL( ( c + ch ).unit( g ) - ch );
          const float scl = clipper::Util::twopi2() * hkl.invresolsq(cell);
          cfftx.cplx_data(c) = scl * cfftx.cplx_data(c);
        }

    // calculate gradient maps
    cfftx.fft_h_to_x( 1.0 / cmap.cell().volume() );

    // copy map
    for ( iu = i0; iu.coord().u() < grid.nu(); iu.next_u() )
      for ( iv = iu; iv.coord().v() < grid.nv(); iv.next_v() )
        for ( iw = iv; iw.coord().w() < grid.nw(); iw.next_w() ) {
        	x1map[iw] = cfftx.real_data( iw.coord() );
        }

    /*
    // output intermediate maps
    clipper::CCP4MAPfile mapout;
    mapout.open_write( "grad1.map" );
    mapout.export_xmap( x1map );
    mapout.close_write();
    */

  } // end map preparation

  // make xmap
  clipper::Xmap<float> mmap( spgr, cell, grid );
  clipper::Xmap<float> y0map( spgr, cell, grid );
  clipper::Xmap<float> y1map( spgr, cell, grid );  // const clipper::Xmap<float> y1map = x1map;
  clipper::Xmap<float> x00map( spgr, cell, grid );
  clipper::Xmap<float> x01map( spgr, cell, grid );
  clipper::Xmap<float> x11map( spgr, cell, grid );

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) mmap[ix]   = mask[ix];

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = 1.0      *1.0      ;
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = 1.0      *x1map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = x1map[ix]*x1map[ix];

  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = 1.0      *ymap[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = x1map[ix]*ymap[ix];

  // filter maps
  clipper::MapFilterFn_step      f0(rad);
  clipper::MapFilterFn_linear    f1(rad);
  clipper::MapFilterFn_quadratic f2(rad);
  clipper::MapFilter_fft<float> fltr( f2, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 0 ) fltr = clipper::MapFilter_fft<float>( f0, 1.0, clipper::MapFilter_fft<float>::Relative );
  if ( filter == 1 ) fltr = clipper::MapFilter_fft<float>( f1, 1.0, clipper::MapFilter_fft<float>::Relative );

  // mask maps
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y0map[ix]  = mmap[ix]*y0map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) y1map[ix]  = mmap[ix]*y1map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x00map[ix] = mmap[ix]*x00map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x01map[ix] = mmap[ix]*x01map[ix];
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) x11map[ix] = mmap[ix]*x11map[ix];

  fltr( mmap, mmap );
  fltr( y0map, y0map );
  fltr( y1map, y1map );
  fltr( x00map, x00map );
  fltr( x01map, x01map );
  fltr( x11map, x11map );

  // calculate U shifts
  for ( MRI ix = ymap.first(); !ix.last(); ix.next() ) {
    std::vector<float> v(2);
    v[0] = y0map[ix];
    v[1] = y1map[ix];
    clipper::Matrix<float> m(2,2);
    m(0,0) = x00map[ix];
    m(1,1) = x11map[ix];
    m(0,1) = m(1,0) = x01map[ix];
    v = m.solve(v);
    x1map[ix] = v[1];
  }

  return true;
}

