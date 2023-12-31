// Clipper app to perform shift field refinement
/* Copyright 2018 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>

extern "C" {
  #include <stdlib.h>
}

#include "shiftfield.h"
#include "pseudoregularize.h"


struct results_by_cycle{ int cyclerr; int cycle; float resolution; float radius; float rwork; float rfree; int nrefl; int nreflp1; };



int main( int argc, char** argv )
{
  CCP4Program prog( "csheetbend", "0.5.0", "$Date: 2021/04/18" );
  prog.set_termination_message( "Failed" );

  std::cout << std::endl << "Copyright 2018-2020 Kevin Cowtan and University of York." << std::endl << std::endl;
  prog.summary_beg();
  std::cout << "$TEXT:Reference: $$ Please reference $$" << std::endl << std::endl << " 'Shift-field refinement of macromolecular atomic models'" << std::endl << "K. Cowtan, S. Metcalfe and P. Bond (2020) Acta Cryst. D76, 1192-1200." << std::endl << std::endl << "$$" << std::endl;
  prog.summary_end();

  // defaults
  enum ANISO { NONE, FOBS, FCAL };
  clipper::String title;
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo= "NONE";
  clipper::String ipcolfree = "NONE";
  clipper::String ipcolfc= "NONE";
  clipper::String pdbfile= "NONE";
  clipper::String mapfile= "NONE";
  clipper::String pdbmask= "NONE";
  clipper::String pdboutfile = "sheetbend.pdb";
  clipper::String mapoutfile = "NONE";
  clipper::String mtzoutfile = "NONE";
  clipper::String opcol = "shifted";
  clipper::String opxml = "NONE";
  int ncyc   = 1;
  int ncycrr = 1;
  int freeflag = 0;
  bool refxyz    = false;
  bool refuiso   = false;
  bool refuaniso = false;
  bool postrefxyz    = false;
  bool postrefuiso   = false;
  bool postrefuaniso = false;
  bool inclconst = false;
  bool pseudoreg = false;
  double rad    = -1.0;
  double radscl =  5.0;
  double res    = -1.0;
  double u_lo_abs = clipper::Util::b2u(  0.1);
  double u_hi_abs = clipper::Util::b2u(999.9);
  std::vector<double> resbycyc;
  ANISO aniso = NONE;
  int filter = 2;
  int n_refln = 1000;
  int n_param = 10;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcolfree = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcolfc = args[arg];
    } else if ( args[arg] == "-resolution" ) {
      if ( ++arg < args.size() ) res    = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-radius" ) {
      if ( ++arg < args.size() ) rad    = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-radius-scale" ) {
      if ( ++arg < args.size() ) radscl = clipper::String(args[arg]).f();
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) pdbfile = args[arg];
    } else if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) mapfile = args[arg];
    } else if ( args[arg] == "-pdbout" ) {
      if ( ++arg < args.size() ) pdboutfile = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) mapoutfile = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) mtzoutfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-xmlout" ) {
      if ( ++arg < args.size() ) opxml  = args[arg];
    } else if ( args[arg] == "-free-flag" ) {
      if ( ++arg < args.size() ) freeflag = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-coord" ) {
      refxyz = true;
    } else if ( args[arg] == "-u-iso" ) {
      refuiso = true;
    } else if ( args[arg] == "-u-aniso" ) {
      refuaniso = true;
      refuiso = false;
    } else if ( args[arg] == "-postrefine-coord" ) {
      postrefxyz    = true;
    } else if ( args[arg] == "-postrefine-u-iso" ) {
      postrefuiso   = true;
    } else if ( args[arg] == "-postrefine-u-aniso" ) {
      postrefuaniso = true;
      postrefuiso   = false;
    } else if ( args[arg] == "-include-const" ) {
      inclconst = true;
    } else if ( args[arg] == "-pseudo-regularize" ) {
      pseudoreg = true;
    } else if ( args[arg] == "-cycles" ) {
      if ( ++arg < args.size() ) ncyc  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-refine-regularize-cycles" ) {
      if ( ++arg < args.size() ) ncycrr  = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-resolution-by-cycle" ) {
      if ( ++arg < args.size() ) {
        std::vector<clipper::String> rc = clipper::String(args[arg]).split(",");
        for ( int i = 0; i < rc.size(); i++ ) resbycyc.push_back( clipper::String(rc[i]).f() );
      }
    } else if ( args[arg] == "-b-iso-range" ) {
      if ( ++arg < args.size() ) {
        std::vector<clipper::String> br = clipper::String(args[arg]).split(",");
        if ( br.size() > 0 ) u_lo_abs = clipper::Util::b2u(clipper::String(br[0]).f());
        if ( br.size() > 1 ) u_hi_abs = clipper::Util::b2u(clipper::String(br[1]).f());
      }
    } else if ( args[arg] == "-aniso-obs" ) {
      aniso = FOBS;
    } else if ( args[arg] == "-aniso-cal" ) {
      aniso = FCAL;
    } else if ( args[arg] == "-pdbin-mask" ) {
      if ( ++arg < args.size() ) pdbmask = args[arg];
    } else if ( args[arg] == "-filter" ) {
      if ( ++arg < args.size() ) filter = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csheetbend\n\t-mtzin <filename>\n\t-colin-fo <colpath>\n\t-colin-free <colpath>\n\t-resolution <reso>\n\t-radius <radius>\n\t-radius-scale <scale>\n\t-pdbin <pdbin>\n\t-pdbout <pdbout>\n\t-free-flag <flag>\n\t-coord\n\t-u-iso\n\t-u-aniso\n\t-postrefine-coord\n\t-postrefine-u-iso\n\t-postrefine-u-aniso\n\t-cycles <cycles>\n\t-resolution-by-cycle <reso,reso,...>\n\t-aniso-obs\n\t-aniso-cal\n\t-b-iso-range <lo>,<high>\n\t-pseudo-regularize\n\t-refine-regularize-cycles <cycles>\n.\n";
    exit(1);
  }


  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::CCP4MAPfile xfile;
  clipper::MMDBfile mfile;
  clipper::MiniMol mmol;
  clipper::Xmap<float> map_in;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );

  // preliminary data opjects
  clipper::HKL_data<clipper::data32::F_sigF> fo0;
  clipper::HKL_data<clipper::data32::Flag>   free;
  clipper::HKL_data<clipper::data32::F_phi>  fc0;

  // read model
  if ( pdbfile != "NONE" ) {
    const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
    mfile.SetFlag( mmdbflags );
    mfile.read_file( pdbfile );
    mfile.import_minimol( mmol );
  }

  // read map
  if ( mapfile != "NONE" ) {
    xfile.open_read( mapfile );
    xfile.import_xmap( map_in );
    xfile.close_read();
  }

  // read X-ray data
  mtzin.open_read( ipfile );
  clipper::Spacegroup spgr0 = mtzin.spacegroup();
  clipper::Cell       cell = mtzin.cell();
  clipper::Resolution reso = mtzin.resolution();
  mtzin.import_hkl_data( fo0, ipcolfo );
  if ( ipcolfree != "NONE" ) mtzin.import_hkl_data( free, ipcolfree );
  if ( ipcolfc   != "NONE" ) mtzin.import_hkl_data( fc0,  ipcolfc   );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // defaults
  if ( !refxyz && !refuiso && !refuaniso) refxyz = true;
  if ( res <= 0.0 ) res = reso.limit();
  if ( resbycyc.size() == 0 ) resbycyc.push_back( res );
  if ( free.is_null() ) { free.init( fo0 ); }
  if ( fc0.is_null()  ) { fc0.init( fo0 );  }
  if ( refuiso || postrefuiso ) std::cout << std::endl << "B factor bounds " << clipper::Util::u2b(u_lo_abs) << " < B < " << clipper::Util::u2b(u_hi_abs) << std::endl;

  // results
  std::vector<results_by_cycle> results;


  // CASE 1:
  // Refine model against X-ray
  if ( pdbfile != "NONE" ) {
    std::cout << std::endl << "Refine MODEL against X-ray observations" << std::endl;

    // copy input molecule
    clipper::MiniMol mmol_orig = mmol;

    // loop over refine-regularize cycles
    for ( int cycrr = 0; cycrr < ncycrr; cycrr++ ) {
      std::cout << std::endl << std::endl << "Refine-regularize cycle: " << cycrr+1 << std::endl;

      // loop over cycles
      for ( int cyc = 0; cyc < ncyc; cyc++ ) {
        // check for final cycle
        bool lastcyc = ( cyc == ncyc-1 );

        // set resolution
        double fcyc = double(cyc) / std::max(double(ncyc-1),1.0);
        double fres = fcyc * double(resbycyc.size()-1);
        int ires0 = int( fres );
        int ires1 = std::min( ires0+1, int(resbycyc.size()-1) );
        double dres = fres - double(ires0);
        double rcyc = resbycyc[ires0] + dres*(resbycyc[ires1]-resbycyc[ires0]);

        // set radius if not user specified
        double radcyc = rad;
        if ( radcyc <= 0.0 ) radcyc = radscl * rcyc;
        std::cout << std::endl << "Cycle: " << cyc+1 << "  Resolution: " << rcyc << "  Radius: " << radcyc << std::endl;

        // truncate resolution
        clipper::Resolution rescyc( rcyc );
        clipper::HKL_sampling hklsam( cell, rescyc );
        clipper::HKL_data<clipper::data32::F_sigF> fo_all( spgr0, cell, hklsam ), fo( spgr0, cell, hklsam );
        for ( HRI ih = fo.first(); !ih.last(); ih.next() ) fo_all[ih] = fo0[ih.hkl()];
        for ( HRI ih = fo.first(); !ih.last(); ih.next() ) if (free[ih.hkl()].flag() != freeflag) fo[ih] = fo0[ih.hkl()];

        // calculate structure factors
        clipper::Atom_list atoms = mmol.atom_list();
        clipper::HKL_data<clipper::data32::F_phi> fc( fo );
        clipper::SFcalc_obs_bulk<float> sfcb;
        sfcb( fc, fo, atoms );
        //std::cerr << "Done sfcalc" << std::endl;

        // now calc R and R-free
        std::vector<double> params( 2, 1.0 );
        clipper::BasisFn_log_gaussian basisfn; // we're using a Gaussian because using a spline can result in negative square roots {see below}
        clipper::TargetFn_scaleLogF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn( fc, fo );

        clipper::ResolutionFn rfn( fo.hkl_info(), basisfn, targetfn, params);
        double r1w, f1w, r1f, f1f, Fo, Fc;
        r1w = f1w = r1f = f1f = 0.0;
        for ( HRI ih = fo_all.first(); !ih.last(); ih.next() )
          if ( !fo_all[ih].missing() ) {
            float sfac = exp(0.5*basisfn.f(ih.hkl(),cell,rfn.params()));
            Fo = fo_all[ih].f();
            Fc = sfac * fc[ih].f();

            if ( free[ih].flag() == freeflag ) {
                    r1f += fabs( Fo - Fc );
                    f1f += Fo;
            } else {
                    r1w += fabs( Fo - Fc );
                    f1w += Fo;
            }
          }

        r1w /= clipper::Util::max( f1w, 0.1 );
        r1f /= clipper::Util::max( f1f, 0.1 );

        std::cout << "R-factor      : " << r1w << std::endl
                  << "Free R-factor : " << r1f << " (at current resolution)" << std::endl;
      
        // do anisotropic scaling
        if ( aniso != NONE )  {
          clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
          clipper::SFscale_aniso<float> sfscl;
          if ( aniso == FOBS ) sfscl( fo, fc );  // scale Fobs
          if ( aniso == FCAL ) sfscl( fc, fo );  // scale Fcal
          std::cout << "\nAnisotropic scaling:\n"
                    << sfscl.u_aniso_orth(F).format() << std::endl;
        }

        // now do sigmaa calc
        clipper::HKL_data<clipper::data32::F_phi>   fb( fo ), fd( fo );
        clipper::HKL_data<clipper::data32::Phi_fom> phiw( fo );
        clipper::HKL_data<clipper::data32::Flag>    flag( fo );
        for ( HRI ih = flag.first(); !ih.last(); ih.next() )
          flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;

        // do sigmaa calc
        clipper::SFweight_spline<float> sfw( n_refln, n_param );
        sfw( fb, fd, phiw, fo, fc, flag );
        //std::cerr << "Done sigmaa" << std::endl;

        // expand to P1
        clipper::Spacegroup spgr1 = clipper::Spacegroup::p1();
        clipper::HKL_data<clipper::data32::F_phi> fphi( spgr1, cell, hklsam );
        clipper::HKL_data<clipper::data32::F_phi> dphi( fphi );
        for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
          fphi[ih] = fc[ih.hkl()];
          dphi[ih] = fd[ih.hkl()];
        }
        std::cout << "Reflections: " << fo0.base_hkl_info().num_reflections() << " P1: " << fo.base_hkl_info().num_reflections() << std::endl;

        //for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << fphi[ih].f() << " " << fphi[ih].phi() << std::endl;
        //for ( HRI ih = dphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << dphi[ih].f() << " " << dphi[ih].phi() << std::endl;

        // apply U value
        //fphi.compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(1.0,-uvalue) );

        // make grid if necessary
        clipper::Grid_sampling grid( spgr1, cell, rescyc );

        // make xmaps
        clipper::Xmap<float> cmap( spgr1, cell, grid );
        clipper::Xmap<float> dmap( spgr1, cell, grid );
        clipper::Xmap<float> mmap( spgr1, cell, grid );
        clipper::Xmap<float> x1map( spgr1, cell, grid );
        clipper::Xmap<float> x2map( spgr1, cell, grid );
        clipper::Xmap<float> x3map( spgr1, cell, grid );
        clipper::Xmap<float> x4map( spgr1, cell, grid );
        clipper::Xmap<float> x5map( spgr1, cell, grid );
        clipper::Xmap<float> x6map( spgr1, cell, grid );

        cmap.fft_from( fphi );
        dmap.fft_from( dphi );
        mmap = 1;

        // read pdb mask file
        if ( pdbmask != "NONE" ) {
          mmap = 0;
          clipper::MMDBfile mmfile;
          clipper::MiniMol mmolmsk;
          mmfile.read_file( pdbmask );
          mmfile.import_minimol( mmolmsk );
          clipper::EDcalc_mask<float> maskcalc( 2.5 );
          maskcalc( mmap, mmolmsk.atom_list() );
          clipper::Map_stats m( mmap );
          std::cout << "MASK " << m.min() << " " << m.mean() << " " << m.max() << std::endl;
        }

        // xyz refinement
        if ( refxyz || ( lastcyc && postrefxyz ) ) {
          std::cout << "REFINE XYZ" << std::endl;

          // make shift field
          if ( inclconst ) Shift_field_refine::shift_field_coord_const( cmap, dmap, mmap, x1map, x2map, x3map, radcyc, filter );
          else             Shift_field_refine::shift_field_coord      ( cmap, dmap, mmap, x1map, x2map, x3map, radcyc, filter );

          // read pdb and update
          for ( int p = 0; p < mmol.size(); p++ )
            for ( int m = 0; m < mmol[p].size(); m++ )
              for ( int a = 0; a < mmol[p][m].size(); a++ ) {
                const clipper::Coord_frac cf = mmol[p][m][a].coord_orth().coord_frac(cell);
                const float du = 2.0*x1map.interp<clipper::Interp_cubic>( cf );
                const float dv = 2.0*x2map.interp<clipper::Interp_cubic>( cf );
                const float dw = 2.0*x3map.interp<clipper::Interp_cubic>( cf );
                clipper::Coord_orth dxyz = clipper::Coord_frac(du,dv,dw).coord_orth(cell);
                mmol[p][m][a].set_coord_orth( mmol[p][m][a].coord_orth() + dxyz );
              }
        }

        // isotropic u refinement
        if ( refuiso || ( lastcyc && postrefuiso ) ) {
          std::cout << "REFINE U ISO" << std::endl;

          // make shift field
          if ( inclconst ) Shift_field_refine::shift_field_u_iso_const( cmap, dmap, mmap, x1map, radcyc, filter);
          else             Shift_field_refine::shift_field_u_iso      ( cmap, dmap, mmap, x1map, radcyc, filter);

          // read pdb and update
          for ( int p = 0; p < mmol.size(); p++ )
            for ( int m = 0; m < mmol[p].size(); m++ )
              for ( int a = 0; a < mmol[p][m].size(); a++ ) {
                const clipper::Coord_frac cf = mmol[p][m][a].coord_orth().coord_frac(cell);
                const float du = 1.0*x1map.interp<clipper::Interp_cubic>( cf );
                mmol[p][m][a].set_u_iso( clipper::Util::bound( u_lo_abs, mmol[p][m][a].u_iso() - du, u_hi_abs ) );
              }
        }

        // anisotropic u refinement
        if ( refuaniso || ( lastcyc && postrefuaniso ) ) {
          std::cout << "REFINE U ANISO" << std::endl;

          // make shift field
          if ( inclconst ) Shift_field_refine::shift_field_u_aniso_const( cmap, dmap, mmap, x1map, x2map, x3map, x4map, x5map, x6map, radcyc, rescyc, filter );
          else             Shift_field_refine::shift_field_u_aniso      ( cmap, dmap, mmap, x1map, x2map, x3map, x4map, x5map, x6map, radcyc, rescyc, filter );

          // read pdb and update
          for ( int p = 0; p < mmol.size(); p++ )
            for ( int m = 0; m < mmol[p].size(); m++ )
              for ( int a = 0; a < mmol[p][m].size(); a++ ) {
                const clipper::Coord_frac cf = mmol[p][m][a].coord_orth().coord_frac(cell);
                if ( mmol[p][m][a].u_aniso_orth().is_null() == true ) { // if necessary, create new anisotropic entries for the pdb if none exist
                       const float db1 = mmol[p][m][a].u_iso();
                       const float db2 = db1;
                       const float db3 = db1;
                       clipper::U_aniso_orth db(db1,db2,db3,0,0,0);
                       mmol[p][m][a].set_u_aniso_orth( db );
                     }
                   float dx1 = -1.0*x1map.interp<clipper::Interp_cubic>( cf );
                   float dx2 = -1.0*x2map.interp<clipper::Interp_cubic>( cf );
                   float dx3 = -1.0*x3map.interp<clipper::Interp_cubic>( cf );
                   float dx4 = -1.0*x4map.interp<clipper::Interp_cubic>( cf );
                   float dx5 = -1.0*x5map.interp<clipper::Interp_cubic>( cf );
                   float dx6 = -1.0*x6map.interp<clipper::Interp_cubic>( cf );
                   clipper::U_aniso_orth dx(dx1,dx2,dx3,dx4,dx5,dx6);
                   mmol[p][m][a].set_u_aniso_orth( mmol[p][m][a].u_aniso_orth() + dx );
                   mmol[p][m][a].set_u_iso( mmol[p][m][a].u_aniso_orth().u_iso() );

                   if (mmol[p][m][a].u_iso() != mmol[p][m][a].u_iso()) { // if this produces a NaN value, reverse the shift (equivalent to not applying it)

                   dx = -dx;

                   mmol[p][m][a].set_u_aniso_orth( mmol[p][m][a].u_aniso_orth() + dx );
                   mmol[p][m][a].set_u_iso( mmol[p][m][a].u_aniso_orth().u_iso() );}

              }
        }

        // store results for XML
        results_by_cycle result;
        result.cyclerr = cycrr;
        result.cycle = cyc;
        result.resolution = rcyc;
        result.radius = radcyc;
        result.rwork = r1w;
        result.rfree = r1f;
        result.nrefl = fo0.base_hkl_info().num_reflections();
        result.nreflp1 = fo.base_hkl_info().num_reflections();
        results.push_back( result );

      } // end of cycle loop

      // pseudo regularization
      if ( pseudoreg ) {
        std::cout << "PSEUDOREGULARIZE" << std::endl;
        PseudoRegularize::regularize( mmol, mmol_orig );
      }

    } // end of cycle loop

  } // end of case 1

  // CASE 2:
  // refine map against X-ray
  if ( pdbfile == "NONE" ) {
    std::cout << std::endl << "Refine MAP against X-ray observations" << std::endl;

    // 1. Don't read model, read map
    // 2. Before first cycle, make shift maps
    // 3. In first cycle, make shifted map using shifts
    // 4. Calc structure factors from this map
    // 5. At end of cycle, update shift maps
    if ( map_in.is_null() ) {
      map_in.init( spgr0, cell, clipper::Grid_sampling( spgr0, cell, reso ) );
      map_in.fft_from( fc0 );
    }

    typedef clipper::Xmap_base::Map_reference_index MRI;
    clipper::Spacegroup    spgrm = map_in.spacegroup();
    clipper::Cell          cellm = map_in.cell();
    clipper::Grid_sampling gridm = map_in.grid_sampling();

    clipper::Xmap<float> map_shifted( spgrm, cellm, gridm );
    clipper::Xmap<float> dumap( spgrm, cellm, gridm );
    clipper::Xmap<float> dvmap( spgrm, cellm, gridm );
    clipper::Xmap<float> dwmap( spgrm, cellm, gridm );
    dumap = dvmap = dwmap = 0.0;

    // loop over cycles
    for ( int cyc = 0; cyc < ncyc; cyc++ ) {

      // set resolution
      double fcyc = double(cyc) / std::max(double(ncyc-1),1.0);
      double fres = fcyc * double(resbycyc.size()-1);
      int ires0 = int( fres );
      int ires1 = std::min( ires0+1, int(resbycyc.size()-1) );
      double dres = fres - double(ires0);
      double rcyc = resbycyc[ires0] + dres*(resbycyc[ires1]-resbycyc[ires0]);

      // set radius if not user specified
      double radcyc = rad;
      if ( radcyc <= 0.0 ) radcyc = radscl * rcyc;
      std::cout << std::endl << "Cycle: " << cyc+1 << "  Resolution: " << rcyc << "  Radius: " << radcyc << std::endl;

      // truncate resolution
      clipper::Resolution rescyc( rcyc );
      clipper::HKL_sampling hklsam( cell, rescyc );
      clipper::HKL_data<clipper::data32::F_sigF> fo_all( spgr0, cell, hklsam ), fo( spgr0, cell, hklsam );
      for ( HRI ih = fo.first(); !ih.last(); ih.next() ) fo_all[ih] = fo0[ih.hkl()];
      for ( HRI ih = fo.first(); !ih.last(); ih.next() ) if (free[ih.hkl()].flag() != freeflag) fo[ih] = fo0[ih.hkl()];

      // calculate structure factors
      for ( MRI ix = map_shifted.first(); !ix.last(); ix.next() ) {
        clipper::Coord_frac cf = ix.coord().coord_frac(gridm);
        cf -= clipper::Coord_frac( dumap[ix], dvmap[ix], dwmap[ix] );
        map_shifted[ix] = map_in.interp<clipper::Interp_cubic>( cf );
      }
      clipper::HKL_data<clipper::data32::F_phi> fc( fo );
      map_shifted.fft_to( fc );
      //std::cerr << "Done sfcalc" << std::endl;

      // now calc R and R-free
      std::vector<double> params( 2, 1.0 );
      clipper::BasisFn_log_gaussian basisfn; // we're using a Gaussian because using a spline can result in negative square roots {see below}
      clipper::TargetFn_scaleLogF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> targetfn( fc, fo );

      clipper::ResolutionFn rfn( fo.hkl_info(), basisfn, targetfn, params);
      double r1w, f1w, r1f, f1f, Fo, Fc;
      r1w = f1w = r1f = f1f = 0.0;
      for ( HRI ih = fo_all.first(); !ih.last(); ih.next() )
        if ( !fo_all[ih].missing() ) {
          float sfac = exp(0.5*basisfn.f(ih.hkl(),cell,rfn.params()));
          Fo = fo_all[ih].f();
          Fc = sfac * fc[ih].f();

          if ( free[ih].flag() == freeflag ) {
                  r1f += fabs( Fo - Fc );
                  f1f += Fo;
          } else {
                  r1w += fabs( Fo - Fc );
                  f1w += Fo;
          }
        }

      r1w /= clipper::Util::max( f1w, 0.1 );
      r1f /= clipper::Util::max( f1f, 0.1 );

      std::cout << "R-factor      : " << r1w << std::endl
                << "Free R-factor : " << r1f << " (at current resolution)" << std::endl;
    
      // do anisotropic scaling
      if ( aniso != NONE )  {
        clipper::SFscale_aniso<float>::TYPE F = clipper::SFscale_aniso<float>::F;
        clipper::SFscale_aniso<float> sfscl;
        if ( aniso == FOBS ) sfscl( fo, fc );  // scale Fobs
        if ( aniso == FCAL ) sfscl( fc, fo );  // scale Fcal
        std::cout << "\nAnisotropic scaling:\n"
                  << sfscl.u_aniso_orth(F).format() << std::endl;
      }

      // now do sigmaa calc
      clipper::HKL_data<clipper::data32::F_phi>   fb( fo ), fd( fo );
      clipper::HKL_data<clipper::data32::Phi_fom> phiw( fo );
      clipper::HKL_data<clipper::data32::Flag>    flag( fo );
      for ( HRI ih = flag.first(); !ih.last(); ih.next() )
        flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;

      // do sigmaa calc
      clipper::SFweight_spline<float> sfw( n_refln, n_param );
      sfw( fb, fd, phiw, fo, fc, flag );
      //std::cerr << "Done sigmaa" << std::endl;

      // expand to P1
      clipper::Spacegroup spgr1 = clipper::Spacegroup::p1();
      clipper::HKL_data<clipper::data32::F_phi> fphi( spgr1, cell, hklsam );
      clipper::HKL_data<clipper::data32::F_phi> dphi( fphi );
      for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
        fphi[ih] = fc[ih.hkl()];
        dphi[ih] = fd[ih.hkl()];
      }
      std::cout << "Reflections: " << fo0.base_hkl_info().num_reflections() << " P1: " << fo.base_hkl_info().num_reflections() << std::endl;

      //for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << fphi[ih].f() << " " << fphi[ih].phi() << std::endl;
      //for ( HRI ih = dphi.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << dphi[ih].f() << " " << dphi[ih].phi() << std::endl;

      // apply U value
      //fphi.compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(1.0,-uvalue) );

      // make grid if necessary
      clipper::Grid_sampling grid( spgr1, cell, rescyc );

      // make xmaps
      clipper::Xmap<float> cmap( spgr1, cell, grid );
      clipper::Xmap<float> dmap( spgr1, cell, grid );
      clipper::Xmap<float> mmap( spgr1, cell, grid );
      clipper::Xmap<float> x1map( spgr1, cell, grid );
      clipper::Xmap<float> x2map( spgr1, cell, grid );
      clipper::Xmap<float> x3map( spgr1, cell, grid );
      clipper::Xmap<float> x4map( spgr1, cell, grid );
      clipper::Xmap<float> x5map( spgr1, cell, grid );
      clipper::Xmap<float> x6map( spgr1, cell, grid );

      cmap.fft_from( fphi );
      dmap.fft_from( dphi );
      mmap = 1;

      std::cout << "REFINE XYZ" << std::endl;

      // make shift field
      Shift_field_refine::shift_field_coord( cmap, dmap, mmap, x1map, x2map, x3map, radcyc, filter );

      // read pdb and update
      for ( MRI ix = dumap.first(); !ix.last(); ix.next() ) {
        const clipper::Coord_frac cf = ix.coord().coord_frac(gridm);
        // do we need to add the current shifts to cf here?
        dumap[ix] += 2.0*x1map.interp<clipper::Interp_cubic>( cf );
        dvmap[ix] += 2.0*x2map.interp<clipper::Interp_cubic>( cf );
        dwmap[ix] += 2.0*x3map.interp<clipper::Interp_cubic>( cf );
      }

      // store results for XML
      results_by_cycle result;
      result.cyclerr = 0;
      result.cycle = cyc;
      result.resolution = rcyc;
      result.radius = radcyc;
      result.rwork = r1w;
      result.rfree = r1f;
      result.nrefl = fo0.base_hkl_info().num_reflections();
      result.nreflp1 = fo.base_hkl_info().num_reflections();
      results.push_back( result );

    } // end of cycle loop

    // calculate final map
    for ( MRI ix = map_shifted.first(); !ix.last(); ix.next() ) {
      clipper::Coord_frac cf = ix.coord().coord_frac(gridm);
      cf -= clipper::Coord_frac( dumap[ix], dvmap[ix], dwmap[ix] );
      map_shifted[ix] = map_in.interp<clipper::Interp_cubic>( cf );
    }
    map_shifted.fft_to( fc0 );
    map_in = map_shifted;

  } // end of case 2


  // write output files
  if ( pdbfile != "NONE" ) {
    mfile.export_minimol( mmol );
    mfile.write_file( pdboutfile );
  }
  if ( mapoutfile != "NONE" ) {
    xfile.open_write( mapoutfile );
    xfile.import_xmap( map_in );
    xfile.close_write();
  }
  if ( mtzoutfile != "NONE" ) {
    mtzout.open_append( ipfile, mtzoutfile );
    mtzout.export_hkl_data( fc0, opcol );
    mtzout.close_append();
  }

  // write xml
  if ( opxml != "NONE" ) {
    FILE *f = fopen( opxml.c_str(), "w" );
    if ( f == NULL ) clipper::Message::message( clipper::Message_fatal( "Error: Could not open xml temporary file: "+opxml ) );
    fprintf( f, "<SheetbendResult>\n" );
    fprintf( f, " <Title>%s</Title>\n", title.c_str() );
    fprintf( f, " <Cycles>\n" );
    for ( int c = 0; c < results.size(); c++ ) {
      fprintf( f, "  <Cycle>\n" );
      fprintf( f, "   <RegularizeNumber>%i</RegularizeNumber>\n", results[c].cyclerr+1 );
      fprintf( f, "   <Number>%i</Number>\n", results[c].cycle+1 );
      fprintf( f, "   <Resolution>%f</Resolution>\n", results[c].resolution );
      fprintf( f, "   <Radius>%f</Radius>\n",         results[c].radius );
      fprintf( f, "   <Rwork>%f</Rwork>\n",           results[c].rwork );
      fprintf( f, "   <Rfree>%f</Rfree>\n",           results[c].rfree );
      fprintf( f, "   <NumberOfReflections>%i</NumberOfReflections>\n",     results[c].nrefl   );
      fprintf( f, "   <NumberOfReflectionsP1>%i</NumberOfReflectionsP1>\n", results[c].nreflp1 );
      fprintf( f, "  </Cycle>\n" );
    }
    fprintf( f, " </Cycles>\n" );
    fprintf( f, " <Final>\n" );
    int c = results.size()-1;
    fprintf( f, "   <RegularizeNumber>%i</RegularizeNumber>\n", results[c].cyclerr+1 );
    fprintf( f, "   <Number>%i</Number>\n", results[c].cycle+1 );
    fprintf( f, "   <Resolution>%f</Resolution>\n", results[c].resolution );
    fprintf( f, "   <Radius>%f</Radius>\n",         results[c].radius );
    fprintf( f, "   <Rwork>%f</Rwork>\n",           results[c].rwork );
    fprintf( f, "   <Rfree>%f</Rfree>\n",           results[c].rfree );
    fprintf( f, "   <NumberOfReflections>%i</NumberOfReflections>\n",     results[c].nrefl   );
    fprintf( f, "   <NumberOfReflectionsP1>%i</NumberOfReflectionsP1>\n", results[c].nreflp1 );
    fprintf( f, " </Final>\n" );
    fprintf( f, "</SheetbendResult>\n" );
    fclose(f);
  }

  prog.set_termination_message( "Normal termination" );
}

