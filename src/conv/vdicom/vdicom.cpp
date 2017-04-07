//
// vdicom.cc
// created as d2v.cc at Thu Jul 17 16:49:23 2014
// in /home/karda/src/i2o/
// by \bk: bernd.kardatzki@med.uni-tuebingen.de
//
// purpose: create "vista" files (.v) from dicom files
//
// status: WIP
//

#include <iostream>
#include <string>
#include <list>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <cmath>

#include <sys/stat.h>

#include <viaio/Vlib.h>
#include <viaio/VImage.h>
#include <viaio/mu.h>
#include <viaio/option.h>

#include "i2o_dcm.h"

#define nDEBUG

using namespace std;


void UsageAndDie( char* argv0, string issue )
{
  cerr << endl
       << "Create one or more \"vista\" files (.v) from dicom files." << endl
       << "Input may be a directory name (no recursion) or a DICOMDIR file or a list of filenames." << endl
       << endl
       << " Usage:"
       << endl << endl
       << " \""
       << basename( argv0 ) << " [-in directory|dicomdir|list of files] [-v] [-f] [-t directory] [-out prefix]\"\n\n"
    "-in : Inputfiles(s),\n"
    "      list of individual dicom files, directory containing the dicom files or a single DICOMDIR file.\n"
    "-out: Outputfilename,\n"
    "      how to name the resulting files (e.g. \"-out out\" creates \"out.v\" in the most simple case).\n"
    "      Otherwise the series number and the series description will be used for filename creation (.v appended).\n"
    "-t  : Target directory,\n"
    "      where to put the resulting \".v\" file(s) (will be created if necessary, default is \".\").\n"
    "-f  : Force,\n"
    "      overwrite a preexisting file instead of appending a number to the newly created file(s).\n"
    "-v  : Verbose,\n"
    "      if omitted you will see only errors and warnings and \"done.\" at the end.\n"
    // "-r : Readonly,\n"
    // "     just pretend, don't write anything.\n"
    ;

  if( issue.size() )
    cerr << "\n*** " << issue << endl;

  exit(1);
}


bool  // free "s" from pollution, e.g. for filenames use "_" instead of special characters
purify( string& s )
{
  bool purified = false;

  for( unsigned i = 0; i < s.length(); i++ )
    if( ( s[i] > 0x7A && s[i] < 0x7F ) // leave umlauts untouched
#if 0	// more explicitely:
	|| s[i] == ' ' 	|| s[i] == '/' 	|| s[i] == '\\' || s[i] == '\''	|| s[i] == '('
	|| s[i] == ')'	|| s[i] == '`'	|| s[i] == '~'	|| s[i] == '*'	|| s[i] == '>'
	|| s[i] == '<'	|| s[i] == ':'	|| s[i] == '&'	|| s[i] == '$'	|| s[i] == '"'
	|| s[i] == '#'	|| s[i] == '!'	|| s[i] == '§'	|| s[i] == '['	|| s[i] == ']'
	|| s[i] == '{'	|| s[i] == '}'	|| s[i] == '?'	|| s[i] == '´'
#else
	// see e.g. http://en.wikipedia.org/wiki/ASCII for an ASCII table
	|| (s[i] < 0x30 && s[i] != '-' && s[i] != '.' )
	|| (s[i] > 0x39 && s[i] < 0x41 )
	|| (s[i] > 0x5A && s[i] < 0x61 && s[i] != '_' )
#endif
	)
      {
	s[i] = '_';
	purified = true;
      }
  return purified;
}

bool // true if directory "d" already exists or if creation was successfull
mkdir_p( string d )
{
  struct stat sb;

  if( stat(d.c_str(), &sb) == -1 )
    {
      if( mkdir( d.c_str(), 0777 ) == 0 )
	return true; // that was easy

      if( errno != ENOENT )
	{ // pointless to continue
	  return false;
	}

      size_t pos = 0;
      while( (pos = d.find_first_of('/', pos )) != string::npos )
	{
	  string t = d.substr( 0, ++pos );
	  if( ! mkdir_p( t ) )
	    return false;
	}
      if( mkdir( d.c_str(), 0777 ) == 0 )
	return true;
    }
  else
    if( S_ISDIR( sb.st_mode ) )
      return true;

  return false;
}


// little housekeeping, remember which files we just created
list<string> MyOwnFiles;


FILE*
CreateOutputFileOrDie ( string tgtDir,
			string tgtPrefix,
			i2o_dcm* dcm, // to get an idea for the filename
			bool overwrite,
			string& name )
{
  FILE*  vOut      = nullptr;
  string vOutName0 = tgtDir + "/";

  if( ! mkdir_p( tgtDir ) )
    {
      cerr << "target directory problem, could not create \"" << tgtDir << "\"." << endl;
      exit(1);
    }

  if( tgtPrefix.size() )
    vOutName0.append( tgtPrefix );
  else
    {
      string str;

      if( dcm->SeriesNr(str) )
	vOutName0.append( str );
      else
	vOutName0.append( "no_series_nr" );
      vOutName0.append( "_" );

      if( dcm->SeriesDescription(str) )
	vOutName0.append( str );
      else
	if( dcm->ProtocolName(str) )
	  vOutName0.append( str );
	else
	  vOutName0.append( "no_name" );
    }

  struct stat sb;
  string vOutName = vOutName0 + ".v";
  int    n        = 1;

  while( stat( vOutName.c_str(), &sb ) == 0 )
    { // oopsi, does already exist, find a solution ...
      if( overwrite )
	{ // fine but don't shoot ourselve
	  bool ownFile = false;

	  for( auto &f : MyOwnFiles )
	    if( vOutName == f )
	      {
		ownFile = true;
		break; // leave the for loop, create a new filename and try again
	      }
	  if( ! ownFile )
	    break; // not our business, leave while loop and overwrite
	}
      vOutName = vOutName0 + "-" + to_string(n++) + ".v";
    }

  vOut = VOpenOutputFile ( vOutName.c_str(), TRUE);

  if( vOut )
    MyOwnFiles.push_back( vOutName );
  else
    {
      cerr << "VOpenOutputFile(\"" << vOutName << "\") failed: " << strerror(errno) << endl;
      exit(1);
    }

  name = vOutName;

  return vOut;
}

////////////////////////////////////////////
// shamelessly stolen from nifti_io.c :  //
// (slightly modified)                  //
/////////////////////////////////////////
typedef struct { float m00,m01,m02,m03,m10,m11,m12,m13,m20,m21,m22,m23,m30,m31,m32,m33; } mat44 ;
typedef struct { float m00,m01,m02,m10,m11,m12,m20,m21,m22,m230,m31,m32; } mat33 ;

float mat33_determ( mat33& R )   /* determinant of 3x3 matrix */
{
  return ( R.m00 * (R.m11*R.m22 - R.m21*R.m12) -
	   R.m10 * (R.m01*R.m22 + R.m21*R.m02) +
	   R.m20 * (R.m01*R.m12 - R.m11*R.m02) ) ;
}

float mat33_rownorm( mat33& A )  /* max row norm of 3x3 matrix */
{
  float r1 = abs(A.m00) + abs(A.m01) + abs(A.m02);
  float r2 = abs(A.m10) + abs(A.m11) + abs(A.m12);
  float r3 = abs(A.m20) + abs(A.m21) + abs(A.m22);
  if( r1 < r2 )
    r1 = r2 ;
  if( r1 < r3 )
    r1 = r3 ;
  return r1 ;
}

float mat33_colnorm( mat33& A )  /* max column norm of 3x3 matrix */
{
  float r1 = abs(A.m00) + abs(A.m10) + abs(A.m20);
  float r2 = abs(A.m01) + abs(A.m11) + abs(A.m21);
  float r3 = abs(A.m02) + abs(A.m12) + abs(A.m22);
  if( r1 < r2 )
    r1 = r2 ;
  if( r1 < r3 )
    r1 = r3 ;
  return r1 ;
}

mat33 mat33_inverse( mat33& R )   /* inverse of 3x3 matrix */
{
  double r11,r12,r13,r21,r22,r23,r31,r32,r33;
  double det ;
  mat33 Q ;
                                           /*  INPUT MATRIX:  */
  r11 = R.m00; r12 = R.m01; r13 = R.m02;  /* [ r11 r12 r13 ] */
  r21 = R.m10; r22 = R.m11; r23 = R.m12;  /* [ r21 r22 r23 ] */
  r31 = R.m20; r32 = R.m21; r33 = R.m22;  /* [ r31 r32 r33 ] */
  
  det = r11*r22*r33 - r11*r32*r23 - r21*r12*r33 + r21*r32*r13 + r31*r12*r23 - r31*r22*r13 ;
  
  if( det != 0.0l )
    det = 1.0l / det ;
  
  Q.m00 = (float)( det*( r22*r33-r32*r23) ) ;
  Q.m01 = (float)( det*(-r12*r33+r32*r13) ) ;
  Q.m02 = (float)( det*( r12*r23-r22*r13) ) ;
  
  Q.m10 = (float)( det*(-r21*r33+r31*r23) ) ;
  Q.m11 = (float)( det*( r11*r33-r31*r13) ) ;
  Q.m12 = (float)( det*(-r11*r23+r21*r13) ) ;
  
  Q.m20 = (float)( det*( r21*r32-r31*r22) ) ;
  Q.m21 = (float)( det*(-r11*r32+r31*r12) ) ;
  Q.m22 = (float)( det*( r11*r22-r21*r12) ) ;
  
  return Q ;
}

mat33 mat33_polar( mat33& A )
{
  mat33 X , Y , Z ;
  float alp,bet,gam,gmi , dif=1.0f ;
  int k=0 ;
  
  X = A ;
  
  /* force matrix to be nonsingular */

  gam = mat33_determ(X) ;
  while( gam == 0.0 )
    {        /* perturb matrix */
      gam = (float)( 0.00001 * ( 0.001 + mat33_rownorm(X) ) ) ;
      X.m00 += gam ; X.m11 += gam ; X.m22 += gam ;
      gam = mat33_determ(X) ;
    }
  
  while(1)
    {
      Y = mat33_inverse(X) ;
      if( dif > 0.3 ){     /* far from convergence */
	alp = (float)( sqrt( mat33_rownorm(X) * mat33_colnorm(X) ) ) ;
	bet = (float)( sqrt( mat33_rownorm(Y) * mat33_colnorm(Y) ) ) ;
	gam = (float)( sqrt( bet / alp ) ) ;
	gmi = (float)( 1.0 / gam ) ;
      }
      else
	{
	  gam = gmi = 1.0f ;  /* close to convergence */
	}
      Z.m00 = (float)( 0.5 * ( gam*X.m00 + gmi*Y.m00 ) ) ;
      Z.m01 = (float)( 0.5 * ( gam*X.m01 + gmi*Y.m10 ) ) ;
      Z.m02 = (float)( 0.5 * ( gam*X.m02 + gmi*Y.m20 ) ) ;
      Z.m10 = (float)( 0.5 * ( gam*X.m10 + gmi*Y.m01 ) ) ;
      Z.m11 = (float)( 0.5 * ( gam*X.m11 + gmi*Y.m11 ) ) ;
      Z.m12 = (float)( 0.5 * ( gam*X.m12 + gmi*Y.m21 ) ) ;
      Z.m20 = (float)( 0.5 * ( gam*X.m20 + gmi*Y.m02 ) ) ;
      Z.m21 = (float)( 0.5 * ( gam*X.m21 + gmi*Y.m12 ) ) ;
      Z.m22 = (float)( 0.5 * ( gam*X.m22 + gmi*Y.m22 ) ) ;
      
      dif = (float)( abs(Z.m00-X.m00)+abs(Z.m01-X.m01)
		     +abs(Z.m02-X.m02)+abs(Z.m10-X.m10)
		     +abs(Z.m11-X.m11)+abs(Z.m12-X.m12)
		     +abs(Z.m20-X.m20)+abs(Z.m21-X.m21)
		     +abs(Z.m22-X.m22) );
    
      k = k+1 ;
      if( k > 100 || dif < 3.e-6 )
	break ;  /* convergence or exhaustion */
      X = Z ;
    }
  
  return Z ;
}

void mat44_to_quatern( mat44& R,
		       float* qb, float* qc, float* qd,
		       float* qx, float* qy, float* qz,
		       float* dx, float* dy, float* dz, float* qfac )
{
  double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
  double xd,yd,zd , a,b,c,d ;
  mat33 P,Q ;
  
  /* offset outputs are read write out of input matrix  */
  if( qx )
    *qx = R.m03;
  if( qy )
    *qy = R.m13;
  if( qz )
    *qz = R.m23;

  /* load 3x3 matrix into local variables */
  
  r11 = R.m00 ; r12 = R.m01 ; r13 = R.m02 ;
  r21 = R.m10 ; r22 = R.m11 ; r23 = R.m12 ;
  r31 = R.m20 ; r32 = R.m21 ; r33 = R.m22 ;
  
  /* compute lengths of each column; these determine grid spacings  */
  
  xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
  yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
  zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;
  
  /* if a column length is zero, patch the trouble */
  
  if( xd == 0.0l )
    {
      r11 = 1.0l ;
      r21 = r31 = 0.0l ;
      xd = 1.0l ;
    }
  if( yd == 0.0l )
    {
      r22 = 1.0l ;
      r12 = r32 = 0.0l ;
      yd = 1.0l ;
    }
  if( zd == 0.0l )
    {
      r33 = 1.0l ;
      r13 = r23 = 0.0l ;
      zd = 1.0l ;
    }
  
  /* assign the output lengths */
  if( dx )
    *dx = (float)xd;
  if( dy )
    *dy = (float)yd;
  if( dz )
    *dz = (float)zd;
  
  /* normalize the columns */
  
  r11 /= xd ; r21 /= xd ; r31 /= xd ;
  r12 /= yd ; r22 /= yd ; r32 /= yd ;
  r13 /= zd ; r23 /= zd ; r33 /= zd ;
  
  /* At this point, the matrix has normal columns, but we have to allow
     for the fact that the hideous user may not have given us a matrix
     with orthogonal columns.
     
     So, now find the orthogonal matrix closest to the current matrix.
     
     One reason for using the polar decomposition to get this
     orthogonal matrix, rather than just directly orthogonalizing
     the columns, is so that inputting the inverse matrix to R
     will result in the inverse orthogonal matrix at this point.
     If we just orthogonalized the columns, this wouldn't necessarily hold. */
  
  Q.m00 = (float)r11 ; Q.m01 = (float)r12 ; Q.m02 = (float)r13 ; /* load Q */
  Q.m10 = (float)r21 ; Q.m11 = (float)r22 ; Q.m12 = (float)r23 ;
  Q.m20 = (float)r31 ; Q.m21 = (float)r32 ; Q.m22 = (float)r33 ;
  
  P = mat33_polar(Q) ;  /* P is orthog matrix closest to Q */
  
  r11 = P.m00 ; r12 = P.m01 ; r13 = P.m02 ; /* unload */
  r21 = P.m10 ; r22 = P.m11 ; r23 = P.m12 ;
  r31 = P.m20 ; r32 = P.m21 ; r33 = P.m22 ;
  
  /*                            [ r11 r12 r13 ]               */
  /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
  /*                            [ r31 r32 r33 ]               */
  
  /* compute the determinant to determine if it is proper */
  
  zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
    +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */
  
  if( zd > 0 )
    {             /* proper */
      if( qfac )
	*qfac = 1.0;
    }
  else
    {             /* improper ==> flip 3rd column */
      if( qfac )
	*qfac = -1.0;
      
      r13 = -r13 ;
      r23 = -r23 ;
      r33 = -r33 ;
    }
  
  /* now, compute quaternion parameters */
  
  a = r11 + r22 + r33 + 1.0l ;
  
  if( a > 0.5l )
    {                /* simplest case */
      a = 0.5l * sqrt(a) ;
      b = 0.25l * (r32-r23) / a ;
      c = 0.25l * (r13-r31) / a ;
      d = 0.25l * (r21-r12) / a ;
    }
  else
    {                       /* trickier case */
      xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
      yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
      zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
      if( xd > 1.0 )
	{
	  b = 0.5l * sqrt(xd) ;
	  c = 0.25l* (r12+r21) / b ;
	  d = 0.25l* (r13+r31) / b ;
	  a = 0.25l* (r32-r23) / b ;
	}
      else
	if( yd > 1.0 )
	  {
	    c = 0.5l * sqrt(yd) ;
	    b = 0.25l* (r12+r21) / c ;
	    d = 0.25l* (r23+r32) / c ;
	    a = 0.25l* (r13-r31) / c ;
	  }
	else
	  {
	    d = 0.5l * sqrt(zd) ;
	    b = 0.25l* (r13+r31) / d ;
	    c = 0.25l* (r23+r32) / d ;
	    a = 0.25l* (r21-r12) / d ;
	  }
      if( a < 0.0l )
	{
	  b = -b ;
	  c = -c ;
	  d = -d;
	  a = -a;
	}
  }
  if( qb )
    *qb = (float)b;
  if( qc )
    *qc = (float)c;
  if( qd )
    *qd = (float)d;

  return ;
}

// there's a shift, make the coordinates consistant
void FixMosaicImgPositionError( i2o_dcm* dcm, int slice, vector<double>& v3d )
{
  int    cols                   = dcm->Cols();
  int    rows                   = dcm->Rows();
  int    totalCols              = dcm->TotalCols();
  int    totalRows              = dcm->TotalRows();
  double voxHeight              = dcm->VoxelHeight();
  double voxWidth               = dcm->VoxelWidth();
  double spacing                = dcm->SliceSpacing();
  const vector<double>& vecRow  = dcm->ImgOrientRow();
  const vector<double>& vecCol  = dcm->ImgOrientCol();
  const vector<double>& vecNorm = dcm->ImgOrientNormal();
  const vector<double>& imgPos  = dcm->ImgPosition();

  double f_c = cols * voxWidth  * ( 0.5 * totalCols / cols  - 0.5 );
  double f_r = rows * voxHeight * ( 0.5 * totalRows / rows - 0.5 );

  v3d[0] = imgPos[0] + f_c * vecRow[0] + f_r * vecCol[0] + slice * spacing * vecNorm[0];
  v3d[1] = imgPos[1] + f_c * vecRow[1] + f_r * vecCol[1] + slice * spacing * vecNorm[1];
  v3d[2] = imgPos[2] + f_c * vecRow[2] + f_r * vecCol[2] + slice * spacing * vecNorm[2];
}


// try to "mimic" some part of a nifti header
VAttrList CreateGeoInfo( i2o_dcm* dcm,
			 int dimension,
			 int slices,
			 int timePoints )
{
  VAttrList vGeoInfo = VCreateAttrList();
  string    str;
  float     qfac;
  short     dimInfo;
  short     freq_dim,phase_dim,slice_dim;

  VSetAttr( vGeoInfo, "intent_name", nullptr, VStringRepn, (VString)"" );

  if( dcm->InPlanePhaseEncodingDir( str ) )
    {
      // nifti puts freq_dim, phase_dim, slice_dim all in one byte
      if( str == "ROW" )
	{
	  dimInfo  = 2 | (1 << 2) | (3 << 4);
	  freq_dim  = 2;
	  phase_dim = 1;
	  slice_dim = 3;
	}
      else  if( str == "COL" )
	{
	  dimInfo = 1 | (2 << 2) | (3 << 4);
	  freq_dim  = 1;
	  phase_dim = 2;
	  slice_dim = 3;
	}
      else
	{
	  dimInfo = (3 << 4);
	  freq_dim  = 0;
	  phase_dim = 0;
	  slice_dim = 3;
	}
      // leave it for reference comparison
      VSetAttr( vGeoInfo,"dim_info", nullptr, VShortRepn, dimInfo );
      // and put it out explicitly
      VSetAttr( vGeoInfo, "freq_dim", nullptr, VShortRepn, freq_dim );
      VSetAttr( vGeoInfo, "phase_dim",nullptr, VShortRepn, phase_dim );
      VSetAttr( vGeoInfo, "slice_dim",nullptr, VShortRepn, slice_dim );

#ifdef DEBUG
      cerr << "dbg: dim_info = " << dimInfo << endl;
#endif
    }

  float* dim = (float*) VMalloc( 8 * sizeof(float) );

  dim[0] = dimension;
  dim[1] = dcm->Cols();
  dim[2] = dcm->Rows();
  dim[3] = slices;
  dim[4] = timePoints;
  dim[5] = 1;
  dim[6] = 4; /*DT_SIGNED_SHORT*/ /* correct would be: DT_UINT_16*/;
  dim[7] = 1;
#ifdef DEBUG
  cerr << "dbg: dim = " << dim[0] << "," << dim[1] << "," << dim[2] << "," << dim[3] << "," << dim[4] << endl;
#endif
  VAttrList eList = VCreateAttrList();
  VBundle eBundle = VCreateBundle ("bundle",eList,8*sizeof(float),(VPointer)dim);
  VSetAttr(vGeoInfo,"dim",nullptr,VBundleRepn,eBundle);

  const vector<double>& imgOrientRow    = dcm->ImgOrientRow();
  const vector<double>& imgOrientCol    = dcm->ImgOrientCol();
  const vector<double>& imgOrientNormal = dcm->ImgOrientNormal();
  vector<double>& imgPosition           = dcm->ImgPosition();

  mat44 R0;

  if( imgOrientRow.size() == 3 )
    {
      R0.m00 = -imgOrientRow[0];
      R0.m10 = -imgOrientRow[1];
      R0.m20 =  imgOrientRow[2];
      R0.m30 =  0;

      if( imgOrientCol.size() == 3 )
	{
	  R0.m01 = -imgOrientCol[0];
	  R0.m11 = -imgOrientCol[1];
	  R0.m21 =  imgOrientCol[2];
	  R0.m31 =  0;

	  if( imgOrientNormal.size() == 3 )
	    {
	      R0.m02 = -imgOrientNormal[0];
	      R0.m12 = -imgOrientNormal[1];
	      R0.m22 =  imgOrientNormal[2];
	      R0.m32 = 0;

	      if( imgPosition.size() == 3 )
		{
		  if( dcm->Mosaic() )
		    FixMosaicImgPositionError( dcm, 0, imgPosition );

		  R0.m03 = -imgPosition[0];
		  R0.m13 = -imgPosition[1];
		  R0.m23 =  imgPosition[2];
		  R0.m33 =  1;

		  float qb,qc,qd,qx,qy,qz,dx,dy,dz;

		  mat44_to_quatern( R0, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac );
		  //nifti_mat44_to_quatern( R0, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac );

		  float* quat = (float*) VMalloc( 6 * sizeof(float) );

		  quat[0] = qb; quat[1] = qc; quat[2] = qd;
		  quat[3] = qx; quat[4] = qy; quat[5] = qz;

		  VAttrList qList   = VCreateAttrList();
		  VBundle   qBundle = VCreateBundle( "bundle", qList, 6*sizeof(float), (VPointer)quat );

		  VSetAttr( vGeoInfo, "qform_code",nullptr,VShortRepn,(VShort) 1 /*NIFTI_XFORM_SCANNER_ANAT*/);
		  VSetAttr( vGeoInfo, "qform", nullptr, VBundleRepn, qBundle );

#ifdef DEBUG
		  cerr << "dbg:" << endl
		       << " quatern_b = " << qb << endl
		       << " quatern_c = " << qc << endl
		       << " quatern_d = " << qd << endl
		       << " qoffset_x = " << qx << endl
		       << " qoffset_y = " << qy << endl
		       << " qoffset_z = " << qz << endl
		       << " dx = " << dx << " dy=" << dy << " dz=" << dz << endl
		       << " qfac=" << qfac << endl;
#endif
		}
	    }
	}
    }

  VAttrList dList   = nullptr;
  VBundle   dBundle = nullptr;
  if( dcm->VoxelWidth() > 0 && dcm->VoxelHeight() > 0 && dcm->SliceThickness() > 0 )
    {
      float* pixDim = (float*)VCalloc(8,sizeof(float));

      pixDim[0] = qfac;
      pixDim[1] = dcm->VoxelWidth();
      pixDim[2] = dcm->VoxelHeight();
      pixDim[3] = dcm->SliceSpacing();
      if( pixDim[3] == -1 )
	{ // If spacing information is unavailable we should better calculate 
	  // the value from coordinates of subsequent slices !
	  // But for the moment this is fine most of the time:
	  pixDim[3] = dcm->SliceThickness();
	}

      pixDim[4] = dcm->TR() / 1000.0;
      pixDim[5] = pixDim[6] = pixDim[7] = 0;

      dList   = VCreateAttrList();
      dBundle = VCreateBundle ("bundle",dList,8*sizeof(float),(VPointer)pixDim);
      VSetAttr(vGeoInfo,"pixdim",nullptr,VBundleRepn,dBundle);
#ifdef DEBUG
      cerr << "dbg: pix dim = " << pixDim[0] << "," << pixDim[1] << "," << pixDim[2] << "," << pixDim[3] << ","
	   << pixDim[4]  << "," << pixDim[5]  << endl;
#endif
    }

  VImage sform = VCreateImage( 1, 4, 4, VFloatRepn );
  VFillImage ( sform, VAllBands, 0);
  VPixel( sform, 0, 0, 0, VFloat) = 1.0;
  VPixel( sform, 0, 1, 1, VFloat) = 1.0;
  VPixel( sform, 0, 2, 2, VFloat) = 1.0;
  VSetAttr( vGeoInfo, "sform_code", nullptr, VShortRepn, (VShort) 0 );
  VSetAttr( vGeoInfo, "sform", nullptr, VImageRepn, sform);

  // and some things not so much related to geometry and coordinates:

  if( dcm->PatientName(str) )
    VSetAttr( vGeoInfo,"patient",nullptr,VStringRepn,(VString)str.c_str());
  // if( dcm->PatientBirthday(str) )
  // 	VSetAttr( vGeoInfo,"birthday",nullptr,VStringRepn,(VString)str.c_str());
  if( dcm->Gender(str) )
    VSetAttr( vGeoInfo,"sex",nullptr,VStringRepn,(VString)str.c_str());
  if( dcm->StudyDescription(str) )
    VSetAttr( vGeoInfo,"study",nullptr,VStringRepn,(VString)str.c_str());
  if( dcm->SeriesDescription(str) )
    VSetAttr( vGeoInfo,"series",nullptr,VStringRepn,(VString)str.c_str());
  if( dcm->ProtocolName(str) )
    VSetAttr( vGeoInfo,"protocol",nullptr,VStringRepn,(VString)str.c_str()); 
  if( dcm->StudyDate(str) )
    VSetAttr( vGeoInfo,"date",nullptr,VStringRepn,(VString)str.c_str());
  if( dcm->ImgOrientation(str) )
    VSetAttr( vGeoInfo,"orientation",nullptr,VStringRepn,(VString)str.c_str());

  double dbl = dcm->TR();
  if( dbl > 0 )
    {
      str = to_string(dbl);
      VSetAttr( vGeoInfo,"repetition_time",nullptr,VStringRepn,(VString)str.c_str());
    }
  dbl = dcm->TE();
  if( dbl > 0 )
    {
      str = to_string(dbl);
      VSetAttr( vGeoInfo,"echoTime",nullptr,VStringRepn,(VString)str.c_str());
    }
  dbl = dcm->FlipAngle();
  if( dbl > 0 )
    {
      str = to_string(dbl);
      VSetAttr( vGeoInfo,"flipAngle",nullptr,VStringRepn,(VString)str.c_str());
    }

  return vGeoInfo;
}


int
main( int argc, char* argv[] )
{
  //
  // check arguments
  //
  if( argc < 2 )
    {
      UsageAndDie( argv[0], "not enough arguments." );
    }

  list<string> inputFiles;	// cmd line specified input file(s)
  string  dicomdirFile;		// we accept one DICOMDIR file
  string  srcDir;		// xor a folder to scan for input files
  bool    verbose   = false;	//
  bool    overwrite = false;	//
  //bool    readOnly  = false;	//
  string  tgtDir = ".";		// folder to write the resulting .v file(s) into
  string  tgtPrefix;		// optional .v file prefix
  i2o_dcm dcm;			// a dicom file object

  struct stat sb;

  int i = 0;

  while ( ++i < argc && argv[i][0] == '-' )
    {
      if( strcmp( argv[ i ], "-v" ) == 0 ) // be verbose
	{
	  verbose = true;
	}
      else if( strcmp( argv[ i ], "-f" ) == 0 ) // be brave and overwrite things if necessary
	{
	  overwrite = true;
	}
      // else if( strcmp( argv[ i ], "-r" ) == 0 ) // dry run, just pretend to do s.th.
      // 	{
      // 	  readOnly = true;
      // 	}
      else if( strcmp( argv[ i ], "-t" ) == 0 ) // the next argument shall be the target directory
	{
	  if( ++i < argc )
	    {
	      tgtDir = argv[i];
	      if( tgtDir.back() == '/' )
		tgtDir.erase( tgtDir.size() - 1 );
	    }
	  else
	    {
	      UsageAndDie( argv[0], "missing argument for \"-t\".");
	    }
	}
      else if( strcmp( argv[ i ], "-out" ) == 0 ) // target prefix
	{
	  if( ++i < argc && argv[i][0] != '-' )
	    {
	      tgtPrefix = argv[i];
	      if( purify( tgtPrefix ) )
		{
		  string issue = "invalid character in outputname \"";
		  issue.append( argv[i] );
		  issue.append( "\"" );
		  UsageAndDie( argv[0], issue );
		}
	    }
	  else
	    {
	      UsageAndDie( argv[0], "insufficient argument for \"-out\".");
	    }
	}
      else if( strcmp( argv[ i ], "-in" ) == 0 )
	{
	  while( ++i < argc && argv[i][0] != '-' )
	    {
	      if( stat( argv[i], &sb ) == 0 )
		{
		  if( S_ISDIR( sb.st_mode ) )
		    {
		      if( srcDir.length() ) // srcDir already set ?
			{
			  UsageAndDie(argv[0], "multiple directories specified.");
			}
		      else
			{
			  if( ! inputFiles.empty() ) // at least one input file already specified ?
			    {
		 	      UsageAndDie( argv[0], "input directory and input files are mutually exclusive." );
			    }
			  else
			    {
			      srcDir = argv[i];
			    }
			}
		    }
		  else
		    {
		      if( srcDir.length() ) //
			{
			  UsageAndDie( argv[0], "input directory and input files are mutually exclusive." );
			}
		      else
			{
			  if( dcm.Read( argv[i] ) ) // is it a dicom ?
			    {
			      if( verbose )
				cerr << "\"" << argv[i] << "\" accepted as a dicom file." << endl;

			      if( dcm.IsDicomDir() ) // is it a DICOMDIR ?
				{
				  if( verbose )
				    cerr << "\"" << argv[i] << "\" really looks like a DICOMDIR file." << endl;

				  if( dicomdirFile.empty() && inputFiles.empty() )
				    {
				      dicomdirFile = argv[i];
				    }
				  else
				    {
				      UsageAndDie( argv[0],
						   "only accepting a single DICOMDIR *or* a list of ordinary dicom files."
						   );
				    }
				}
			      inputFiles.push_back( argv[i] );
			    }
			  else
			    {
			      cerr << "-- \"" << argv[i] << "\" not accepted as a dicom file." << endl;
			      exit(1);
			    }
			}
		    }
		}
	      else
		{
		  cerr << "stat(\"" << argv[i] << "\") failed: " << strerror(errno) << endl;
		}
	    }
	  // if input list is in the middle of the road '-' breaks the while loop
	  // and the argument counter must be adjusted
	  if( i < argc && argv[i][0] == '-' )
	    --i;
	}
      else
	{
	  string issue =  "unknown parameter \"";
	  issue.append( argv[i] );
	  issue.append( "\"" );
	  UsageAndDie( argv[0], issue );
	}
    }

  //
  // so much to the command line, see what we got
  //

  if( srcDir.empty() && inputFiles.empty() )
    {
      UsageAndDie( argv[0], "don't know what to do." );
    }

  //
  // arguments ok, short roundup if requested
  if( verbose )
    {
      if( srcDir.length() )
	cerr << "Source directory is \"" << srcDir << "\"." << endl;
      else
	cerr << inputFiles.size() << " source file(s) specified." << endl;
      cerr << "Target directory is \"" << tgtDir << "\"." << endl;
      cerr << "File overwriting \"" << ( overwrite ? "on" : "off" ) << "\"." << endl;
      cerr << "Outputfile prefix is \"" << (tgtPrefix.length() ? tgtPrefix : "" ) << "\"." << endl;
    }

  //
  // collect input files
  //

  DicomDirHead dcmDir;
  unsigned     totalFileCnt;

  // three cases:

  // we have a directory name or a DICOMDIR file
  if( srcDir.length() || dicomdirFile.length() )
    {
      string src = srcDir.length() ? srcDir : dicomdirFile;

      totalFileCnt = ScanDicomDir( src, dcmDir );
      if( verbose )
	cerr << "... found " << totalFileCnt << " file(s) in \"" << src << "\"" << endl;
    }
  else
    // or a list of files
    {
      totalFileCnt = ScanDicomFiles( inputFiles, dcmDir );

      if( totalFileCnt != inputFiles.size() )
	{
	  cerr << "-- " << totalFileCnt << " file(s) found, you specified " << inputFiles.size() << endl;
	  exit(1);
	}
      if( verbose )
	cerr << "... " << totalFileCnt << " file(s) found." << endl;
    }

  // either way, we got it sorted
  if( totalFileCnt > 0 )
    {
      list<DicomPatient>& patList = dcmDir.patientList;
      
      //
      // PATIENT
      //
      for( auto &patient : patList )
	{
	  if( verbose )
	    cerr << "Patient \"" << patient.patientName << "\": nr of studies = "
		 << patient.studyList.size() << endl;

	  //
	  // STUDY
	  //
	  for( auto &stdy : patient.studyList )
	    {
	      if( verbose )
		cerr << " study ID = " << stdy.studyID << " \"" << stdy.description << "\": nr of series = "
		     << stdy.seriesList.size() << endl;

	      //
	      // SERIES
	      //
	      for( auto series : stdy.seriesList )
		{
		  if( verbose )
		    cerr << "  series #" << series.seriesNr << " \"" << series.description
			 << "\": nr of files/images = " << series.fileList.size() << endl;

		  //************************************************************************************
		  // loop over all files in this series to collect the data and the header informations,
		  // name and create one .v, write header and data.
		  //************************************************************************************

		  int  cnt    = 0;
		  bool mosaic = false;
		  bool orientationMismatch = false;

		  i2o_dcm* dcmFiles[ series.fileList.size() ];

		  //
		  // FILES
		  //
		  for( auto &dcmFile : series.fileList )
		    {
		      string   fname = dcmFile.path;

		      i2o_dcm* dcm  = new i2o_dcm;

		      dcmFiles[cnt] = dcm;

		      if( verbose )
			cerr << "   #" << dcmFile.imgNr << ": \"" << fname << "\"" << endl;

		      if( dcm->Read( fname ) )
			{
			  
			  // if( dcm->IsDicomDir() ) // ??? is this still necessary ??? Check i2o_dcm.cc.
			  //   {
			  //     cerr << "skipped DICOMDIR file \"" << fname << "\" ..." << endl;
			  //     continue;
			  //   }

			  if( dcm->DataSize() == 0 || dcm->Cols() <=0 || dcm->Rows() <= 0 )
			    {
			      cerr << "Series " << dcm->SeriesNr()
				   << ": \"" << fname
				   << "\" apparently has no data or invalid size value(s), skipped ..." << endl;
			      continue;
			    }

			  if( dcm->SamplesPerPixel() != 1 )
			    {
			      cerr << "Series " << dcm->SeriesNr()
				   << ": \"" << fname
				   << "\" apparently contains an unhandled data type, skipped ..." << endl;
			      continue;
			    }

			  cnt++;

			  if( dcm->Mosaic())
			    {
			      mosaic = true;
			      if( verbose )
				{
				  cerr << "\t MOSAIC: total "
				       << dcm->TotalCols() << "x" << dcm->TotalRows() << " : "
				       << dcm->Slices() << " slices " << dcm->Cols() << "x" << dcm->Rows();
				}
			    }
			  else
			    {
			      if( mosaic ) // paranoid ???
				{
				  cerr << "Mosaic status changed inside series. Ignoring..." << endl;
				  mosaic = false; // BUT HOW CAN WE PROCEED ???
				}
			      if( verbose )
				cerr << "\t "
				     << dcm->TotalCols() << "x" << dcm->TotalRows() << "x" << dcm->Slices();
			      
			    }
			  if( verbose )
			    {
			      cerr << ", location=" << dcm->SliceLocation()
				   << ", orientation=" << dcm->ImgOrientation()
				   << ", " << dcm->OrderOfSlices()
				   << ", distance=" << dcm->DistanceFromOrigin()
				   << ", spacing=" << dcm->SliceSpacing()
				   << endl;
			      
#ifdef DEBUG
			      string str;
			      vector<double> v3d = dcm->ImgOrientRow();

			      str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			      cerr << "row: " << str;

			      v3d = dcm->ImgOrientCol();
			      str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			      cerr << ", col: " << str;

			      v3d = dcm->ImgOrientNormal();
			      str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			      cerr << ", normal: " << str;

			      v3d = dcm->ImgPosition();
			      str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			      cerr << ", pos: " << str
				   << endl;
#endif
			    }

			  if( cnt > 1 )
			    {
			      i2o_dcm* previousDcm = dcmFiles[ cnt - 2 ];
			      
			      // now we take the chance to check attribute consistancy ...
			      // But there should certainly be a little bit more about it here.
			      
			      // and this should better be done in the DicomDirHead code anyway ...
			      
			      // if( dcm->ImgPosition() == previousDcm->ImgPosition() )
			      //   { // has to happen for functional datasets ...
			      //     cerr << "!!! Identical coordinates for subsequent images. "
			      // 	// << dcm->SliceLocation() << " == "
			      // 	// << previousDcm->SliceLocation()
			      // 	   << endl;
			      //     splitSeries = true;
			      //   }
			      
			      if( dcm->PhaseEncodingDirPositive() != previousDcm->PhaseEncodingDirPositive() )
				{
				  cerr << "Series " << dcm->SeriesNr()
				       << ", \"" << fname <<"\": PhaseEncodingDirectionPositive mismatch. "
				       << dcm->PhaseEncodingDirPositive() << " <> "
				       << previousDcm->PhaseEncodingDirPositive() << endl;
				}


			      if( dcm->ImgOrientation() != previousDcm->ImgOrientation() )
				{ // might be a localizer series ...
				  cerr << "Series " << dcm->SeriesNr()
				       << ", \"" << fname << "\": Orientation mismatch (localizer?). "
				       << dcm->ImgOrientation() << " <> "
				       << previousDcm->ImgOrientation()
				       << ". Ignored."<< endl;
				  orientationMismatch = true;
				}
			      
			      const vector<double> t0 = previousDcm->AcquisitionRefTimes();
			      
			      if( t0.size() )
				{
				  const vector<double> t1 = dcm->AcquisitionRefTimes();
				  
				  if ( t0.size() == t1.size() )
				    {
				      if( t0 != t1 )
					{
					  bool significantly = false;
					  const float epsilon = 2.5 + 0.0001; // ms, the MR's time raster
					  
					  for( unsigned i = 0; i < t0.size(); i++ )
					    if( fabs( t0[i] - t1[i] ) > epsilon )
					      significantly = true;
					  
					  if( significantly )
					    {
					      // this may happen erroneously, go ahead, no exit
					      // but report the details
					      cerr << cnt
						   << ":\"AcquisitionRefTimes\" differ for subsequent files:"
						   << endl;
					      
					      for( unsigned i = 0; i < t0.size(); i++ )
						if( fabs( t0[i] - t1[i] ) > epsilon )
						  cerr << i << ": " << t0[i] << " <> " << t1[i] << endl;
					    }
					}
				    }
				  else
				    {
				      cerr << cnt
					   << ": number of \"AcquisitionRefTimes\" differ for subsequent files: "
					   << t0.size() << " <> " << t1.size() << endl;
				      exit(1);
				    }
				}

			      int n0 = dcm->DataSize();
			      int n1 = previousDcm->DataSize();
			      if( n0 != n1 )
				{
				  cerr << cnt << ": \"data size\" differ for subsequent files,"
				       << n0 << "!=" << n1 << endl;
				  exit(1);
				}
			      
			      n0 = dcm->Slices();
			      n1 = previousDcm->Slices();
			      
			      if( n0 != n1 )
				{
				  cerr << cnt << ": number of slices differ for subsequent files,"
				       << n0 << "!=" << n1 << endl;
				  exit(1);
				}
			    }
			}
		      else
			{
			  cerr << "S.th. is wrong with \"" << fname << "\".";
			  if( dcm->Msg().size())
			    cerr << dcm->Msg() << endl;
			  else
			    cerr << endl;
			  exit(1);
			}
		    }
		  
		  if( cnt == 0 )
		    continue;
		  
		  // everything read in and "checked",
		  // now try to do s.th. useful:
	      
		  //
		  if( dcmFiles[0]->ImgOrientation() == "sagittal" )
		    {
		      if( verbose && cnt > 1 && !orientationMismatch )
			cerr << "reversing order of sagittal slices" << endl;
		      for( int i = 0; i < cnt/2; ++i )
			{
			  i2o_dcm* t = dcmFiles[ i ];
			  dcmFiles[ i ] = dcmFiles[ cnt - 1 - i ];
			  dcmFiles[ cnt - 1 - i ] = t;
			}
		    }
	      
		  i2o_dcm* dcm = dcmFiles[0];
	      
		  int    rows  = dcm->Rows();
		  int    cols  = dcm->Cols();
		  double tr    = dcm->TR();
		  string voxel_str;
	      
		  if( dcm->VoxelWidth() > 0 && dcm->VoxelHeight() > 0 && dcm->SliceThickness() > 0 )
		    {
		      voxel_str = to_string( dcm->VoxelWidth() ) + " "
			+ to_string( dcm->VoxelHeight() ) + " ";
		      if( dcm->SliceSpacing() == -1 )
			voxel_str += to_string( dcm->SliceThickness() );
		      else
			voxel_str += to_string( dcm->SliceSpacing() );
		    }
	      
		  string vistaFilename;
	      
		  FILE*  vOut = CreateOutputFileOrDie( tgtDir, tgtPrefix, dcm, overwrite, vistaFilename );
	      
		  VAttrList vAttrList = VCreateAttrList();
	      
		  //
		  // right now, we handle two cases in principle
		  //
	      
		  if( mosaic && cnt > 1 ) // cnt == 1 is a special case handled in the else part
		    {
		      const unsigned timePoints = cnt;
		      const unsigned slices     = dcm->Slices();
		  
		      VAttrList vGeoInfo = CreateGeoInfo( dcm, 4, slices, timePoints );
		      VAppendAttr( vAttrList, "geoinfo", nullptr, VAttrListRepn, vGeoInfo);
		  
		      VImage* vSliceImages = (VImage *) VCalloc( slices, sizeof( VImage ) );
		  
		      const vector<double> acquisitionRefTime = dcmFiles[cnt-1]->AcquisitionRefTimes();
		  
		      for( unsigned slice = 0; slice < slices; slice++ )
			{
			  vSliceImages[slice] = VCreateImage( timePoints, rows, cols, VShortRepn );

			  VFillImage( vSliceImages[slice], VAllBands, 0 );
			  VSetAttr( vSliceImages[slice]->attributes, "repetition_time", nullptr, VDoubleRepn, tr );
		      
			  if( acquisitionRefTime.size() > slice )
			    VSetAttr( vSliceImages[slice]->attributes, "slicetime", nullptr, VDoubleRepn,
				      acquisitionRefTime[slice] );
		      
			  VSetAttr( vSliceImages[slice]->attributes, "voxel", nullptr, VStringRepn,
				    (VString)voxel_str.c_str());
		      
			  char* tgtSliceData = (char*)vSliceImages[slice]->data;
		      
			  for( unsigned t = 0; t < timePoints; ++t )
			    {
			      int  xdim,ydim;
			      int  totalSliceSize;
			  
			      dcmFiles[t]->GetSliceData( slice, xdim, ydim, totalSliceSize, tgtSliceData );
			      tgtSliceData += totalSliceSize;
			    }
		      
			  VAppendAttr(vAttrList,"image",nullptr,VImageRepn,vSliceImages[slice]);
		      
			}
		    }
		  else
		    {
		      // might be either one multislice multifile dataset ( e.g. anatomical volume )
		      // or a number of single slice files in one series ( e.g. localizer )
		      // ...

		      // ... or the special case: a single mosaic file
		      const int slices = ( cnt == 1 ? dcm->Slices() : cnt );
		      
		      string    str;

		      VAttrList vGeoInfo = CreateGeoInfo( dcm, 3, slices, 1 );
		  
		      VAppendAttr( vAttrList, "geoinfo", nullptr, VAttrListRepn, vGeoInfo);
		  
		      //
		      VImage vImg = VCreateImage( slices, rows, cols, VShortRepn );
		      if (!vImg)
			{
			  cerr << "VCreateImage( "
			       << slices << "," << rows << "," << cols << ", VShortRepn ) failed." << endl;
			  exit(1);
			}
		  
		      if( dcm->PatientName(str) )
			VSetAttr(vImg->attributes,"patient",nullptr,VStringRepn,(VString)str.c_str());
		      // if( dcm->PatientBirthday(str) )
		      // 	VSetAttr(vImg->attributes,"birthday",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->Gender(str) )
			VSetAttr(vImg->attributes,"sex",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->StudyDescription(str) )
			VSetAttr(vImg->attributes,"study",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->SeriesDescription(str) )
			VSetAttr(vImg->attributes,"series",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->ProtocolName(str) )
			VSetAttr(vImg->attributes,"protocol",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->StudyDate(str) )
			VSetAttr(vImg->attributes,"date",nullptr,VStringRepn,(VString)str.c_str());
		      if( dcm->ImgOrientation(str) )
			VSetAttr(vImg->attributes,"orientation",nullptr,VStringRepn,(VString)str.c_str());
		  
		      if( dcm->VoxelWidth() > 0 && dcm->VoxelHeight() > 0 && dcm->SliceThickness() > 0 )
			{
			  str = to_string( dcm->VoxelWidth() ) + " "
			    + to_string( dcm->VoxelHeight() ) + " ";
			  if( dcm->SliceSpacing() == -1 )
			    str += to_string( dcm->SliceThickness() );
			  else
			    str += to_string( dcm->SliceSpacing() );
			  VSetAttr(vImg->attributes,"voxel",nullptr,VStringRepn,(VString)str.c_str());
			}

		      vector<double> v3d = dcm->ImgOrientRow();
		      if( v3d.size() == 3 )
			{
			  str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			  //VSetAttr(vImg->attributes,"rowVec",nullptr,VStringRepn,(VString)str.c_str());
			  VSetAttr(vImg->attributes,"ImgOrientationRow",nullptr,VStringRepn,(VString)str.c_str());
			}
		      v3d = dcm->ImgOrientCol();
		      if( v3d.size() == 3 )
			{
			  str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			  //VSetAttr(vImg->attributes,"columnVec",nullptr,VStringRepn,(VString)str.c_str());
			  VSetAttr(vImg->attributes,"ImgOrientationCol",nullptr,VStringRepn,(VString)str.c_str());
			}
		      v3d = dcm->ImgOrientNormal();
		      if( v3d.size() == 3 )
			{
			  str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			  //VSetAttr(vImg->attributes,"sliceVec",nullptr,VStringRepn,(VString)str.c_str());
			  VSetAttr(vImg->attributes,"ImgOrientationNormal",nullptr,VStringRepn,(VString)str.c_str());
			}
		  
		      v3d = dcm->ImgPosition();
		      if( v3d.size() == 3 )
			{
			  str = to_string(v3d[0]) + " " + to_string(v3d[1]) + " " + to_string(v3d[2]);
			  //VSetAttr(vImg->attributes,"indexOrig",nullptr,VStringRepn,(VString)str.c_str());
			  VSetAttr( vImg->attributes, "ImagePositionPatient", nullptr, VStringRepn, (VString)str.c_str() );
			}
		      double dbl = dcm->TR();
		      if( dbl > 0 )
			{
			  str = to_string(dbl);
			  VSetAttr(vImg->attributes,"repetition_time",nullptr,VStringRepn,(VString)str.c_str());
			}
		      dbl = dcm->TE();
		      if( dbl > 0 )
			{
			  str = to_string(dbl);
			  VSetAttr(vImg->attributes,"echoTime",nullptr,VStringRepn,(VString)str.c_str());
			}
		      dbl = dcm->FlipAngle();
		      if( dbl > 0 )
			{
			  str = to_string(dbl);
			  VSetAttr(vImg->attributes,"flipAngle",nullptr,VStringRepn,(VString)str.c_str());
			}
		  
		      char* tgtSliceData = (char*)vImg->data ;
		  
		      if( mosaic ) // the special case, one timepoint, a single volume
			{
			  for( unsigned slice = 0; slice < slices; slice++ )
			    {
			      int xdim,ydim;
			      int totalSliceSize;
			      
			      dcm->GetSliceData( slice, xdim, ydim, totalSliceSize, tgtSliceData );
			      tgtSliceData += totalSliceSize;
			    }			  
			}
		      else
			for( int slice = 0; slice < slices; slice++ )
			  {
			    dcmFiles[slice]->GetRawData( tgtSliceData );
			    tgtSliceData += dcmFiles[slice]->DataSize();
			  }
		      VAppendAttr( vAttrList, "image", nullptr, VImageRepn, vImg );
		  
		    }
	      
		  VWriteFile( vOut, vAttrList );
		  fclose( vOut );

		  // OPEN ISSUE
		  // vista "memory management" is just a mess and there is definitely much more to do
		  // than this:
		  VDestroyAttrList( vAttrList );
	      
		  //if( verbose )
		  cerr << "\"" << vistaFilename << "\" created." << endl;
		}
	    }
	}
    }
  else
    cerr << "Conversion list empty." << endl;
  
  cerr << basename( argv[0] ) << ": done." << endl;
  
}
