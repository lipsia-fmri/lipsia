//
// i2o_dcm.hh
// created at Mon Jan  7 14:06:10 2013 (as a renewal from older material)
// in /home/karda/src/i2o/
// by \bk: bernd.kardatzki@med.uni-tuebingen.de
//
#ifndef I2O_DCM_HH
#define I2O_DCM_HH

#include <iostream>
#include <string>
#include <list>
#include <vector>

using namespace std;

//
// simple infrastructure for either scanning a single "DICOMDIR" file
// or a directory of dicom files or a list of dicom files:
struct slicePreview { 
  slicePreview(): rows(0),columns(0),data(NULL){};
  ~slicePreview(){ if( data ) delete [] data; };
  
  int rows; int columns; unsigned char* data; 
};

struct DicomFile { 
  DicomFile(int nr, string p):imgNr(nr),path(p){};
  int imgNr; string path; slicePreview preview;
};

struct DicomSeries {
  DicomSeries(int nr, string descr):seriesNr(nr), description(descr){};
  int seriesNr; string description;
  list<DicomFile> fileList;  // this list is sorted
};

struct DicomStudy {
  DicomStudy(int nr, string descr):studyID(nr),description(descr){};
  int studyID; string description;
  list<DicomSeries> seriesList; // this list is sorted
};

struct DicomPatient {
  DicomPatient(string pn):patientName(pn){}; 
  string patientName;
  list<DicomStudy> studyList; 
};

struct DicomDirHead {
  string baseDir;
  list<DicomPatient> patientList; 
};

//
int				// return total number of dicom images found
ScanDicomDir( string path,	// DICOMDIR or directory name to scan
	      DicomDirHead& dcmdir  // structure to return the list of patients (see above)
	      );

int				// return total number of dicom images found
ScanDicomFiles( list<string>& inputFiles, // list of files to put together into a details struct:
		DicomDirHead& dcmdir  // structure to return the list of patients (see above)
		);


// ... and now to s.th. completely different:
// this is about reading a dicom file fast 
// and collecting the relevant dicom tags into class data members
//

enum OrderOfSlices_T
  {
    UndefinedOrder = 0,
    AscendingOrder = 1,
    DescendingOrder = 2 ,
    InterleavedOrder = 4,
    NoOrder
  };

class i2o_dcm
{
public:
  i2o_dcm();
  i2o_dcm( string fn );

  ~i2o_dcm();

  string Msg(); // if sth goes wrong, Msg() hopefully tells us why
  
  string Info();		// some (very rudimentary) information about the dicom file
  
  bool   Read( string fn, bool deferred=true );	// open the file and read it in
  bool   Write( string fn );	// write the dicom file
  
  bool   IsDicomDir();		// either directly read from a dicom tag or concluded heuristically
  // in case of IsDicomDir() returns true, this creates a a list of lists:
  int    DicomDir( DicomDirHead& ); // returns total # of images referenced/listed

  // space and time savers
  //
  // member string variable and access functions
#define I2O_STR(Name)							\
  private:		/* PRIVATE DATA MEMBER */			\
  string m_##Name;							\
public:			/* PUBLIC ACCESS METHODS */			\
 const string& Name( void ){						\
   if( m_##Name.size() == 0 )						\
     if( !m_fullScanDone)						\
       FullScan();							\
   return m_##Name;	/* might be empty ! */				\
 }									\
 bool Name( string& str ){						\
   Name();								\
   if( m_##Name.size() ){						\
     str = m_##Name;							\
     return true;							\
   }									\
   else									\
     return false;							\
 }
  
  I2O_STR(ImageType);
  // == "private string m_ImageType";
  //    "public  const string& ImageType(void);"
  //    "public  bool ImageType( string& );"
  I2O_STR(PatientName);
  I2O_STR(Gender);
  I2O_STR(PatientID);
  I2O_STR(PatientAge);
  I2O_STR(PatientWeight);
  I2O_STR(PatientBirthday);	// yyyymmdd
  I2O_STR(PatientComment);	// this is the place to put the "@" or the "_" 
  I2O_STR(PatientPosition);
  I2O_STR(StudyDate);		// yyyymmdd
  I2O_STR(StudyTime);		// hhmmss.ffffff
  I2O_STR(SeriesDate);		// yyyymmdd
  I2O_STR(SeriesTime);		// hhmmss.ffffff
  I2O_STR(ContentDate);		// yyyymmdd
  I2O_STR(ContentTime);		// hhmmss.ffffff
  I2O_STR(AcquisitionDate);	// yyyymmdd
  I2O_STR(AcquisitionTime);	// hhmmss.ffffff
  I2O_STR(ManufacturerModelName);
  I2O_STR(StudyDescription);
  I2O_STR(SeriesDescription);
  I2O_STR(ProtocolName);
  I2O_STR(ImgOrientation);
  I2O_STR(Modality);
  I2O_STR(InPlanePhaseEncodingDir);
  I2O_STR(StudyInstanceUID);
  I2O_STR(SeriesInstanceUID);
  I2O_STR(SWVersion);
#undef I2O_STR
  
#define I2O_DBL_INT(TypeDblOrInt,Name)					\
  private:		/* PRIVATE DATA MEMBER */			\
  TypeDblOrInt m_##Name/*=-1*/;						\
public:			/* PUBLIC ACCESS METHODS */			\
 TypeDblOrInt Name( void ){						\
   if( m_##Name < 0 )							\
     if( !m_fullScanDone)						\
       FullScan();							\
   return m_##Name;     /* might return -1 */				\
 }									\
 bool Name( string& s ){						\
   Name();								\
   if( Name() > 0 ){							\
     s = to_string(m_##Name);						\
     return true;							\
   }									\
   else									\
     return false;							\
 }

  I2O_DBL_INT(int,Cols);
  I2O_DBL_INT(int,TotalCols);
  I2O_DBL_INT(int,Rows);
  I2O_DBL_INT(int,TotalRows);
  I2O_DBL_INT(int,Slices);
  I2O_DBL_INT(int,FoVx);
  I2O_DBL_INT(int,FoVy);
  I2O_DBL_INT(int,SamplesPerPixel);
  
  I2O_DBL_INT(int,StudyID);
  I2O_DBL_INT(int,SeriesNr);
  I2O_DBL_INT(int,AcquisitionNr);
  I2O_DBL_INT(int,InstanceNr); // formerly known as image number
  I2O_DBL_INT(int,PhaseEncodingDirPositive);
  I2O_DBL_INT(double,TE);
  I2O_DBL_INT(double,TR);
  I2O_DBL_INT(double,TI);
  I2O_DBL_INT(double,FlipAngle);
  I2O_DBL_INT(double,VoxelWidth); // in mm
  I2O_DBL_INT(double,VoxelHeight);
  I2O_DBL_INT(double,SliceThickness);
  I2O_DBL_INT(double,SliceSpacing); // center to center
  I2O_DBL_INT(double,SliceLocation); // usefull for non mosaics only
#undef I2O_DBL_INT

#define I2O_VECTOR(Name)						\
  private:		/* PRIVATE DATA MEMBER */			\
  vector<double> m_##Name;						\
public:			/* PUBLIC ACCESS METHODS */			\
 vector<double>& Name( void )						\
  {									\
    if( m_##Name.size() == 0 )						\
      if( !m_fullScanDone )						\
	FullScan();							\
    return m_##Name;							\
  }							

  struct I2O_POS { double x,y,z; };
#define I2O_VECTOR3(Name)						\
  private:		/* PRIVATE DATA MEMBER */			\
  vector<I2O_POS> m_##Name;						\
public:			/* PUBLIC ACCESS METHODS */			\
 vector<I2O_POS>& Name( void )						\
  {									\
    if( m_##Name.size() == 0 )						\
      if( !m_fullScanDone )						\
	FullScan();							\
    return m_##Name;							\
  }							

  I2O_VECTOR(SliceNormal);
  
  I2O_VECTOR(ImgPosition); // center of the first voxel transmitted, (0x20,0x32) ImagePositionPatient
  I2O_VECTOR(ImgOrientRow);    // direction cosine
  I2O_VECTOR(ImgOrientCol);    // direction cosine
  I2O_VECTOR(ImgOrientNormal); // crossproduct ot the two vectors above (calculated)

  I2O_VECTOR(AcquisitionRefTimes); // AcquisitionRefTime for each slice
  I2O_VECTOR3(ImgPositionSlice); // ImagePositionPatient for each slice

  I2O_VECTOR(DiffusionGradDirection);
  
#undef I2O_VECTOR
#undef I2O_VECTOR3
  
  bool  Mosaic() { return ( Slices() > 1 && TotalCols() > Cols() ); }

  bool  Ascending()  { return m_orderOfSlices == AscendingOrder; }
  bool  Interleaved(){ return m_orderOfSlices == InterleavedOrder; }
  bool  Descending() { return m_orderOfSlices == DescendingOrder; }

  string OrderOfSlices();

  double DistanceFromOrigin();
  
  //
  // payload access
  
  int DataSize();		// in bytes

  //
  // unique_ptr<> seems to be appropriate for the following once -std=c++11 becomes legal for lipsia
  // so things will definitely change here:
  
  // it's the users responsibility to free[] the data later because GetData() does the allocation 
  // ... and always check "xdim", "ydim", "zdim", "size" !
  void*  GetRawData( int& xdim, int& ydim, int& zdim /*(#slices)*/, int& size /* in bytes */);
  // and in case you prefer to do the allocation all by yourself:
  bool   GetRawData( void* dest ); // just a memcpy basically, use Cols(), Rows() and Slices()
  // direct pointer to the payload, no memcpy, CAUTION: valid only during object lifetime
  void*  RawDataDirect( int& xdim, int& ydim, int& zdim /*(#slices)*/, int& totalSize /* in bytes */);
  // return pointer to 3d data slice, 'dest' will be allocated if necessary 
  void*  GetSliceData( int sliceNr, int& xdim, int& ydim, int& totalSliceSize, void* dest = NULL );
  // return pointer to 3d data volume, 'dest' will be allocated if necessary
  void*  Mosaic2Vol( int& xdim, int& ydim, int& zdim, int& destSize, void* dest = NULL );  
  
  // Direct injection, handle with care ...
  bool    SetData( void*, int totalByteSize );

#if __cplusplus < 201103L
#warning Compile with -std=c++11 and consider updating your C++ Compiler
#endif
  size_t Datatype_hash_code();
  string Datatype_C_Name();
  
  // ... again, use it with care ( and always respect "size" )
  char*   Grep( int grp, int ele, int& size );
  
 private:
  
  // methods
  void    Init(); // trashed as soon as compiling with -std=c++11 becomes legal under the lipsia laws
    
  bool    Open( string fn );
  void    Close( void );
  
  void    Dump();
  void    FullScan(); // the complete file read is usually deferred, until FullScan is called
  
  void    FigureOutDataType();
  void    FigureOutImgOrientation();

  // separate the dicom group parsing
  void    Grp08(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp10(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp18(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp19(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp20(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp28(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp29(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp51(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp7FE0(char* b, u_int16_t element, int32_t len, char VR[3]);
  void    Grp7FE1(char* b, u_int16_t element, int32_t len, char VR[3]);
  
// private data members:

  string  m_filename;
  off_t   m_filesize/*=-1*/;

  bool    m_fullScanDone/*=false*/;  // remembers if we're already gone through the whole thing or not
  char*   m_readStoppedHere/*=nullptr*/;  // and if not, where did we stop

  char*   m_firstDcmTag/*=nullptr*/;  // this is where the fun starts, usually @128+4

  void*   m_mm/*=nullptr*/;		// base adress, the file is mmaped() instead of read()
  // quick and simple mechanism to assure that only a small part of the file must been touched
  // to accept or reject a file as a valid dicom
  // DOBEDO:
  // void* zz;	        // in case of a compressed dicom file, mmap is not applicable
  // void* rr;          // mmap() might fail, read() is our friend
  // void* firstByte;   // regardless, this is a pointer to the very beginning

  int     m_isDicomdir/*=-1*/;	// s.th.special, the triple truth

  string  m_msgTxt; // short term event memory, plaintext status description
  
  void*   m_payload /*=nullptr*/;   // i.e. the DATA ( will be a unique_ptr in the very near future ... )
  int32_t m_payloadSize /*=0*/;     // in bytes, int32_t as the maximal data size a single dicom file can carry
  bool    m_littleEndian /*=true*/; // hopefully true, otherwise we're in trouble and have to care about endianess

  OrderOfSlices_T m_orderOfSlices/*=UndefinedOrder*/;

  int     m_echoNr/*=-1*/;  
  int     m_averages/*=-1*/;
  int     m_bValue/*=-1*/;
  int     m_acquisitionMatrix/*=-1*/;

  size_t  m_dataTypeHashCode/*=0*/;
  bool    m_signedDataValues/*=false*/;
  short   m_bitsAllocated/*=-1*/;
  short   m_bitsStored/*=-1*/;
  short   m_highBit/*=-1*/;
  short   m_planarConfiguration/*=-1*/;
  short   m_numberOfFrames/*=-1*/;

};



// little helper
inline ostream& operator << (ostream& os, i2o_dcm& dcm)
{
  return os << dcm.Info();
}

#endif
