//
// i2o_dcm.cc
// created at Mon Jan  7 13:52:26 2013
// in /home/karda/src/i2o/
// by \bk: bernd.kardatzki@med.uni-tuebingen.de
//

#define nDEBUG
#define nDEBUGCSA

#include <iostream>
#include <climits>
#include <cerrno>
#include <cstring>
#include <vector>
#include <typeinfo>
#include <complex>

#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <libgen.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <stdint.h>

#include "i2o_dcm.h"


using namespace std;


//
// some little helpers:

// remove trailing spaces
char* trim(char* c, char t = ' ' )
{
  char* p;

  for( p=c ; *p; p++)
    ;
  for(--p; *p == t; p--)
    *p = '\0';
  return c;
}

// same with a string object
void trim( string& s, char t = ' ' )
{
  if( unsigned i = s.length() )
    {
      while( s[i-1] == t )
	--i;
      if( s.length() > i )
	s.resize( i );
    }
}

//
// Interpret byte stream as a dicom element.
// (we dont care about the byte order anymore, the universe is x86 everywhere ...)
//
char*		  // return address of data part or nullptr on failure
DcmElement(
	   char*       v,	// start here
	   u_int16_t&  grp,	// return grp nr
	   u_int16_t&  el,	// return ele nr
	   int32_t&    len,	// return data size (value might be negative)
	   char        VR[3]	// return VR
	   )
{
  char* pData = nullptr; // pointer to the data which follows the dicom meta stuff
  
  VR[0] = '\0';
  
  // special dicom perversion witnessed in a philips file
  // if( *(u_int32_t*) v == 0xf0f0f0f0 )
  //   {			     // but this should be ancient history ...
  //     cerr << " !!! " << __FILE__ << __LINE__
  // 	   << ": dicom format oddity, expect inconsistancies and heavy weather." << endl;
  //     v += 4;
  //   }

  // assertions
  if( ( grp = *(u_int16_t*)( v ) ) < 0 ) // zero is a valid group number
    return nullptr;
  if( ( el = *(u_int16_t*)( v + 2 ) ) < 0 ) // element zero is the group length (rarely found)
    return nullptr;

  //
  // explicit value representation ?
  //
#define VR_IS( a, b ) ( v[4] == a && v[5] == b )

  if( VR_IS(    'O', 'B' )
      || VR_IS( 'O', 'W' )
      || VR_IS( 'O', 'F' )
      || VR_IS( 'S', 'Q' )
      || VR_IS( 'U', 'N' )
      || VR_IS( 'U', 'T' )
      )
    {
      len = *(int32_t*)( v + 8 );
      pData = v + 12;
      VR[0] = v[4];
      VR[1] = v[5];
    }
  else if( VR_IS(    'A', 'E')
	   || VR_IS( 'A', 'S' ) || VR_IS( 'A','T' ) || VR_IS( 'C', 'S' )
	   || VR_IS( 'D', 'A' ) || VR_IS( 'D','S' ) || VR_IS( 'D', 'T' )
	   || VR_IS( 'F', 'L' ) || VR_IS( 'F','D' ) || VR_IS( 'I', 'S' )
	   || VR_IS( 'L', 'O' ) || VR_IS( 'L','T' ) || VR_IS( 'P', 'N' )
	   || VR_IS( 'S', 'H' ) || VR_IS( 'S','L' ) || VR_IS( 'S', 'S' )
	   || VR_IS( 'S', 'T' ) || VR_IS( 'T','M' ) || VR_IS( 'U', 'I' )
	   || VR_IS( 'U', 'L' ) || VR_IS( 'U','S' )
	   )
    {
      len = *(int16_t*)( v + 6 );
      pData = v + 8;
      VR[0] = v[4];
      VR[1] = v[5];
    }
  else
    {
      if( ( len = *(int32_t*)( v + 4 ) ) >= 0 ){	// length 0 is valid
	pData = v + 8;
	if( len == 2){
#undef VR_IS
#define VR_IS(a,b) ( pData[0]==a && pData[1]==b )
	  if( VR_IS(    'O', 'B' )
	      || VR_IS( 'O', 'F' )
	      || VR_IS( 'O', 'W' )
	      || VR_IS( 'S', 'Q' )
	      || VR_IS( 'U', 'N' )
	      || VR_IS( 'U', 'T' )
	      )
	    {
	      len = *(int32_t*)( pData + 2 );
	      pData += 6;
	      VR[0] = pData[0];
	      VR[1] = pData[1];
	    }
	  else if( VR_IS(    'A', 'E' )
		   || VR_IS( 'A', 'S' ) || VR_IS( 'A', 'T') || VR_IS( 'C', 'S')
		   || VR_IS( 'D', 'A' ) || VR_IS( 'D', 'S') || VR_IS( 'D', 'T')
		   || VR_IS( 'F', 'L' ) || VR_IS( 'F', 'D') || VR_IS( 'I', 'S')
		   || VR_IS( 'L', 'O' ) || VR_IS( 'L', 'T') || VR_IS( 'P', 'N')
		   || VR_IS( 'S', 'H' ) || VR_IS( 'S', 'L') || VR_IS( 'S', 'S')
		   || VR_IS( 'S', 'T' ) || VR_IS( 'T', 'M') || VR_IS( 'U', 'I')
		   || VR_IS( 'U', 'L' ) || VR_IS( 'U', 'S')
		   )
#undef VR_IS
	    {
	      len = *(int32_t*)( pData + 2 );
	      pData += 6;
	      VR[0] = pData[0];
	      VR[1] = pData[1];
	    }
	}
      }
      else {
	//
	if ( len == -1 ){ // might be an "unknown sequence length" or a delimiter or whatever
	  pData = v + 8;
	}
      }
    }

  if( pData
      && ( pData[0] == (char)0xf0 )
      && ( pData[1] == (char)0xf0 )
      && ( pData[2] == (char)0xf0 )
      && ( pData[3] == (char)0xf0 )
      )
    pData += 4;

  //
  return pData;
}


// Siemens "private": _C_ommon _S_iemens _A_rchitecture
// A substructured binary dicom data element, reverse engineered so take with care (see below)
struct CSAItem
{
  string       name;
  string       vr;      // value representation
  list<string> value;   // array of values, apparently everything is stored in ASCII

  CSAItem() {}
  CSAItem( const CSAItem& other )
  {
    if( &other != this ){
      if( other.name.length() ){
	name    = other.name;
	vr      = other.vr;
	if( value.size() )
	  value   = other.value;
      }
    }
  }
  ~CSAItem(){
  }
};

enum ValueLength { NotApplicable=0, Fixed, Maximum };

const struct DicomValueRepresentation
{
  const char*  Name;
  const char*  Description;
  unsigned     length;		// in bytes
  ValueLength  lengthIs;
  bool         isAscii;
} DicomVR[] = {
  { "AE", "Application Entity",    16, Fixed, true },
  { "AS", "Age String",             4, Maximum, true }, // 0-9 D W M Y
  { "AT", "Attribute Tag",          4, Fixed, false },
  { "CS", "Code String",           16, Maximum, true },
  { "DA", "Date",                  10, Maximum, true }, // "yyyymmdd" or "yyyy.mm.dd"
  { "DS", "Decimal String",        16, Maximum, true },
  { "DT", "Date Time",             26, Maximum, true },
  { "FL", "Floating Point Single",  4, Fixed, false }, // IEEE 754
  { "FD", "Floating Point Double",  8, Fixed, false }, // IEEE 754
  { "IS", "Integer String",        12, Maximum, true },
  { "LO", "Long String",           64, Maximum, true },
  { "LT", "Long Text",          10240, Maximum, true },
  { "OB", "Other Byte String",      0, NotApplicable, false },
  { "OF", "Other Float String",     0, NotApplicable, false }, // IEEE 754 32Bit
  { "OW", "Other Word String",      0, NotApplicable, false },
  { "PN", "Person Name",           64, Maximum, true }, // family names^given names^middle names^name prefixes^name suffixes
  { "SH", "Short String",          16, Maximum, true },
  { "SL", "Signed Long",            4, Fixed, false },
  { "SQ", "Sequence of Items",      0, NotApplicable, false },
  { "SS", "Signed Short",           2, Fixed, false },
  { "ST", "Short Text",          1024, Maximum, true },
  { "TM", "Time",                  16, Maximum, true }, // hhmmss.ffffff
  { "UI", "Unique Identifier",     64, Maximum, true },
  { "UL", "Unsigned Long",          4, Fixed, false },
  { "UN", "Unknown",                0, NotApplicable, false },
  { "US", "U_Int16_T",              2, Fixed, false },
  { "UT", "Unlimited Text", 0xffffffff - 2, Maximum, true },
  { nullptr, nullptr, 0, NotApplicable, false }
};


bool
FetchCSA( char*            b,       // input data stream
	  int              len,     // input data length
	  vector<CSAItem>& itemList // return result as a list
	  )
{
  const int MinCSASize = 128;
  
  const char*  p            = b; // work on a copy of the start address
  unsigned int foundEntries = 0; // counter
  unsigned int expectedEntries;
  unsigned int subelements;
  unsigned int vl;
#ifdef DEBUGCSA
  unsigned int vm;
  u_int32_t    unused; // there is still s.th. left to investigate
#endif

  if( len <= MinCSASize ){
#ifdef DEBUGCSA
    cout << "FetchCSA(): not enough data." << endl;
#endif
    return false;
  }
#ifdef DEBUGCSA
  cout << "FetchCSA(): total data length=" << len << endl;
#endif

  if( strncmp ( p, "SV10", 4 ) == 0 ){
    // CSA2
    p += 4;
#ifdef DEBUGCSA
    cout << "\"CSA2\"" << endl;
#endif
  }
  else
    if( strncmp( p, "STEP;", 5) == 0 )
      {
	cout << "skipping \"STEP format\" while parsing CSA, still don't know enough about it ..." << endl;
	return false;
      }
    else
      cout << "presumably \"CSA1\", untested code follows..." << endl;

#ifdef DEBUGCSA
  unused = *(u_int32_t*) p;
  cout << unused << " (???)." << endl;
#endif
  p+=4;

  expectedEntries = *(u_int32_t*) p ;  p += 4;
#ifdef DEBUGCSA
  cout << "FetchCSA(): " << expectedEntries << " expected entries." << endl;
#endif

#ifdef DEBUGCSA
  unused = *(u_int32_t*) p;
  cout << unused << " (???)." << endl;
#endif
  p+=4;

  while( expectedEntries > foundEntries )
    {
      CSAItem entry;

      // the name is the first thing to get
      entry.name = p;
      p += 64;			// fixed length and null terminated
#ifdef DEBUGCSA
      cout << "____________________________" << foundEntries << "/" << expectedEntries << ":" << endl;
      cout << "CSA NAME=\"" << entry.name << "\"" << endl;
#endif

#ifdef DEBUGCSA
      vm = *(u_int32_t*) p ;
      cout << "VM=" << vm << endl;
#endif
      p += 4;

      // the second item should give the Value Representation
      const DicomValueRepresentation* dcmVR = DicomVR; // null terminated list
      // check against our list of known names
      while( dcmVR->Name )
	{
	  // find out which one
	  if( ( dcmVR->Name[0] == p[0] ) && ( dcmVR->Name[1] == p[1] ) )
	    {
	      // VR accepted
	      p += 4;
	      entry.vr = dcmVR->Name;
	      
	      ++foundEntries;
	      
#ifdef DEBUGCSA
	      cout << "VR=\"" << entry.vr << "\" " << foundEntries << "/" << expectedEntries << endl;
#endif
	      break;
	    }
	  // try next
	  dcmVR++;
	}
      if( ! dcmVR->Name )
	{
	  cerr << endl;
	  cerr << " --- No known VR entry found in assumed CSA, break after "
	       << p - b << " bytes ( out of " << len << " )" << endl;
	  break;
	}

#ifdef DEBUGCSA
      unused = *(u_int32_t*) p;
      cout << unused << " (???)." << endl;
#endif
       p+=4;

      subelements = *(u_int32_t*) p; p+=4;
#ifdef DEBUGCSA
      cout << "subelements=" << subelements << "" << endl;
#endif

#ifdef DEBUGCSA
      unused = *(u_int32_t*) p ;
      cout << unused << " (???)." << endl;
#endif
      p += 4;

      for( unsigned i = 0; i < subelements; i++ )
	{
	  vl = *(u_int32_t*)p ; p += 4;
#ifdef DEBUGCSA
	  int unused0 = *(u_int32_t*) p ;  p += 4;
	  int unused1 = *(u_int32_t*) p ;  p += 4;
	  int unused2 = *(u_int32_t*) p ;  p += 4;
	  cout << "vl=" << vl << "(" << unused0 << "," << unused1 << ","<< unused2
	       << ")" << endl;
#else
	  p += 4 + 4 + 4;
#endif
	  if( vl )
	    { // always assume ASCII here  - even 'FD' is written in ASCII, so ...
	      //                           - even 'UN' is used with a whole lot of ... things ...
	      entry.value.push_back(p);
#ifdef DEBUGCSA
	      cout << p << endl;
#endif
	      p += vl; // go ahead the value length
	    }
	  // step to the next 4byte boundary
	  while( (long)( p - b ) % 4 )
	    p++;
	}

#ifdef DEBUGCSA
      cout << "\"" << entry.name << "\", ";
      if( entry.value.size() )
	{
	  cout << entry.value.size() << " items: ";
	  for( auto& value : entry.value )
	    cout << "\"" << value << "\" ";
	  cout << endl;
	}
      else
	cout << " no items." << endl;
#endif

      itemList.push_back( entry );  // put it at the end of the list
    }

#ifdef DEBUGCSA
  cout <<  "FetchCSA(): expected " << expectedEntries << ", found " << foundEntries << endl;
#endif
  if( expectedEntries != foundEntries ){
    return false;
  }
  // still here ?
  return true;
}



//////////////////////////////////////////////////////////////////////////////
// i2o_dcm ///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

i2o_dcm::i2o_dcm() {/*empty*/ Init(); }

i2o_dcm::i2o_dcm( string fn ) {  Init(); Read( fn ); }

i2o_dcm::~i2o_dcm() {  Close(); }

void i2o_dcm::Init() // c++11 will get rid of this
{
  m_Cols = -1;
  m_TotalCols = -1;
  m_Rows = -1;
  m_TotalRows = -1;
  m_Slices = -1;
  m_FoVx = -1;
  m_FoVy = -1;

  m_StudyID = -1;
  m_SeriesNr = -1;
  m_AcquisitionNr = -1;
  m_InstanceNr = -1;
  m_TE = -1;
  m_TR = -1;
  m_TI = -1;
  m_FlipAngle = -1;
  m_VoxelWidth = -1;
  m_VoxelHeight = -1;
  m_SliceThickness = -1;
  m_SliceSpacing = -1;
  m_SliceLocation=-1;
  m_SamplesPerPixel = -1;

  m_filesize = 0;
  m_fullScanDone = false;
  m_readStoppedHere =nullptr;
  m_firstDcmTag=nullptr;
  m_mm = nullptr;
  m_isDicomdir = -1;
  m_payload = nullptr;
  m_payloadSize = 0;
  m_littleEndian = true;
  m_orderOfSlices = UndefinedOrder;
  m_echoNr = -1;
  m_averages = -1;
  m_bValue = -1;
  m_acquisitionMatrix = -1;

  m_dataTypeHashCode = 0;
  m_signedDataValues = false;
  m_bitsAllocated = -1;
  m_bitsStored = -1;
  m_highBit = -1;
  m_planarConfiguration = -1;
  m_numberOfFrames = -1;
  m_PhaseEncodingDirPositive = -1;
}

// in case of a failed method call, that might give an idea about what went wrong ( apart from stderr )
string i2o_dcm::Msg() { return m_msgTxt; }

string i2o_dcm::OrderOfSlices()
{
  if( !m_fullScanDone )
    FullScan();
  
  switch(m_orderOfSlices)
    {
    case UndefinedOrder:
      return "undefined slice order";
    case AscendingOrder:
      return "ascending";
    case DescendingOrder:
      return "descending";
    case InterleavedOrder:
      return "interleaved";
    case NoOrder:
      return "no slice order";
    }
  return "slice order ???"; // make the compiler happy
}

int
i2o_dcm::DataSize()
{
  if( m_payloadSize == 0 )
    if( !m_fullScanDone )
      FullScan();
  return m_payloadSize;
}

bool
i2o_dcm::SetData( void* newData, int totalByteSize )
{
  if( !m_fullScanDone )
    FullScan();
  
  if( totalByteSize == m_payloadSize )
    {
      memcpy( m_payload, newData, totalByteSize );
      return true;
    }
  else
    {
      m_msgTxt = "SetData(): sizes do not match but they have to, "
	+ to_string(totalByteSize) + "!=" + to_string(m_payloadSize) + ".";
      return false;
    }
}


#if __cplusplus >= 201103L
size_t
i2o_dcm::Datatype_hash_code()
{
  if( m_dataTypeHashCode == 0 )
    if( !m_fullScanDone )
      FullScan();
  return m_dataTypeHashCode;
}

string
i2o_dcm::Datatype_C_Name()
{
  size_t hash = Datatype_hash_code();

  if( hash == typeid(signed char).hash_code())
    return "signed char";
  else if( hash == typeid(unsigned char).hash_code())
    return "unsigned char";
  else if( hash == typeid(signed short).hash_code())
    return "signed short";
  else if( hash == typeid(unsigned short).hash_code())
    return "unsigned short";
  else if( hash == typeid(signed int).hash_code())
    return "signed int";
  else if( hash == typeid(unsigned int).hash_code())
    return "unsigned int";
  else if( hash == typeid(float).hash_code())
    return "float";
  else if( hash == typeid(complex<float>).hash_code())
    return "complex<float>";
  else
    return "unknown data type";
}
#endif


bool
i2o_dcm::GetRawData( void* dest )
{
  if( dest )
    {
      if( !m_fullScanDone )
	FullScan();
      if( m_payload ){
	memcpy( dest, m_payload, m_payloadSize );
	return true;
      }
    }
  return false;
}


void*
i2o_dcm::GetRawData( int& xdim, int& ydim, int& zdim, int& size )
{
  unsigned char* dest = nullptr;

  if( !m_fullScanDone )
    FullScan();

  if( m_payload && m_Rows && m_Cols && m_Slices )
    {
      xdim = m_Cols;
      ydim = m_Rows;
      zdim = m_Slices;

      dest = new (nothrow) unsigned char[ m_payloadSize ];

      if( dest == nullptr )
	{
	  m_msgTxt = "Could not allocate memory for " +
	    to_string(m_Cols) + "x" + to_string(m_Rows) + "x" + to_string(m_Slices) + " values.";
	  cerr << m_msgTxt << endl;
	}
      else
	{
	  memcpy( dest, m_payload, m_payloadSize );
	  size = m_payloadSize;
	}
    }
  else
    {
      m_msgTxt = "Invalid or missing dimension or missing data ( dim=\"" +
	to_string(m_Cols) + "x" + to_string(m_Rows) + "x" + to_string(m_Slices) + "\" ).";
      cerr << m_msgTxt << endl;
    }

  return dest;
}


void*
i2o_dcm::RawDataDirect( int& xdim, int& ydim, int& zdim, int& size )
{
  void* dest = nullptr;

  if( !m_fullScanDone )
    FullScan();

  if( m_payload && m_payloadSize && m_Rows && m_Cols && m_Slices )
    {
      xdim = m_Cols;
      ydim = m_Rows;
      zdim = m_Slices;
      size = m_payloadSize;
      dest = m_payload;		// <--
    }
  else
    {
      m_msgTxt = "Invalid or missing dimension or missing data ( dim=\""
	+ to_string(m_Cols) + "x" + to_string(m_Rows) + "x" + to_string(m_Slices) + "\" ).";
      cerr << m_msgTxt << endl;
    }

  return dest;
}


void*
i2o_dcm::GetSliceData( int slice, int& xdim, int& ydim, int& totalSliceSize, void* dest )
{
  if( !m_fullScanDone )
    FullScan();

  if( m_payload && m_Rows && m_Cols && m_Slices )
    {
      if( ( slice >= 0 ) && ( slice < m_Slices ) )
	{
	  int bytesPerValue = m_bitsAllocated/8;

	  xdim = m_Cols;
	  ydim = m_Rows;

	  if( dest == nullptr )
	    dest = new (nothrow) unsigned char[ xdim * ydim * bytesPerValue ];

	  if( dest == nullptr )
	    {
	      m_msgTxt = "Could not allocate memory for " +
		to_string(m_Cols) + "x" + to_string(m_Rows) + " values.";
	      cerr << m_msgTxt << endl;
	    }
	  else
	    {
	      int tilesPerRow = m_TotalCols / m_Cols;
	      int offset      = bytesPerValue*( (slice / tilesPerRow) * m_Rows * m_TotalCols +
						(slice % tilesPerRow) * m_Cols );

	      // type casts to suppress pointer arithmetic warnings
	      unsigned char* src = (unsigned char*)m_payload + offset;
	      unsigned char* tgt = (unsigned char*)dest;

	      for( int i = 0; i < m_Rows; i++)
		{
		  memcpy( tgt, src, m_Cols * bytesPerValue );
		  tgt += m_Cols * bytesPerValue;
		  src += m_TotalCols * bytesPerValue;
		}
	      totalSliceSize = m_Cols * m_Rows * bytesPerValue;
	    }
	}
      else
	{
	  m_msgTxt = "Invalid slice nr " +  to_string(slice) + ".";
	  cerr << m_msgTxt << endl;
	}
    }
  else
    {
      m_msgTxt = "Invalid or missing dimension or missing data ( dim=\""
	+ to_string(m_Cols) + "x" + to_string(m_Rows) + "x" + to_string(m_Slices) + "\" ).";
      cerr << m_msgTxt << endl;
    }

  return dest;
}

void*
i2o_dcm::Mosaic2Vol( int& xdim, int& ydim, int& zdim, int& destSize, void* dest )
{
  if( dest )
    {
      if( destSize != m_payloadSize )
	{
	  m_msgTxt = "Mosaic2Vol: data size values do not match, " +
	    to_string( destSize ) + " != " + to_string( m_payloadSize ) + ".";
	  cerr << m_msgTxt << endl;
	  return nullptr;
	}
    }
  else
    {
      destSize = m_payloadSize;

      dest = new (nothrow) unsigned char[ destSize ];
      if( dest == nullptr )
	{
	  m_msgTxt = "Could not allocate memory for "
	    + to_string(m_Cols) + "x" + to_string(m_Rows) + " values.";
	  cerr << m_msgTxt << endl;
	  return nullptr;
	}
    }

  xdim = m_Cols;
  ydim = m_Rows;

  unsigned char* tgt = (unsigned char*)dest;
  for( zdim = 0; zdim < m_Slices; zdim++ )
    {
      int totalSliceSize;

      GetSliceData(zdim, xdim, ydim, totalSliceSize, tgt);
      tgt += totalSliceSize;
    }

  return dest;
}


string
i2o_dcm::Info()
{
  string info;

  if( m_mm )
    {
      if( IsDicomDir() )
	{
	  info = "DICOMDIR file";
	}
      else
	{
	  info = "regular dicom file";
	}
    }
  else
    info = "NA";
  return info;
}


//
// try to read an integer value from a dicom data field
//
int
GetValue(char* b, char* vr, int len)
{
  if( vr[0] == '\0' ) // implicit VR, might be binary or ascii
    {
      bool  isAscii = true;
      short  digits = 0;
      
      for( int i=0; i < len && b[i] != ' '; i++ )
	if( b[i] < '0' || b[i] > '9' )
	  {
	    isAscii = false;
	    break;
	  }
	else
	  digits++;
      
      if( isAscii )
	if( digits )
	  return atoi( b );
	else
	  {
	    cerr << "INTERNAL ERROR ! GetValue() empty string: len=" << len << ", VR=\"" << vr << "\"" << endl;
	    exit(1);
	  }
      else
	{
	  if( len == 1 )
	    return *(int8_t*)b;
	  if( len == 2 )
	    return *(int16_t*)b;
	  if (len == 4)
	    return *(int32_t*)b;
	}
    }
  else // explicit VR, use it
    { 
      if( vr[0] == 'S' && vr[1] == 'L' ) // signed long
	return *(int32_t*)b;
      if( vr[0] == 'S' && vr[1] == 'S' ) // signed short
	return *(int16_t*)b;
      if( vr[0] == 'U' && vr[1] == 'L' ) // unsigned long
	return *(u_int32_t*)b;
      if( vr[0] == 'U' && vr[1] == 'S' ) // unsigned short
	return *(u_int16_t*)b;
      
      return atoi( b );
    }
  cerr << "INTERNAL ERROR ! GetValue() unhandled case: len=" << len << ", VR=\"" << vr << "\"" << endl;
  exit(1);
}

void i2o_dcm::Grp08(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x08:
      m_ImageType.assign( b, len );
      trim( m_ImageType );
      // might be e.g.: "ORIGINAL\PRIMARY\M\ND\NORM"
      //           or : "ORIGINAL\PRIMARY\M\ND\NORM\MOSAIC"
#ifdef DEBUG
      cout << "len=" << len << ", VR=" << VR << ", ImageType=" << m_ImageType << endl;
#endif
      break;
      
    case 0x20:
      m_StudyDate.assign( b, len );
      trim( m_StudyDate );
      break;
      
    case 0x21:
      m_SeriesDate.assign( b, len );
      trim( m_SeriesDate );
      break;
      
    case 0x22:
      m_AcquisitionDate.assign( b, len );
      trim( m_SeriesDate );
      break;
      
    case 0x23:			// formerly known as "image date"
      m_ContentDate.assign( b, len );
      trim( m_ContentDate );
      break;
      
    case 0x30:
      m_StudyTime.assign( b, len ); // hhmmss.ffffff
      trim(m_StudyTime);
      break;
      
    case 0x31:
      m_SeriesTime.assign( b, len ); // hhmmss.ffffff
      trim( m_SeriesTime);
      break;
      
    case 0x32:				  // time at which the raw image data started to be aquired
      m_AcquisitionTime.assign( b, len ); // hhmmss.ffffff
      trim( m_AcquisitionTime );
      break;
      
    case 0x33:			      // image pixel data creation
      m_ContentTime.assign( b, len ); // hhmmss.ffffff
      trim( m_ContentTime );
      break;
      
    case 0x50:			// AccessionNumber
      break;
      
    case 0x60:
      m_Modality.assign( b, len );
      trim( m_Modality );
      break;
      
    case 0x1030:
      m_StudyDescription.assign( b, len );
      trim( m_StudyDescription );
      break;
      
    case 0x103e:
      m_SeriesDescription.assign( b, len );
      trim( m_SeriesDescription );
      break;
      
    case 0x1090:
      m_ManufacturerModelName.assign( b, len );
      trim( m_ManufacturerModelName );
      break;
      
    default:
      break;
    }
}

void i2o_dcm::Grp10(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x10:
      m_PatientName.assign( b, len );
      trim( m_PatientName );
      break;
      
    case 0x20:
      m_PatientID.assign( b, len );
      trim( m_PatientID );
      break;
      
    case 0x30:
      m_PatientBirthday.assign( b, len );
      trim( m_PatientBirthday );
      break;
      
    case 0x40:
      m_Gender.assign( b, len );
      trim( m_Gender );
      break;
      
    case 0x1010:
      m_PatientAge.assign( b, len );
      trim( m_PatientAge );
      break;
      
    case 0x1030:
      m_PatientWeight.assign( b, len );
      trim( m_PatientWeight );
      break;
      
    case 0x4000:
      m_PatientComment.assign( b, len );
      trim( m_PatientComment );
      break;
      
    default:
      break;
    }
}


void i2o_dcm::Grp18(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x20: // scanning sequence
      break;
      
    case 0x21: // scanning variant
      break;
      
    case 0x22: // scan options
      break;
      
    case 0x23: // MR Acquisitiontype
      break;
      
    case 0x24: // Sequence name
      break;
      
    case 0x25: // Angio flag
      break;
      
    case 0x50: //  Slice Thickness
      m_SliceThickness = atof(b);
      break;
      
    case 0x80: //  Repetition time,
      m_TR = atof(b);
      break;
      
    case 0x81: //  Echo time
      m_TE = atof(b);
      break;
      
    case 0x82: //  Inversion time
      m_TI = atof(b);
      break;
      
    case 0x83: //  Number of averages
      break;
      
    case 0x84: //  Imaging frequency
      break;
      
    case 0x85: //  Imaged Nucleus
      break;
      
    case 0x86: //  Imaged Nucleus
      break;
      
    case 0x87: //  Magnetic field strength
      break;
      
    case 0x88: // Spacing Between Slices
      m_SliceSpacing = atof(b);
      break;
      
    case 0x89: //  Number of phase encoding steps
      // int nrOfPhaseEncodingSteps = (u_int16_t*)b;
      break;
      
    case 0x91: //  echo train length
      break;
      
    case 0x93: //  percent sampling
      break;
      
    case 0x94: //  percent phase field of view
      break;
      
    case 0x95: //  pixel bandwidth
      break;
      
    case 0x1000: // device serial number
      // DeviceSerial.assign(b,len);
      break;
      
    case 0x1020:
      m_SWVersion.assign(b,len);
      trim( m_SWVersion );
      break;
      
    case 0x1030:
      m_ProtocolName.assign(b,len);
      trim( m_ProtocolName);
      break;
      
    case 0x1251:
      //transmitCoilName.assign(b,len);
      break;
      
    case 0x1310:
      // acquisionMatrix, multiplicity=4,
      // e.g.: 4000 0000 0000 4000 for 64\0\0\64 if row axis was the freq encoding axis
      // e.g.: 0000 4000 4000 0000 for 0\64\64\0 if row axis was the phase encoding axis
      // "frequency rows\frequency columns\phase rows\phase columns"
      // usefull for mosaic size detection if everything else fails
      {

	u_int16_t* p = (u_int16_t*)b;
	u_int16_t freqEncRows,freqEncCols,phaseEncRows,phaseEncCols;

	freqEncRows  = p[0];
	freqEncCols  = p[1];
	phaseEncRows = p[2];
	phaseEncCols = p[3];
#ifdef DEBUG
	cout << "acquisionMatrix: " << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << endl;
#endif
	int rows = freqEncRows ? freqEncRows:phaseEncRows;
	int cols = freqEncCols ? freqEncCols:phaseEncCols;

	if( m_Rows == -1 )
	  m_Rows = rows;
	if( m_Cols == -1 )
	  m_Cols = cols;
      }
      break;
      
    case 0x1312: // InPlanePhaseEncodingDirection
      m_InPlanePhaseEncodingDir.assign( b, len );
      trim( m_InPlanePhaseEncodingDir );
      // "COL" phase encoded in columns
      // "ROW" phase encoded in rows
      break;
      
    case 0x1314:
      m_FlipAngle = GetValue(b, VR, len);
      break;
      
    case 0x1315: // variable flip angle flag
      break;
      
    case 0x1316: // SAR
      // m_Sar = atof(b);
      break;
      
    case 0x1318: // dBdt
      //m_dBdt = atof(b);
      break;
      
    case 0x5100:
      // patient's position, eg "HFS" = Head-First/Supine
      m_PatientPosition.assign( b, len );
      trim( m_PatientPosition );
      break;
      
    }
}

void i2o_dcm::Grp19(char* b, u_int16_t element, int32_t len, char VR[3])
{
#ifdef DEBUG
  string s_b;
#endif
  switch(element)
    {
    case 0x10:			// Private Creator
      // "SIEMENS MR HEADER"
      break;

    case 0x1008:		// CSA Image Header Type
      // "IMAGE NUM 4"
      break;

    case 0x1009:		// CSA Image Header Version
      // LO "1.0" 
      break;
      
    case 0x100a:		// nrOfSlicesInMosaic
      m_Slices = *(u_int16_t*) b;
      break;
      
    case 0x100b:	 // slice measurement duration
      // double sliceMeasurementDuration = atof(b);
      break;
      
    case 0x100c:		// B Value 
      m_bValue = GetValue(b, VR, len);
      break;
      
    case 0x100d:		// diffusion directionality
      // CS "DIRECTIONAL" "NONE"
      break;
      
    case 0x100e:		// diffusion grad direction
      // FD |0.114766\-0.264845\0.918128|
      {
      	double* p = (double*)b;
	int     n = len/sizeof(double);

	m_DiffusionGradDirection.resize( n );
	
	for( int i = 0; i < n; i ++ )
	  m_DiffusionGradDirection[i] = p[i];

#ifdef DEBUG	
	cout << "DBG: Diffusion Direction \"" << m_DiffusionGradDirection[0];
	for( int i = 1; i < n; i ++ )
	  cout << "," << m_DiffusionGradDirection[i];
	cout << "\"" << endl;
#endif
      }
      break;
      
    case 0x100f:		// gradient mode
#ifdef DEBUG
      s_b.assign(b,len);
      trim(s_b);
      cout << "DBG: Gradient Mode \"" << s_b << "\"" << endl;
#endif
      break;

    case 0x1012:		// TablePositionOrigin
      // SL 0\0\-1270
      break;

    case 0x1013:		// ImaAbsTablePosition
      // SL 0\0\-1270 
      break;

    case 0x1014:		// ImaRelTablePosition
      // IS [0\0\0]
      break;

    case 0x1015:		// SlicePosition_PCS
      // FD -877.59036135999997\-890.84337330000005\-30.527109150000001
      break;

    case 0x1016:		// TimeAfterStart
      // DS [14.99]
      break;

    case 0x1017:		// SliceResolution
      // DS [1]
      break;

    case 0x1018:		// RealDwellTime
      // IS [2000] 
      break;

    case 0x1027:		// diff/B_matrix
      // FD |13\-30\104\73\-254\901|
      {
#ifdef DEBUG
	double* p = (double*)b;
	int     n = len/sizeof(double);
	
	cout << "DBG: diff/B_matrix \"" << p[0];
	for( int i = 1; i < n; i ++ )
	  cout << "," << p[i];
	cout << "\"" << endl;
#endif
      }
      
      break;
      
    case 0x1028:		// Bandwidth per pixel phase encode
#ifdef DEBUG
      {
	double* p = (double*)b;
	cout << "DBG: Bandwidth per pixel \"" << p << "\"" << endl;
      }
#endif
      break;
      
    case 0x1029:      // MosaicRefAcqTimes starting with 0.0
      // FD, len/(sizeof(double) == nrOfImagesInMosaic
      {
	int     n = len/sizeof(double);
	double* p = (double*)b;

	m_AcquisitionRefTimes.resize( n );
	
	for( int i=0; i < n; i++ )
	  m_AcquisitionRefTimes[i] = *p++;
	
#ifdef DEBUG
	cout << "DBG: MosaicRefAcqTimes, n==" << n << endl;
	for(int i = 0; i < n; i++)
	  cout << "MosaicRefAcqTimes[" << i << "]: " << m_AcquisitionRefTimes[i] << endl;
#endif
      }
      break;
      
    default:
      break;
    }
}

void i2o_dcm::Grp20(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x0D:
      m_StudyInstanceUID.assign( b,len );
      trim( m_StudyInstanceUID );
      break;

    case 0x0E:
      m_SeriesInstanceUID.assign( b,len );
      trim( m_SeriesInstanceUID );
      break;

    case 0x10:
      m_StudyID = GetValue(b, VR, len);
      break;

    case 0x11:
      m_SeriesNr = GetValue(b, VR, len);
      break;

    case 0x12:
      m_AcquisitionNr = GetValue(b, VR, len);
      break;

    case 0x13:
      m_InstanceNr = GetValue(b, VR, len); // aka image number
      break;

    case 0x32: // ImagePositionPatient
      // "The Image Position (0020,0032) specifies the x, y, and z coordinates
      // of the upper left hand corner of the image; it is the center of the
      // first voxel transmitted."
      double posX, posY, posZ;

      if( sscanf( b, "%lf\\%lf\\%lf", &posX, &posY, &posZ ) == 3 )
	{
	  m_ImgPosition.resize(3);

	  m_ImgPosition[0] = posX;
	  m_ImgPosition[1] = posY;
	  m_ImgPosition[2] = posZ;
	}
      else
	cerr << "Error reading dicom tag \"Image Position Patient\" (0020,0032)." << endl;
      break;

    case 0x37: // ImageOrientationPatient
      // from the standard:
      // "Image Orientation (0020,0037) specifies the
      //  direction cosines of the first row and the first column with respect
      //  to the patient. These Attributes shall be provide as a pair.
      //  Row value for the x, y, and z axes respectively followed by the Column
      //  value for the x, y, and z axes respectively.
      //  The direction of the axes is defined fully by the patient's orientation.
      //  The x-axis is increasing to the left hand side of the patient.
      //  The y-axis is increasing to the posterior side of the patient.
      //  The z-axis is increasing toward the head of the patient.
      //  The patient based coordinate system is a right handed system, i.e. the
      //  vector cross product of a unit vector along the positive x-axis and a
      //  unit vector along the positive y-axis is equal to a unit vector along
      //  the positive z-axis.
      //  Note: If a patient lies parallel to the ground, face-up on the table,
      //  with his feet-to-head direction same as the front-to-back direction of
      //  the imaging equipment, the direction of the axes of this patient based
      //  coordinate system and the equipment based coordinate system in previous
      //  versions of this Standard will coincide."
      double rowX,rowY, rowZ;
      double colX,colY, colZ;

      sscanf(b,"%lf\\%lf\\%lf\\%lf\\%lf\\%lf", &rowX, &rowY, &rowZ, &colX, &colY, &colZ);

      m_ImgOrientRow.resize(3);
      m_ImgOrientCol.resize(3);
      m_ImgOrientNormal.resize(3);

      m_ImgOrientRow[0] = rowX;
      m_ImgOrientRow[1] = rowY;
      m_ImgOrientRow[2] = rowZ;

      m_ImgOrientCol[0] = colX;
      m_ImgOrientCol[1] = colY;
      m_ImgOrientCol[2] = colZ;

      // once we have X and Y we can have Z as the cross product:
      m_ImgOrientNormal[0] = rowY*colZ - rowZ*colY;
      m_ImgOrientNormal[1] = rowZ*colX - rowX*colZ;
      m_ImgOrientNormal[2] = rowX*colY - rowY*colX;

#ifdef DEBUG
      cout << b << endl;
      cout << "ImageOriPat,row:" << m_ImgOrientRow[0] << "," << m_ImgOrientRow[1] << "," << m_ImgOrientRow[2] << endl
	   << "ImageOriPat,col:" << m_ImgOrientCol[0] << "," << m_ImgOrientCol[1] << "," << m_ImgOrientCol[2] << endl
	   << "-> [row]X[col] :" << endl
	   << m_ImgOrientRow[1]*m_ImgOrientCol[2] - m_ImgOrientRow[2]*m_ImgOrientCol[1] << endl
	   << m_ImgOrientRow[2]*m_ImgOrientCol[0] - m_ImgOrientRow[0]*m_ImgOrientCol[2] << endl
	   << m_ImgOrientRow[0]*m_ImgOrientCol[1] - m_ImgOrientRow[1]*m_ImgOrientCol[0] << endl;
#endif
      break;

    case 0x1002: // Images in Acquisition (usually missing)
      // cout << "Images in Acquisition " << endl;
      break;

    case 0x1004: // Acquisitions in Study (usually missing)
      // cout << "Acquisitions in Study " << endl;
      break;

    case 0x1041: // SliceLocation
      m_SliceLocation = atof(b);
      break;
    }
}

void i2o_dcm::Grp28( char* b, u_int16_t element, int32_t len, char VR[3] )
{
  switch(element)
    {
    case 0x02:
      m_SamplesPerPixel = GetValue(b, VR,len);
#ifdef DEBUG
      if ( m_SamplesPerPixel > 1 )
	cout << "\""
	     << m_filename
	     << "\": samples per pixel == " << m_SamplesPerPixel << ", UNTESTED. Please report !"
	     << endl;
#endif
      break;
      
    case 0x04: // Photometric Interpretation
      {
#ifdef DEBUG
	// TODOTODOTODOTODOTODOTODOTODOTODOTODOTODOTODOTODOTODO
	char oneAndOnly[]="MONOCHROME2"; // i.e. greyscale
	// "RGB", "MONOCHROME1" (0==white), "YBR_FULL","YBR_FULL_422" (jpeg)

	if( strncmp( b, oneAndOnly, strlen(oneAndOnly) ) )
	  {
	    string s(b,len);
	    cout  << "\""
		  << m_filename
		  << "\": photometric Interpretation (\""
		  << s << "\"). UNTESTED. Please report !" << endl;
	  }
#endif
      }
      break;
      
    case 0x06:
      // interlaced(0)==default, separated(1)
      // (samplesPerPixel must be > 1)
      m_planarConfiguration = GetValue(b, VR, len);
#ifdef DEBUG
      cout  << "\""
	    << m_filename
	    << "\": Planar configuration=" << m_planarConfiguration <<", UNTESTED. Please report !"
	    << endl;
#endif
      break;
      
    case 0x08:
      m_numberOfFrames = GetValue(b, VR, len);
#ifdef DEBUG
      cout  << "\""
	    << m_filename
	    << "\": NumberOfFrames=" << m_numberOfFrames <<", UNTESTED. Please report !"
	    << endl;
#endif
      break;
      
    case 0x10:
      m_TotalRows = GetValue(b, VR, len);
      break;
      
    case 0x11:
      m_TotalCols = GetValue(b, VR, len);
      break;
      
    case 0x30: // pixel spacing
      // Physical distance in the patient between the center of each pixel, specified by a
      // numeric pair - adjacent row spacing (delimiter) adjacent column spacing in mm
      sscanf(b,"%lf\\%lf", &m_VoxelHeight, &m_VoxelWidth); // row spacing, column spacing
      break;
      
    case 0x100:
      // usually 16
      m_bitsAllocated = GetValue(b, VR, len);
      break;
      
    case 0x101:
      // usually 12
      m_bitsStored = GetValue(b, VR, len);
      break;
      
    case 0x102:
      // usually 11
      m_highBit = GetValue(b, VR, len);
      break;
      
    case 0x103:
      // unsigned(0), signed(1)
      m_signedDataValues = GetValue(b, VR, len);
      break;

    default:
      break;
    }
}

void i2o_dcm::Grp29(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x1008: // CSADataType
      // e.g. "SPEC NUM 4"
      break;
      
    case 0x1009: // CSAImageHeaderVersion
      // e.g. "syngo MR B17"
      break;
      
    case 0x1010: // CSAImageHeaderInfo
    case 0x1110: // CSAImageHeaderInfo
      {
	vector< CSAItem > csaList;
#ifdef DEBUG
	cout << "0x10x10 parsing ..."  << endl;
#endif
	if( FetchCSA((char*)b, len, csaList) )
	  {
	    for( auto& csa : csaList ){
	      if( ! csa.value.empty() ) // FILTER OUT EMPTY ENTRIES
		{
#ifdef DEBUG
		  //#if 1
		  cout << ": (29:1010) '" << csa.name
		       << "', itemCnt=" << csa.value.size()
		       << ", VR='" << csa.vr
		       << "', '";
		  for( auto& j : csa.value )
		    cout << j << ";";
		  cout << "'" << endl;
#endif
		  if( csa.name == "EchoLinePosition" )         {}
		  else if( csa.name == "ImagePositionPatient" )
		    {
#if 0
		      cout << "CSA ImagePositionPatient, " << csa.value.size() << " items :"<< endl;
		      for ( int n = 0; n < i->itemCnt; n++)
			cout << "CSA ImagePositionPatient" << csa.value.c_str() << endl;
#endif
		    }
		  else if( csa.name == "ImageOrientationPatient" )
		    {
#if 0
		      cout << "CSA ImageOrientationPatient, " << csa.itemCnt << " items :"<< endl;
		      for ( int n = 0; n < csa.itemCnt; n++)
			cout << "CSA ImageOrientationPatient" << csa.value[n].c_str() << endl;
#endif
		    }
		  else if( csa.name == "FlipAngle" )
		    {
		      if( m_FlipAngle <= 0 )
			m_FlipAngle = atof( csa.value.front().c_str() );
		      else
			if( m_FlipAngle != atof(csa.value.front().c_str()))
			  cerr << "\"" << m_filename << "\": Header values for FlipAngle values differ ! "
			       << m_FlipAngle << " != " << atof(csa.value.front().c_str()) << endl;
		    }
		  else if( csa.name == "SliceThickness" )
		    {
		      if( m_SliceThickness <= 0 )
			m_SliceThickness = atof(csa.value.front().c_str());
		      else
			if( m_SliceThickness != atof(csa.value.front().c_str()))
			  cerr << "\"" << m_filename << "\": Header values for SliceThickness values differ ! " 
			       << m_SliceThickness << " != " << atof(csa.value.front().c_str()) << endl;
		    }
		  else if( csa.name == "DataPointRows" )
		    {
		      int dataPointRows = atoi(csa.value.front().c_str());
		      if( m_TotalRows != -1 && m_TotalRows != dataPointRows )
			{
#ifdef DEBUG
			  cout << "\"DataPointRows\": Resetting m_TotalRows from " << m_TotalRows
			       << "to " << dataPointRows << endl;
#endif
			  m_TotalRows = dataPointRows;
			}
		    }
		  else if( csa.name == "DataPointColumns" )
		    {
		      int dataPointCols = atoi(csa.value.front().c_str());
		      if( m_TotalCols != -1 && m_TotalCols != dataPointCols )
			{
#ifdef DEBUG
			  cout << "\"DataPointCols\": Resetting m_TotalCols from " << m_TotalCols
			       << "to " << dataPointCols << endl;
#endif
			  m_TotalCols = dataPointCols;
			}
		      m_TotalCols = dataPointCols;
		    }
		  else if( csa.name == "DataRepresentation" )
		    {
		      if( csa.value.front() == "COMPLEX" )
			m_dataTypeHashCode = typeid(complex<float>).hash_code();
		    }
		  else if( csa.name == "Rows" )
		    {
		      int rows = atoi(csa.value.front().c_str());
		      if( m_Rows != -1 && m_Rows != rows )
			{
#ifdef DEBUG
			  cout << "\"Rows\": Resetting m_Rows from " << m_Rows
			       << "to " << rows << endl;
#endif
			  m_Rows = rows;
			}
		    }
		  else if( csa.name == "Columns" )
		    {
		      int cols = atoi(csa.value.front().c_str());
		      if( m_Cols != -1 && m_Cols != cols )
			{
#ifdef DEBUG
			  cout << "\"Cols\": Resetting m_Cols from " << m_Cols
			       << "to " << cols << endl;
#endif
			  m_Cols = cols;
			}
		    }
		  else if( csa.name == "SliceLocation" ) {}
		  else if( csa.name == "PixelSpacing" )  {}
		  else if( csa.name == "PixelBandwidth" ){}
		  else if( csa.name == "SpectroscopyAcquisitionPhaseColumns" ){}
		  else if( csa.name == "SpectroscopyAcquisitionPhaseRows" )   {}
		  else if( csa.name == "SpectroscopyAcquisitionOut-of-planePhaseSteps" ){}
		  else if( csa.name == "SpectroscopyAcquisitionDataColumns" ) {}
		  else if( csa.name == "SignalDomainColumns" )  {}
		  else if( csa.name == "EchoColumnPosition" )   {}
		  else if( csa.name == "EchoPartitionPosition" ){}
		  else if( csa.name == "UsedChannelMask" )      {}
		  else if( csa.name == "Actual3DImaPartNumber" ){}
		  else if( csa.name == "ICE_Dims" )             {}
		  else if( csa.name == "Filter1" )              {}
		  else if( csa.name == "Filter2" )              {}
		  else if( csa.name == "ProtocolSliceNumber" )  {}
		  else if( csa.name == "RealDwellTime" )        {}
		  else if( csa.name == "PixelFile" )            {}
		  else if( csa.name == "PixelFileName" )        {}
		  else if( csa.name == "B_value" )
		    {
#ifdef DEBUG
		      cout << "B_value: " << csa.value.front() << endl;
#endif
		      m_bValue = atof(csa.value.front().c_str());
		    }
		  else if( csa.name == "SliceMeasurementDuration" )
		    {
		      //double sliceMeasurementDuration = atof(csa.front().c_str());
		      //cout << "SliceMeasurementDuration=" << sliceMeasurementDuration << endl;
		    }
		  else if( csa.name == "SequenceMask" )         {}
		  else if( csa.name == "AcquisitionMatrixText" )
		    {
		      if( m_acquisitionMatrix <= 0){
			m_acquisitionMatrix = atoi( csa.value.front().c_str() );
		      }
		    }
		  else if( csa.name == "MeasuredFourierLines" ) {}
		  else if( csa.name == "FlowEncodingDirection" ){}
		  else if( csa.name == "FlowVenc" )             {}
		  else if( csa.name == "UsedChannelString" )    {}
		  else if( csa.name == "PhaseEncodingDirection" )
		    {
#ifdef DEBUG
		      cout << "!!! " << csa.name << " " << csa.value.size() << " " << csa.value.front() << endl;
#endif
		    }
		  else if( csa.name == "PhaseEncodingDirectionPositive" )
		    {
#ifdef DEBUG
		      cout << "!!! " << csa.name << " " << csa.value.size() << " " << csa.value.front() << endl;
		      int n = 0;
		      for( auto v : csa.value )
			cout << "PhaseEncodingDirectionPositive[" << n++ << "]: " << v << endl;
#endif
		      m_PhaseEncodingDirPositive = atoi(csa.value.front().c_str());
		    }
		  else if( csa.name == "DiffusionDirectionality" )
		    {
#ifdef DEBUG
		      cout << "!!! " << csa.name << endl;
#endif
		    }
		  else if( csa.name == "DiffusionGradientDirection" )
		    {
#ifdef DEBUG
		      cout << "!!! " << csa.name << endl;
#endif
		    }
		  else if( csa.name == "TimeAfterStart" )     {}
		  else if( csa.name == "ImaAbsTablePosition" ){}
		  else if( csa.name == "NonPlanarImage" )     {}
		  else if( csa.name == "SlicePosition_PCS" )  {}
		  else if( csa.name == "MultistepIndex" )     {}
		  else if( csa.name == "ImaRelTablePosition" ){}
		  else if( csa.name == "ImaCoilString" )      {}
		  else if( csa.name == "RFSWDDataType" )      {}
		  else if( csa.name == "GSWDDataType" )       {}
		  else if( csa.name == "ImaPATModeText" )     {}
		  else if( csa.name == "BandwidthPerPixelPhaseEncode" ){}
		  else if( csa.name == "MRDiffusion" )
		    {
#if 0
		      cout << "MRDiffusion, items=" << csa.itemCnt << ":" << endl;
		      for(int n=0; n < csa.itemCnt; n++ )
			cout << "MRDIffusion[" << n << "]: " << csa.value[n] << endl;
#endif
		    }
		  else if( csa.name == "ImageType4MF" )       {}
		  else if( csa.name == "ImageHistory" )       {}
		  else if( csa.name == "MosaicRefAcqTimes" )
		    {
#ifdef DEBUG
		      // alternate candidate for these values, see (0x10,0x1029)
		      cout << "MosaicRefAcqTimes, items=" << csa.value.size() << ":" << endl;
		      int n = 0;
		      for( auto v : csa.value )
			cout << "MosaicRefAcqTimes[" << n++ << "]: " << v << endl;
#endif
		    }
		  else if( csa.name == "NumberOfImagesInMosaic" )
		    {
		      if( csa.value.size() )
			{
			  m_Slices = atoi(csa.value.front().c_str());
#ifdef DEBUG
			  cout << "CSA: NumberOfSlicesInMosaic = " << m_Slices << endl;
#endif
			}
		    }
		  else if( csa.name == "SliceNormalVector" )
		    {
		      if( csa.value.size() == 3 )
			{
			  m_SliceNormal.resize(3);
			  auto v = csa.value.cbegin();
			  m_SliceNormal[0] = atof( (v++)->c_str() );
			  m_SliceNormal[1] = atof( (v++)->c_str() );
			  m_SliceNormal[2] = atof(  v->c_str() );
#ifdef DEBUG
			  cout << "Image normal=="
			       << m_SliceNormal[0] << ","
			       << m_SliceNormal[1] << ","
			       << m_SliceNormal[2] << endl;
#endif			
			}
		    }
		}
	    }
	  }
      }
      break;

    case 0x1018: // CSASeriesHeaderType
      // CS [MR]
      break;
      
    case 0x1019: // CSASeriesHeaderVersion
      // LO [20101201]
      break;
      
    case 0x1020: // CSASeriesHeaderInfo
      {
	vector<CSAItem> csaList;

	if( FetchCSA( (char*)b, len, csaList) )
	  {
	    for( auto& csa : csaList )
	      {
		if( ! csa.value.empty() ) // FILTER OUT EMPTY ENTRIES
		{
#ifdef DEBUG
		  //#if 1
		  cout << ": (29:1020): '"
		       << csa.name
		       << "', vm="
		       << csa.value.size()
		       << ", VR='" << csa.vr
		       << "', '";
		  for( auto& value : csa.value )
		    cout << value << ";";
		  cout << "'" << endl;
#endif
		  if( csa.name == "UsedPatientWeight" )	            {}
		  else if( csa.name == "ImagePositionPatient" )
		    {
#if 0
		    cout << "CSA ImagePositionPatient, " << csa.itemCnt << " items:"<< endl;
		    for ( int n = 0; n < csa.itemCnt; n++)
		      cout << "CSA ImagePositionPatient" << csa.value[n].c_str() << endl;
#endif
		  }
		  else if( csa.name == "ImageOrientationPatient" )
		    {
#if 0
		    cout << "CSA ImageOrientationPatient, " << csa.itemCnt << " items:"<< endl;
		    for ( int n = 0; n < csa.itemCnt; n++)
		      cout << "CSA ImageOrientationPatient" << csa.value[n].c_str() << endl;
#endif
		  }
		  else if( csa.name == "NumberOfPrescans" )	    {}
		  else if( csa.name == "TransmitterCalibration" )    {}
		  else if( csa.name == "PhaseGradientAmplitude" )    {}
		  else if( csa.name == "ReadoutGradientAmplitude" )  {}
		  else if( csa.name == "SelectionGradientAmplitude" ){}
		  else if( csa.name == "GradientDelayTime" )	    {}
		  else if( csa.name == "RfWatchdogMask" )	    {}
		  else if( csa.name == "RfPowerErrorIndicator" )	    {}
		  else if( csa.name == "SarWholeBody" )		    {}
		  else if( csa.name == "Sed" )		            {}
		  else if( csa.name == "SequenceFileOwner" )	    {}
		  else if( csa.name == "Stim_mon_mode" )		    {}
		  else if( csa.name == "Operation_mode_flag" )	    {}
		  else if( csa.name == "dBdt_max" )		    {}
		  else if( csa.name == "t_puls_max" )		    {}
		  else if( csa.name == "dBdt_thresh" )		    {}
		  else if( csa.name == "dBdt_limit" )		    {}
		  else if( csa.name == "SW_korr_faktor" )	    {}
		  else if( csa.name == "Stim_max_online" )	    {}
		  else if( csa.name == "Stim_max_ges_norm_online" )  {}
		  else if( csa.name == "Stim_lim" )		    {}
		  else if( csa.name == "Stim_faktor" )		    {}
		  else if( csa.name == "CoilForGradient" )	    {}
		  else if( csa.name == "CoilTuningReflection" )	    {}
		  else if( csa.name == "CoilId" )		    {}
		  else if( csa.name == "MiscSequenceParam" )	    {}
		  else if( csa.name == "MrProtocolVersion" )	    {}
		  else if( csa.name == "MrProtocol" )		    {}
		  else if( csa.name == "DataFileName" )		    {}
		  else if( csa.name == "RepresentativeImage" )	    {}
		  else if( csa.name == "PositivePCSDirections" )	    {}
		  else if( csa.name == "RelTablePosition" )	    {}
		  else if( csa.name == "ReadoutOS" )		    {}
		  else if( csa.name == "LongModelName" )		    {}
		  else if( csa.name == "SliceArrayConcatenations" )  {}
		  else if( csa.name == "SliceResolution" )	    {}
		  else if( csa.name == "MrEvaProtocol" )		    {}
		  else if( csa.name == "AbsTablePosition" )	    {}
		  else if( csa.name == "AutoAlignMatrix" )	    {}
		  else if( csa.name == "MeasurementIndex" )	    {}
		  else if( csa.name == "CoilString" )		    {}
		  else if( csa.name == "PATModeText" )	            {}
		  else if( csa.name == "PatReinPattern" )	    {}
		  else if( csa.name == "MosaicRefAcqTimes" )
		    {
#ifdef DEBUG
		      cout << "MosaicRefAcqTimes in CSA, unused/skipped." << endl;
#endif
		    }
		  else if( csa.name == "MrPhoenixProtocol" ) // the <XProtocol>
		    {
		      size_t found = csa.value.front().find("sSliceArray.ucMode");
		      if( found != string::npos )
			{
			  const char* p = csa.value.front().c_str() + found;

			  int n, sliceOrdering;

			  // might be hex or decimal
			  if( ( n = sscanf( p, "sSliceArray.ucMode = %x", &sliceOrdering ) ) == 0 )
			    n = sscanf( p, "sSliceArray.ucMode = %d", &sliceOrdering );

			  if( n > 0 )
			    {
			      /* 1 == Ascending
				 2 == Descending
				 4 == Interleaved
			      */
			      switch(sliceOrdering)
				{
				case 1: m_orderOfSlices = AscendingOrder;   break;
				case 2: m_orderOfSlices = DescendingOrder;  break;
				case 4: m_orderOfSlices = InterleavedOrder; break;
				default:
				  cerr << "\"" << m_filename
				       << "\": sSliceArray_ucMode = " << sliceOrdering
				       << " ???, setting slice ordering to \"undefined\"" << endl;
				  m_orderOfSlices = UndefinedOrder;
				  break;
				}
			    }
			  else
			    cout << "scan for \"sSliceArray_ucMode\" failed. " << endl;
			}

		      found = csa.value.front().find("sSliceArray.lSize");
		      if( found != string::npos )
			{
			  const char* p = csa.value.front().c_str() + found;
			  int n, lSize;
			  if( ( n = sscanf( p, "sSliceArray.lSize = %d", &lSize ) ) )
			    {
			      /* this is what we are looking for:
				 sSliceArray.asSlice[0].dThickness	 = 	4.0
				 sSliceArray.asSlice[0].dPhaseFOV	 = 	192.0
				 sSliceArray.asSlice[0].dReadoutFOV	 = 	192.0
				 sSliceArray.asSlice[0].sPosition.dSag	 = 	2.582277852
				 sSliceArray.asSlice[0].sPosition.dCor	 = 	-34.69273759
				 sSliceArray.asSlice[0].sPosition.dTra	 = 	-70.04539956
				 sSliceArray.asSlice[0].sNormal.dSag	 = 	0.03339980043
				 sSliceArray.asSlice[0].sNormal.dCor	 = 	0.3713450025
				 sSliceArray.asSlice[0].sNormal.dTra	 = 	0.9278940362
				 sSliceArray.asSlice[1].dThickness	 = 	4.0
			      */
#ifdef DEBUG
			      cout << "XProtocol ImgPositions:" << endl;
#endif
			      m_ImgPositionSlice.resize(lSize);
			      for( int slice = 0; slice < lSize; slice++ )
				{
				  double dSag,dCor,dTra;

				  found = csa.value.front().find("sPosition.dSag", found );
				  if( found != string::npos )
				    {
				      p = csa.value.front().c_str() + found;
				      if( ( n = sscanf( p, "sPosition.dSag = %lf", &dSag ) ) )
					{
					  found = csa.value.front().find("sPosition.dCor", found );
					  p = csa.value.front().c_str() + found;
					  if( ( n = sscanf( p, "sPosition.dCor = %lf", &dCor ) ) )
					    {
					      found = csa.value.front().find("sPosition.dTra", found );
					      p = csa.value.front().c_str() + found;
					      if( ( n = sscanf( p, "sPosition.dTra = %lf", &dTra ) ) )
						{
						  m_ImgPositionSlice[slice].x = dSag;
						  m_ImgPositionSlice[slice].y = dCor;
						  m_ImgPositionSlice[slice].z = dTra;
#ifdef DEBUG
						  cout << slice << ": " << dSag << ", " << dCor << ", " << dTra << endl;
#endif
						}
					    }
					}
				    }
				}

			    }
			}

		    }
		  else if( csa.name == "Isocentered" )	          {}
		  else if( csa.name == "CoilForGradient2" )        {}
		  else if( csa.name == "ProtocolChangeHistory" )   {}
		  else if( csa.name == "GradientMode" )            {}
		  else if( csa.name == "FlowCompensation" )        {}
		  else if( csa.name == "PostProcProtocol" )        {}
		  else if( csa.name == "RFSWDOperationMode" )	  {}
		  else if( csa.name == "RFSWDMostCriticalAspect" ) {}
		  else if( csa.name == "SARMostCriticalAspect" )   {}
		  else if( csa.name == "TablePositionOrigin" )     {}
		  else if( csa.name == "VFModelInfo" )  	          {}
		  else if( csa.name == "VFSettings" )		  {}
		  else if( csa.name == "AutoAlignData" )		  {}
		  else if( csa.name == "FmriModelParameters" )	  {}
		  else if( csa.name == "FmriModelInfo" )		  {}
		  else if( csa.name == "FmriExternalParameters" )  {}
		  else if( csa.name == "FmriExternalInfo" )	  {}
		  else if( csa.name == "B1rms" )		          {}
		  else if( csa.name == "B1rmsSupervision" )	  {}
		  else if( csa.name == "TalesReferencePower" )	  {}
		  else if( csa.name == "PhaseSliceOversampling" )  {}
		  else if( csa.name == "SafetyStandard" )	  {}
		  else if( csa.name == "DICOMImageFlavor" )	  {}
		  else if( csa.name == "DICOMAcquisitionContrast" ){}
		  else if( csa.name == "EchoTrainLength" )	  {}
		  else if( csa.name == "RFEchoTrainLength" )	  {}
		  else if( csa.name == "GradientEchoTrainLength" ) {}
		  else if( csa.name == "Laterality4MF" )		       {}
		  else if( csa.name == "ArterialSpinLabelingContrast" ) {}
#ifdef CSADEBUG
		  else
		    cout << "\"" << csa.name << "\" unimplemented."  << endl;
#endif
		}
	      }
	  }
      }
      break;
    }
}

void i2o_dcm::Grp51(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x100b: // mosaic tile matrix size, eg. "64*64"
      // aka AcquisitionMatrixText
      // our primary source for mosaic tile size
      
      if( m_Slices > 1 ) // only relevant for mosaic and we should definitely know about it already now
	{
	  int   cols = atoi(b);
	  int   rows = -1;
	  char  p[ len + 1 ];
	  
	  // sscanf( b, "%d*%d", &cols, &rows); // does not work with s.th.like "84p*84"
	  // so we better scan it by ourselve
	  
	  strncpy( p, b, len ); 	// work on a private copy
	  
	  int i;
	  for( i = 0; i < len; i++ )
	    if( ! isdigit( p[i] ) )
	      break;
	  for( ; i < len; i++ )
	    if( isdigit( p[i] ))
	      {
		rows = atoi(p+i);
		break;
	      }
	  
	  if( m_Cols != cols )
	    {
	      if( m_Cols != -1 )
#ifdef DEBUG
		cout << m_filename << "(0x51:0x100b) Resetting m_Cols from " << m_Cols
		     << " to " << cols << endl;
#endif
	      m_Cols = cols;
	    }
	  
	  if( rows != -1 && m_Rows != rows )
	    {
#ifdef DEBUG
	      if( m_Rows != -1 )
		cout << m_filename << "(0x51:0x100b): Resetting m_Rows from " << m_Rows
		     << " to " << rows << endl;
#endif
	      m_Rows = rows;
	    }
	}
      break;

    case 0x100c:
      sscanf( b, "FoV %d*%d", &m_FoVx, &m_FoVy);
      break;

    case 0x100e: // Image Orientation
      {
	size_t m = min( 3, len );
	
	if( strncasecmp(b,"sag", m ) == 0 )
	  m_ImgOrientation="sagittal";
	else if( strncasecmp(b,"cor", m ) == 0 )
	  m_ImgOrientation="coronal";
	else if( strncasecmp(b,"tra", m ) == 0 )
	  m_ImgOrientation="transversal";
	else
	  {
	    //m_ImgOrientation="oblique"; // or "unknown" ? Need a test file...
	    m_ImgOrientation.assign(b,len);
	    trim(m_ImgOrientation);
	  }
      }
      break;
      
    default:
      break;
    }
  /*
    (0051,0010) LO [SIEMENS MR HEADER ]                               # 18,1 Private Creator
    (0051,1008) CS [IMAGE NUM 4 ]                                     # 12,1 CSA Image Header Type
    (0051,1009) LO [1.0 ]                                             # 4,1 CSA Image Header Version ??
    (0051,100a) LO (SH) [TA 01.28]                                    # 8,1 ? (Time of Acquisition ?)
    (0051,100b) LO (SH) [128*128 ]                                    # 8,1 AcquisitionMatrixText
    (0051,100c) LO (SH) [FoV 2112*2112 ]                              # 14,1 ? (Field of View ?)
    (0051,100d) SH [SP F99.7]                                         # 8,1 Slice Location (?)
    (0051,100e) LO (SH) [Tra ]                                        # 4,1 ? (Image Orientation ?)
    (0051,100f) LO [C:R01-30]                                         # 8,1 CoilString
    (0051,1011) LO [p8]                                               # 2,1 PATModeText
    (0051,1012) SH [TP 0]                                             # 4,1 ? (ImaPATModeText ?)
    (0051,1013) SH [+LPH]                                             # 4,1 PositivePCSDirections
    (0051,1016) LO [p8  ]                                             # 4,1 ? (ImageTypeText ?)
    (0051,1017) SH [SL 1.5]                                           # 6,1 ? (SliceThicknessText ?)
    (0051,1019) LO [A1]                                               # 2,1 ? (ScanOptionsText ?)
  */
}

void i2o_dcm::Grp7FE0(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x10:
      m_payload     = b;
      m_payloadSize = len;
      break;
    }
}

void i2o_dcm::Grp7FE1(char* b, u_int16_t element, int32_t len, char VR[3])
{
  switch(element)
    {
    case 0x1010:
      m_payload     = b;
      m_payloadSize = len;
      break;
    }
}

//
// initial reading, just enough to either identify or reject it
//
bool
i2o_dcm::Read( string fn, bool deferred /*==true*/)
{
  const string failedTxt = "\"" + fn + "\" most probably not a valid dicom file, parsing failed.";

  Close(); // forget about the past

  if( ! Open( fn ) )
    return false;

  char* b = (char*) m_mm;
  char* B = b;

  //bool dicomSignatureFound = false;

  // according to the standard the first 128 bytes should be zero
  // followed by a four letter word:

  // "DICM" should be @128, odd enough but ...
  // as long as all these variants still exist somewhere we try as hard as we can ...
  if(strncmp((char*)( b ), "DICM", 4 ) == 0 )
    {
      cout << fn << ": signature \"DICM\" at non standard offset 0, strange." << endl;
      b += 4;
    }
  else if(strncmp((char*)( b + 128 ), "DICM", 4 ) == 0 )
    {
#ifdef DEBUG
      cout << fn << ": signature \"DICM\" found as expected." << endl;
#endif
      b += 128 + 4;
    }
  else if( strncmp((char*)( b + 256 ), "DICM", 4 ) == 0 )
    {
      b += 256 + 4;
      cout << fn << ": signature \"DICM\" at non standard offset 256, strange." << endl;
    }
  else if( strncmp((char*)( b + 512 ), "DICM", 4 ) == 0 )
    {
      b += 512 + 4;
      cout << fn << ": signature  \"DICM\" at non standard offset 512, strange." << endl;
    }
  else
    for( char* dicm = b; dicm < b + 1024; dicm++)
      if( strncmp((char*)( dicm ), "DICM", 4 ) == 0 )
	{
	  cout << fn << ": signature  \"DICM\" at non standard offset " << dicm - b << ", strange." << endl;
	  b = dicm + 4;
	}

  // if( b != B )
  //    dicomSignatureFound = true;
  // accept a missing "DICM" so old NEMA files might still be readable

  // special dicom perversion found in a philips file
  // if( *(u_int32_t*)b == 0xf3f3fef3 ) // or was it "0x3f3f3f3f" ?
  //   {
  //     cout << fn << " !!! " << __FILE__ << __LINE__
  // 	   << ": dicom format oddity, please report." << endl;
  //     b += 4;
  //   }

  m_firstDcmTag = b; // remember the initial position as the starting point for later use

  long examArea = sysconf(_SC_PAGE_SIZE); // one memory page should be enough to get ..
  // .. all information needed for listing and sorting dicom files

  char* end = (char*)m_mm + ( m_filesize < examArea ? m_filesize : examArea );

  const int MinimumTokens = 24; // if we reach that nice number it really must be a dicom

  int       tokenCnt = 0;
  int32_t   len      = 0;

  do
    {
      u_int16_t  grp;
      u_int16_t  element;
      char       VR[3]  = { '\0', '\0', '\0' };
      string     s_b;

      B = b + len; // the next dicom tag candidate

      b = DcmElement( B, grp, element, len, VR ); // fetch it

      if( b == nullptr )
	{
	  Close();
	  m_msgTxt = failedTxt;

	  return false;
	}
      if( len < 0 ){
	// "unknown sequence length" or delimiters or whatever ... don't care
	len = 0;
      }
      else
	{
	  if( b + len >= end )
	    break;

#ifdef DEBUG
	  cout << "i2o_dcm::Read(): (0x" << hex << grp << ",0x" << element << ") len=" << dec << len << endl;
#endif
	  switch(grp)
	    {
	    case /*GRP*/ 0x2:
	      switch(element)
		{
		case 0x02: // MediaStorageSOPClassUID
		  {
		    if( strcmp("1.2.840.10008.1.3.10",b) == 0 )
		      { // Media Storage Directory Storage
			m_isDicomdir = 1;
#ifdef DEBUG
			cout << "\"" << m_filename << "\" looks like a DICOMDIR file to me." << endl;
#endif
		      }
		    else
		      {
		      }
		  }
		  break;
		case 0x3: // MediaStorageSOPInstanceUID
		  {
		  }
		  break;
		case 0x10: // TransferSyntaxUID
		  {
		    // cout << "Transfer Syntax: " << b << endl;
		    // 1.2.840.10008.1.2 : Raw data, Implicit VR, Little Endian
		    // 1.2.840.10008.1.2.x : Raw data, Eplicit VR  x = 1: Little Endian  x = 2: Big Endian
		    // 1.2.840.10008.1.2.4.xx : JPEG compression   xx = 50-64: Lossy JPEG   xx = 65-70: Lossless JPEG
		    // 1.2.840.10008.1.2.5 : Lossless Run Length Encoding
		    if( strcmp("1.2.840.10008.1.2.1",b) == 0 )
		      { // cout <<  ": Raw data, Eplicit VR, Little Endian." << endl;
			m_littleEndian = true;
		      }
		    else if( strcmp("1.2.840.10008.1.2.2",b) == 0 )
		      {
			m_littleEndian = false;
			cout <<  ": Raw data, Eplicit VR, BIG ENDIAN !!! NOT IMPLEMENTED, GARBAGE FOLLOWS ..." << endl;
		      }
		    else if( strcmp("1.2.840.10008.1.2",b) == 0 )
		      {
			// implicit VR, (hopefully still) handled correctly implicitely in DcmElement().
			m_littleEndian = true;
		      }
		    else if( strcmp("1.2.840.10008.1.2.4.5",b) == 0 )
		      {
			cout <<  ": JPEG lossy compression. !!! NOT IMPLEMENTED, GARBAGE FOLLOWS ..." << endl;
		      }
		    else if( strcmp("1.2.840.10008.1.2.4.6",b) == 0 )
		      {
			cout <<  ": JPEG lossless compression. !!! NOT IMPLEMENTED, GARBAGE FOLLOWS ..." << endl;
		      }
		    else
		      {
			cout << m_filename <<  ":\"" << b << "\". Ignoring unknown Transfer Syntax ..." << endl;
		      }
		  }
		  break;
		case 0x12: // ImplementationClassUID
		  {
		  }
		  break;
		case 0x13: // ImplementationVersionName
		  {
		  }
		  break;
		}
	      break;
	    case /*GRP*/ 0x4:
	      switch(element)
		{
		case 0x1200: // OffsetOfTheFirstDirectoryRecordOfTheRootDirectoryEntity
		case 0x1202: // OffsetOfTheLastDirectoryRecordOfTheRootDirectoryEntity
		case 0x1212: // FileSetConsistencyFlag
		case 0x1220: // DirectoryRecordSequence
		  m_isDicomdir = 1;
#ifdef DEBUG
		  cout << "\"" << m_filename << "\" looks like a DICOMDIR file to me." << endl;
#endif
		}
	      break;
	    case /*GRP*/ 0x8:
	      Grp08( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x10:
	      Grp10( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x18:
	      Grp18( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x19: // Siemens Private
	      Grp19( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x20:
	      Grp20( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x28:
	      Grp28( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x29:
	      Grp29( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x51: // Siemens Private
	      Grp51( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x7fe0: // IMAGE GRAPHICS
	      Grp7FE0( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x7fe1: // CSA NON-IMAGE
	      Grp7FE1( b, element, len, VR );
	      break;
	    }
	  tokenCnt++;
	}
    } while ( true /* b + len < end */);

  m_readStoppedHere = B;

#ifdef DEBUG
  cout << m_readStoppedHere - (char*)m_mm << ", after " << tokenCnt << " dcm tokens." << endl;
#endif

  // one of these is mandatory for a valid dicom file
  if( ! ( tokenCnt > MinimumTokens
	  || ( m_isDicomdir == 1 )
	  || ( m_StudyID != -1
	       && m_InstanceNr != -1
	       && m_SeriesNr != -1
	       && m_PatientName.size() > 0
	       )
	  )
      )
    {
      Close();
      m_msgTxt = failedTxt;
      return false;
    }

  // if we're still here it must be either yes or no, ahmm 1 or 0
  if( m_isDicomdir == -1 )
    m_isDicomdir = 0;

  FigureOutDataType();
  FigureOutImgOrientation();

  if( m_Slices != -1 && m_Slices < 2 ) // no mosaic, adjust cols/rows to avoid confusion
    {
      m_Slices = 1;
      if( m_Rows < 0 )
	m_Rows = m_TotalRows;
      if( m_Cols < 0 )
	m_Cols = m_TotalCols;
#ifdef DEBUG
      cout << "RESETTING m_Cols to " << m_Cols << " and m_Rows to " << m_Rows << endl;
#endif
    }

  // well, maybe we better read everything in right now
  if( ! deferred )
    FullScan();

  return true;
}


// just a SaveAs(othername)
bool
i2o_dcm::Write( string fn )
{
  int fd;

  fd = open( fn.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0666 );

  if( fd == -1 )
    {
      cerr << "open(\"" << fn << "\") failed: " << strerror(errno) << endl;
      return false;
    }

  ssize_t written = write( fd, m_mm, m_filesize );

  close( fd );

  if( written != m_filesize )
    {
      cerr << "write(fd \"" << fn << "\") failed: " << strerror(errno) << endl;
      return false;
    }
  return true;
}


//
// scan the whole file and collect everything we need
//
void
i2o_dcm::FullScan()
{
#ifdef DEBUG
  cout << "FullScan()ning " << m_filename << endl;
#endif
  if( m_firstDcmTag ) // that means the initial Read() already happened
    {
      m_fullScanDone = true; // slightly optimistic ;)

      int32_t   len   = 0;
      u_int16_t grp, element;
      char      VR[3] = { '\0', '\0', '\0' };
      char*     end   = (char*)m_mm + m_filesize;

      char* b = m_readStoppedHere;

#ifdef DEBUG
      string s_b;
#endif

      do {
	
	b = DcmElement( b + len, grp, element, len, VR );
	
#ifdef DEBUG
	cout << "i2o_dcm::FullScan():"
	  " grp=0x" << hex << grp <<
	  ", ele=0x" << element <<
	  ", len=" << dec << len << endl;
#endif
	
	if( len < 0 ) // "unknown sequence length" or delimiters or whatever ...
	  len = 0;
	else
	  switch(grp)
	    {
	    case /*GRP*/ 0x8:
	      Grp08( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x10:
	      Grp10( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x18:
	      Grp18( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x19: // Siemens Private
	      Grp19( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x20:
	      Grp20( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x28:
	      Grp28( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x29:
	      Grp29( b, element, len, VR );
	        break;
	    case /*GRP*/ 0x51: // Siemens Private
	      Grp51( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x7fe0: // IMAGE GRAPHICS
	      Grp7FE0( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x7fe1: // CSA NON-IMAGE
	      Grp7FE1( b, element, len, VR );
	      break;
	    case /*GRP*/ 0x9999:
	      switch(element)
		{
		case 0x01:
		  break;
		}
	      break;
	    }
      }
      while ( b + len < end );

      FigureOutDataType();

      if( m_Slices < 2 ) // No mosaic, adjust cols/rows to avoid confusion.
	{   // We have to trust the values for total data matrix size.
	  m_Slices = 1;
	  m_Rows = m_TotalRows;
	  m_Cols = m_TotalCols;
#ifdef DEBUG
	  cout << "RESETTING m_Cols to " << m_Cols << " and m_Rows to " << m_Rows << ", " << endl;
#endif
	}

#ifdef DEBUG
      cout << "FullScan(): \"" << m_filename << "\", size=" << m_filesize << ", patient name=\""
	   << m_PatientName << "\", " << endl
	   << "        cols=" << m_Cols << ", rows=" << m_Rows << ", slices=" << m_Slices
	   << ", totalCols=" << m_TotalCols << ", totalRows=" << m_TotalRows << endl;
#endif
    }
  else
    cerr << "Internal error: FullScan() without prior Read()." << endl;
} // FullScan()


void
i2o_dcm::FigureOutImgOrientation()
{
  if( m_ImgOrientation.size() == 0)
    {
      if( m_ImgOrientRow.size() && m_ImgOrientCol.size() && m_ImgOrientNormal.size() )
	{
	  //TODO
	  //idea:
	  //      axial    1/0/0/0/1/0
	  //      sagittal 0/1/0/0/0/-1
	  //      coronal  1/0/0/0/0/-1
	  // check scalar product imgOrientNormal()*ideal[Axial|Sagital|Coronal]Normal() < threshold
	  // "transversal","coronal","sagital","oblique"
	}
    }
}

void
i2o_dcm::FigureOutDataType()
{
  // might have happened already during the parsing
  if( m_payload && ( m_dataTypeHashCode == 0 ) )
    {
      switch( m_bitsAllocated )
	{
	case 8:
	  if( m_signedDataValues )
	    m_dataTypeHashCode = typeid(signed char).hash_code();
	  else
	    m_dataTypeHashCode = typeid(unsigned char).hash_code();
	  return;
	case 16:
	  if( m_signedDataValues )
	    m_dataTypeHashCode = typeid(signed short).hash_code();
	  else
	    m_dataTypeHashCode = typeid(unsigned short).hash_code();
	  return;
	case 32:
	  if( m_signedDataValues )
	    m_dataTypeHashCode = typeid(signed int).hash_code();
	  else
	    m_dataTypeHashCode = typeid(unsigned int).hash_code();
	  return;
	default:
	  if( m_isDicomdir != 1 )
	    cerr << __FILE__ << ": cannot figure out data type for file \"" << m_filename << "\""<< endl
		 << "bitsAllocated=" << m_bitsAllocated << ", bitsStored=" << m_bitsStored
		 << ", samplesPerPixel=" << m_SamplesPerPixel << ", rows=" << m_Rows << ", cols=" << m_Cols
		 << endl;
	}
    }
}


void
i2o_dcm::Dump()
{
  if( m_firstDcmTag )
    {
      char* b = m_firstDcmTag;

      int32_t   len = 0;
      u_int16_t grp, element;
      char      VR[3] = { '\0', '\0', '\0' };
      char*     end = (char*)m_mm + m_filesize;

      do
	{
	  b = DcmElement( b + len, grp, element, len, VR );

	  cout << dec << "@" << b - (char*)m_mm
	       << hex << ": (0x" << grp << ",0x" << element << ") " << VR
	       << dec << ", len=" << len << endl;

	  if( len < 0 )
	    len = 0;
	}
      while ( b + len <= end );
    }

}


char*
i2o_dcm::Grep( int g, int e, int& size )
{
  if( m_firstDcmTag )
    {
      char* b = m_firstDcmTag;

      int32_t    len = 0;
      u_int16_t  grp, element;
      char       VR[3]  = { '\0', '\0', '\0' };
      char*      end = (char*)m_mm + m_filesize;

      do
	{
	  b = DcmElement( b + len, grp, element, len, VR );

	  if( (grp == g) && (element == e) )
	    {
	      size = len;
	      return b;
	    }

	  if( len < 0 )
	    len = 0;
	}
      while ( b + len < end );
    }
  return nullptr;
}


void
i2o_dcm::Close()
{
  if( m_mm )
    {
      if( m_filesize && munmap( m_mm, m_filesize ) == -1 )
	cerr << "munmap(\"" << m_mm << "\") failed: " << strerror(errno) << endl;

      m_PatientName.clear();

    }
}


bool
i2o_dcm::Open( string fn )
{
  int         fd;
  struct stat sb;

  fd = ::open( fn.c_str(), O_RDONLY);
  if( fd == -1 )
    {
      m_msgTxt = "open(\"" + fn + "\") failed, because: \"" + strerror(errno) +"\".";
      return false;
    }
  if( fstat(fd, &sb) == -1 )
    {
      ::close( fd );
      m_msgTxt = "fstat() for \"" + fn + "\" failed, because: \"" + strerror(errno) +"\".";
      return false;
    }
  m_filesize = sb.st_size;

  m_mm = mmap( nullptr, m_filesize, PROT_WRITE, MAP_PRIVATE, fd, 0 ); // leave the file untouched

  ::close( fd ); // we dont need fd anymore, the system keeps its own copy

  if ( m_mm == (void*)(-1) )
    {
      m_mm = nullptr;
      m_msgTxt = "mmap() for \"" + fn + "\" failed, because: \"" + strerror(errno) +"\".";
      return false;
    }

  m_filename = fn;

  return true;
}

double
i2o_dcm::DistanceFromOrigin()
{
  if( ! m_fullScanDone )
    FullScan();
  if( m_ImgPosition.size() )
    return sqrt( m_ImgPosition[0] * m_ImgPosition[0] +
		 m_ImgPosition[1] * m_ImgPosition[1] +
		 m_ImgPosition[2] * m_ImgPosition[2] );
  else
    return -1.0;
}

// a new 'patientName' will be added if needed, or returned if already listed
DicomPatient*
AddOrGetPatient( list<DicomPatient>& patList, string patientName )
{
  for( auto& patient : patList )
    if( patient.patientName == patientName )
      {
	return &patient;
      }
  patList.push_back( patientName );

  return &patList.back();
}


DicomStudy*
AddOrGetStudy( DicomPatient* patient, int studyID, string studyDescription )
{
  for( auto& study : patient->studyList )
    if( study.studyID == studyID )
	return &study;

  patient->studyList.push_back( DicomStudy( studyID, studyDescription ) );

  return &(patient->studyList.back());
}


DicomSeries*
AddOrGetSeries( DicomStudy* study, int seriesNr, string seriesDescription )
{
  // check if it is already there
  for( auto& series : study->seriesList )
    if( series.seriesNr == seriesNr )
	return &series;

  // sort it in
  //  for( const auto& series : study->seriesList)
  for( auto series = study->seriesList.begin(), list_end = study->seriesList.end(); series != list_end; ++series)
    if( series->seriesNr > seriesNr )
      {
	return &(*study->seriesList.emplace( series, seriesNr, seriesDescription ));
      }
  // or append
  study->seriesList.push_back( DicomSeries( seriesNr, seriesDescription ) );
  return &study->seriesList.back();
}


DicomFile*
AddImageFile( DicomSeries* series, int imgNr, string path )
{
  // check if it is already there
  //for( auto img = series->fileList.begin(), list_end = series->fileList.end(); img != list_end; ++img)
  for( auto &img : series->fileList )
    if( img.imgNr == imgNr && img.path == path)
      {
	// identical image/instance numbers may happen (e.g. gre_field_mapping)
	// and because we're still not clear about what's the best thing to do then, we complain:
	cerr << "request for adding already listed image \"" << imgNr << "\" \"" << img.path
	     << "\" into series \"" << series->seriesNr << "\"" << endl;
	return &img;
      }
  // sort it in
  for( auto img = series->fileList.begin(), list_end = series->fileList.end(); img != list_end; ++img )
    if( img->imgNr > imgNr)
      {
	// ( it's always an "insert before" )
	return &(*series->fileList.insert( img, DicomFile( imgNr, path ) ) );
      }
  // or append
  series->fileList.push_back( DicomFile( imgNr, path ) );

  return &series->fileList.back();
}


inline 

bool
i2o_dcm::IsDicomDir()
{
  if( m_isDicomdir == -1 )
    if( !m_fullScanDone )
      FullScan();
  return ( m_isDicomdir == 1 );
}


// scan a DICOMDIR file
int		                          // return total number of image/data files referenced.
i2o_dcm::DicomDir( DicomDirHead& dcmDir ) // fill scan result into a dicomDir object.
{
  // remark: the following is based on a low documentation level and a long hex dump meditation
  //         and so it might need some cleanup and lot more testing ... but it worked so far

  // if not already done:
  if( m_firstDcmTag == nullptr )
    if( ! Read( m_filename ) )
      return 0;

  if( ! IsDicomDir() )
    return 0;

  int           imageCnt = 0;

  DicomPatient* patient = nullptr;
  DicomStudy*   study   = nullptr;
  DicomSeries*  series  = nullptr;
  DicomFile*    file    = nullptr;

  string        seriesDescription;
  string        studyDescription;
  string        imgPath;
  
  int           imgNr         = -1;
  int           acquisitionNr = -1;

  int           previewBits   = 0;

  string        dirName;
  char          imgFilename[ PATH_MAX + 1 ];

  char*         cwd = get_current_dir_name();

  // extract the path from dicomdir's name
  strncpy( imgFilename, m_filename.c_str(), PATH_MAX + 1 );
  
  dirName.assign( dirname( imgFilename) );

  if( chdir( dirName.c_str() ) == -1 )
    { // should never happen ..., at least give a warning
      m_msgTxt = "chdir(\"" + dirName + "\") failed: " + strerror(errno);
      cerr << m_msgTxt << endl;
    }

  int32_t   len = 0;
  u_int16_t grp, element;
  char      VR[3]  = { '\0', '\0', '\0' };
  char*     end    = (char*)m_mm + m_filesize;

  const size_t dirRecTypeLen = 32;
  char dirRecType[ dirRecTypeLen + 1 ];

  int dirRecSeq = 0;
  int refImgSeq = 0;
  int imgIconSeq = 0;
  int item = 0;
  int dirRecItem = 0;
  int refImgItem = 0;
  int imgIconItem = 0;
  int seriesCnt = 0;
  int studyCnt = 0;
  int refFileIdCnt = 0;

  char* b = m_firstDcmTag;

  do {

    b = DcmElement( b + len, grp, element, len, VR );

    if( len < 0 )
      len = 0;

    if( ( grp == 0x04) && ( element == 0x1220 ) ) // DirectoryRecordSequence
      {
#ifdef DEBUG
	cout << "Directory Sequence" << endl;
#endif
	dirRecSeq++;
      }
    else if( ( grp == 0xfffe) && ( element == 0xe000 ) ) // Item
      {
	item++;
	if( imgIconSeq )
	  {
#ifdef DEBUG
	    cout << "ImageIcon" << endl;
#endif
	    imgIconItem++;
	  }
	else if ( refImgSeq )
	  {
#ifdef DEBUG
	    cout << "ReferencedImage" << endl;
#endif
	    refImgItem++;
	  }
	else if( dirRecSeq )
	  {
#ifdef DEBUG
	    cout << "DirectoryItem" << endl;
#endif
	    dirRecItem++;
	  }
	else
	  {
	    cout << "ITEM WITHOUT CONTEXT" << endl;
	  }
      }
    else if( ( grp == 0xfffe) && ( element == 0xe00d ) ) // ItemDelimitationItem
      {
	//--item;
	if( imgIconSeq )
	  {
#ifdef DEBUG
	    cout << "ImageIcon end" << endl;
#endif
	    //--imaIconItem;
	  }
	else if ( refImgSeq )
	  {
#ifdef DEBUG
	    cout << "ReferencedImage end" << endl;
#endif
	    //	    --refImaItem;
	  }
	else if( dirRecSeq )
	  {
#ifdef DEBUG
	    cout << "DirectoryItem end" << endl;
#endif
	    //--dirRecItem;
	  }
	else
	  {
	    cerr << "ITEM DELIMITER WITHOUT CONTEXT" << endl;
	  }
      }
    else if( ( grp == 0x04) && ( element == 0x1430 ) ) // DirectoryRecordType
      {
	strncpy( dirRecType, b, min( dirRecTypeLen, (size_t)len) );
	dirRecType[ min( dirRecTypeLen, (size_t)len) ] = '\0';
      }
    else if( ( grp == 0x08) && ( element == 0x1140 ) ) // ReferencedImageSequence
      {
#ifdef DEBUG
	cout << "ReferencedImageSequence" << endl;
#endif
	refImgSeq++;
      }
    else if( ( grp == 0xfffe) && ( element == 0xe0dd ) ) // SequenceDelimitationItem
      {
	if( imgIconSeq )
	  {
#ifdef DEBUG
	    cout << "ImageIconSequence end" << endl;
#endif
	    --imgIconSeq;
	  }
	else if ( refImgSeq )
	  {
#ifdef DEBUG
	    cout << "ReferencedImageSequence end" << endl;
#endif
	    --refImgSeq;
	  }
	else if( dirRecSeq )
	  {
#ifdef DEBUG
	    cout << "Directory Sequence end" << endl;
#endif
	    --dirRecSeq;
	  }
	else
	  {
	    cerr << "SEQUENCE DELIMITER WITHOUT CONTEXT" << endl;
	  }

      }
    else if( ( grp == 0x88) && ( element == 0x0200 ) ) // ImageIconSequence
      {
#ifdef DEBUG
	cout << "ImageIconSequence" << endl;
#endif
	imgIconSeq++;
      }
    else if( len && ( grp == 0x10 ) && ( element == 0x10 ) ) // patient name
      {
	char c[len+1];
	memcpy( c, b, len );
	c[len] = '\0';

	patient = AddOrGetPatient( dcmDir.patientList, trim(c) );
#ifdef DEBUG
	cout << "Patient \"" << patient->patientName << "\" added to list" << endl;
#endif
	study  = nullptr;
	series = nullptr;
	file   = nullptr;
      }
    else if( len && ( grp == 0x08 ) && ( element == 0x1030 ) ) // study description
      {
	studyDescription.assign(b,len);
	trim(studyDescription);
      }
    else if( len && ( grp == 0x08 ) && ( element == 0x103e ) ) // series description
      {
	seriesDescription.assign(b,len);
	trim(seriesDescription);
      }
    else if( len && ( grp == 0x20 ) && ( element == 0x10 ) ) // study ID
	{
	  char c[ len+1 ];
	  
	  memcpy( c, b, len );
	  c[len] = '\0';

	  studyCnt++;

	  if( patient )
	    {
	      study  = AddOrGetStudy( patient, atoi(c), studyDescription );
	      series = nullptr;
	      file   = nullptr;
#ifdef DEBUG
	      cout << "Study \"" << study->studyID << "\" added to list" << endl;
#endif
	    }
	  else
	    cerr << "ERROR: study (\"" << c << "\") without Patient. Ignoring ..." << endl;
	}
    else if( len && ( grp == 0x20 ) && ( element == 0x11 ) ) // series number
      {
	char c[ len+1 ];
	
	memcpy( c, b, len );
	c[len] = '\0';

	if( study )
	  {
	    seriesCnt++;
	    series = AddOrGetSeries( study, atoi(c), seriesDescription );
	    file = nullptr;
#ifdef DEBUG
	    cout << "Series \"" << series->seriesNr << "\" added to list" << endl;
#endif
	  }

	else
	  cerr << "ERROR: series (\"" << c << "\") without Study. Ignoring ..." << endl;

      }
    else
      if( len && ( grp == 0x20 ) && ( element == 0x12 ) ) // acquisition number
	{
	  char c[ len+1 ];
	  
	  memcpy( c, b, len );
	  c[len] = '\0';

	  acquisitionNr = atoi(c);
	  // HACK ! The dirRecType=PRIVATE items do not have an instance number but an
	  // acquisition number instead, so we cheat here to stay happy without further modifications
	  imgNr = acquisitionNr;

	}
    if( len && ( grp == 0x20 ) && ( element == 0x13 ) ) // image/instance number
      {
	char c[len+1];
	memcpy( c, b, len );
	c[len] = '\0';

	imgNr = atoi(c);

      }
    else
      if( len && ( grp == 0x04 ) && ( element == 0x1500 ) ) // referenced file id
	{
	  refFileIdCnt++;

	  strncpy( imgFilename, b, min( len, PATH_MAX + 1 ));
	  imgFilename[ min(  len, PATH_MAX ) ] = '\0';
	  trim(imgFilename);

	  // some pathname cleanups
	  // unix rules ... path separator correction:
	  for( char* p = imgFilename; *p; p++)
	    if( *p == '\\' )
	      *p = '/';
#if 0
	  // we might try to correct for upper/lower case issues
	  // but that would be much more than that:
	  {
	    struct stat sb;
	    if( stat( imgFilename, &sb ) == -1 )
	      { // one shot, try lower case
		strncpy( tName, imgFilename, PATH_MAX-1 ); // save the original
		for( char* p=imgFilename; *p; p++ )
		  *p = tolower( *p );
		// check again
		if( stat( imgFilename, &sb ) == -1 )
		  { // still failing, give up here, restore name
		    strcpy( imgFilename, tName );
		  }
	      }
	  }
#endif
	  imgPath = dirName + "/" + imgFilename;
	  
	}
      else if ( len && ( grp == 0x28 ) && ( element == 0x10 ) )
	{
	  if( file )
	    file->preview.rows = *(u_int16_t*) b;
	  //cout << "rows = " << previewRows << endl;
	}
      else if ( len && ( grp == 0x28 ) && ( element == 0x11 ) )
	{
	  if( file )
	    file->preview.columns = *(u_int16_t*) b;
	  //cout << "cols = " << previewCols << endl;
	}
      else if ( len && ( grp == 0x28 ) && ( element == 0x100 ) )
	{
	  previewBits = *(u_int16_t*) b;
	  if( previewBits > 8 )
	    {
	      file->preview.columns = 0;
	      file->preview.rows = 0;
	      m_msgTxt = "preview data only implemented for at most 8 bits/pixel";
	      //+ to_string( previewBits ) + " requested.";
	    }

	  //cout << "bits = " << previewBits << endl;
	}
      else if ( len && ( grp == 0x7fe0 ) && ( element == 0x010 ) )
	{
	  if( file && previewBits == 8 )
	    file->preview.data = new unsigned char[ file->preview.columns * file->preview.rows ];
	  //cout << "data = " << previewData << endl;
	}

    if( imgNr >= 0 && imgPath.length() )
      {
	if( series )
	  {
	    file    = AddImageFile( series, imgNr, imgPath );
	    imgNr   = -1;
	    imgPath.clear();
	    imageCnt++;
#ifdef DEBUG
	    cout << "Imagefile \"" << file->imgNr << "\" \"" << file->path << "\" added to list" << endl;
#endif
	  }
	else
	  cerr << "ERROR: dicom file reference without Series. Ignoring ..." << endl;
      }
  } while ( b + len < end );

#ifdef DEBUG
  cout << "grp: " << hex << grp << endl;
  cout << "element: " << element << dec << endl;

  cout << "dirRecSeq: " << dirRecSeq << endl;
  cout << "refImaSeq: " << refImgSeq << endl;
  cout << "imaIconSeq: "<< imgIconSeq << endl;
  cout << "item: " << item << endl;
  cout << "dirRecItem: " << dirRecItem << endl;
  cout << "refImaItem: " << refImgItem << endl;
  cout << "imaIconItem: " << imgIconItem << endl;
  cout << "imageCnt: " << imageCnt << endl;
  cout << "seriesCnt: " << seriesCnt << endl;
  cout << "refFileIdCnt: " << refFileIdCnt << endl;
  cout << "studyCnt: " << studyCnt << endl;
#endif

  if( cwd )
    {
      chdir( cwd );
      free( cwd );
    }

  return imageCnt;
}


int			                  // return total number of dicom data files found
ScanDicomFiles( list<string>& inputFiles, // list of files to put together into a details struct:
		DicomDirHead& dcmDir	  // structure to return the list of patients (see above)
		)
{
  int imageCnt = 0;

  for( auto &filename : inputFiles )
    {
      struct stat sb;
      
      if( stat( filename.c_str(), &sb ) == -1 )
	{
	  cerr << "stat(\"" << filename << "\") failed: " << strerror(errno) << endl;
	  continue;
	}

      if( S_ISREG( sb.st_mode ) || S_ISLNK( sb.st_mode ) )
	{
	  i2o_dcm  dcm;
	  
	  if( dcm.Read( filename ) )
	    if( dcm.IsDicomDir() == 0 ) // skip DICOMDIR files 
	      {				// ( we should perhaps skip other meta files too ... or not ? )
		string patientName = dcm.PatientName();
		int    studyID;
		int    seriesNr;
		int    instanceNr;
		
		if( patientName.size() )
		  if( (studyID = dcm.StudyID() ) > 0 )
		    if( (seriesNr = dcm.SeriesNr() ) > 0 )
		      if( (instanceNr = dcm.InstanceNr() ) > 0 )
			{
			  DicomPatient* patient = AddOrGetPatient( dcmDir.patientList, patientName );
			  DicomStudy*   study   = AddOrGetStudy( patient, studyID, dcm.StudyDescription() );
			  DicomSeries*  series  = AddOrGetSeries( study, seriesNr, dcm.SeriesDescription() );
			  AddImageFile( series, instanceNr, filename );
			  
			  imageCnt++;
			}
	      }
	}
    }
  return imageCnt;
}


//
// scan a directory for dicom data files (no recursion)
//
int				// return number of total images found
ScanDir4Dicoms( DicomDirHead& dcmDir ) // return a list of patients
{
  string       dname    = dcmDir.baseDir.c_str();
  list<string> inputFiles;
  DIR*         dir;

  // scan the directory and collect all files into a list

  if( (dir = opendir( dname.c_str() )) == nullptr )
    {
      cerr << "opendir(\"" << dname << "\") failed: " << strerror(errno) << endl;
      return 0;
    }

  struct dirent* dp;

  while( ( dp = readdir( dir ) ) != nullptr )
    {
      if( ( strcmp(".", dp->d_name) == 0 ) || ( strcmp("..", dp->d_name) == 0 ))
	continue;

      string filename = dname + "/" + dp->d_name;

      inputFiles.push_back( filename );
    }
  closedir( dir );

  return ScanDicomFiles( inputFiles, dcmDir );
}


int				 // return number of dicom data files found
ScanDicomDir( string path,	 // DICOMDIR file or directory to scan
	      DicomDirHead& dcmDir ) // return patient list
{
  int    retVal = -1;
  struct stat sb;

  if( path.back()=='/')
    // cosmetics
    path.resize(path.size()-1);

  if( stat( path.c_str(), &sb ) == -1 )
    {
      cerr << "stat(\"" << path << "\") failed: " << strerror(errno) << endl;
    }
  else
    {
      int imageCnt=0;

      if( S_ISDIR( sb.st_mode ) )
	{
	  dcmDir.baseDir = path; // a directory populated with little dicoms

	  imageCnt = ScanDir4Dicoms( dcmDir );
	}
      else
	{
	  i2o_dcm dcm( path ); // a DICOMDIR file hopefully

	  imageCnt = dcm.DicomDir( dcmDir );
	}
      retVal = imageCnt;
    }
  return retVal;
}
