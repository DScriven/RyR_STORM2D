/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, tiff6qt.h, is a header file for tiff6qt.cpp, part of the RyR_STORM2D program.
*
* RyR_STORM2D links to the proprietary Qt system (currently ver 5.15.2) as well as the the free
* CGAL algorithmic library and the free TIFF library. It also requires OpenGL ver 4.3 or higher.
* On Windows it links to the Visual C++ redistributable
*
* RyR_STORM2D is free software: you can redistribute it and/or modify it under the terms of the
* GNU General Public License as published by the Free Software Foundation, either version 3 of
* the License, or any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with RyR_STORM2D.  If not, see <https://www.gnu.org/licenses/>.
*
**********************************************************************************************/#if !defined (_tiff6qt_h)
#define _tiff6qt_h

#include <sstream>
#include <QString>

using namespace std;

#include "tiffio.h"

enum {TIFF_Mono_8bit,  TIFF_RGB_8bit,  TIFF_Map_8bit, TIFF_RGBA_8bit,
      TIFF_Mono_16bit, TIFF_RGB_16bit, TIFF_Map_16bit, TIFF_RGBA_16bit};


class TiffImage
{
protected:
     TIFF  *tif;

     bool    bRGB;            // True for 3 or 4 component colour TIFF
     bool    bMap;            // True for colourmapped TIFF
     bool    bByteImage;         // True for 8 bit TIFF image
     bool    bReversedImage;     // True if BLACK is max, & white is min
     bool    bSwitchEndian;
     bool    bWrite;
     bool    bMinMax;
     bool    bReadHeaderOnly;

     uint8_t*  ByteData;
     uint8_t*  tempData;

     uint32_t   nPixelsPerPlane;
     uint32_t   nPixelsPerImage;
     uint32_t   nBytesPerPixel;
     uint32_t   nBytesPerImage;
     uint32_t   nBytesPerPlane;

     uint32_t  nXdim;    // Image X dim
     uint32_t  nYdim;    // Image Y dim
     uint32_t  nNoZPlanes;
     int     nX0;
     int     nY0;

     uint16_t* red;    // colourmap pointers
     uint16_t* green;
     uint16_t* blue;

     uint16_t* nFmin;
     uint16_t* nFmax;
     
     uint16_t  nResolutionUnit;
     uint16_t  nBitsPerSample;
     uint16_t  nCompression;
     uint16_t  nSamplesPerPixel;
     uint16_t  nPhotometric;
     uint16_t  nPlanarConfiguration;
     int16_t   nOrientation;
     

     int     nByteMultiplier;

     float   fXResolution;
     float   fYResolution;
     float   fXPixelSize;
     float   fYPixelSize;

     QString  Software;
     QString  Description;
     QString  DateTime;
     QString  ErrInfo;
     
     char    cChannel;

public:
     TiffImage();
     ~TiffImage();
     
     void    CopyParameters(const TiffImage& orgfile);
     bool    CopyData(const TiffImage& orgfile);

     QString  ImageDescription() {return Description;}
     QString  SoftwareName() {return Software;}
     QString  DateandTime() {return DateTime;}
     void    addDescription(const char* szLine);
     void    addDescription(QString info) {Description.append(info);}
     bool    Segment(int xstart, int xend, int ystart, int yend);
     void    CalculateMinMax();

     bool    WhiteisMin() {return bReversedImage;}
     uint16_t  Photometric() {return nPhotometric;}
     uint16_t  BitsPerSample() {return nBitsPerSample;}
     uint16_t  SamplesPerPixel() {return nSamplesPerPixel;}
     uint    PntsPerPlane(){return nPixelsPerPlane;}
     uint    PntsPerImage(){return nPixelsPerImage;}
     void    ColourMap(uint16_t* &r, uint16_t* &g, uint16_t* &b)
               {r = red; g = green; b = blue;}

     uint8_t*   BytePointer() { return ByteData;}
     int16_t*     ShortPointer(){ return reinterpret_cast<short*> (ByteData);}
     uint16_t*  UShortPointer(){ return reinterpret_cast<uint16_t*> (ByteData);}
     
     void     setBytePointer(uint8_t* pb) { bByteImage=true; ByteData=pb;}
     void     setUShortPointer(uint16_t* ps){ bByteImage = false; ByteData=reinterpret_cast<uint8_t*> (ps);}
     void     setShortPointer(int16_t* ps){ bByteImage=false; ByteData=reinterpret_cast<uint8_t*> (ps);}

     float     XPixelSize() {return fXPixelSize;}
     float     YPixelSize() {return fYPixelSize;}
     int       X0()  { return nX0;}
     int       Y0()  { return nY0;}
     uint32_t    Xdim(){ return nXdim;}
     uint32_t    Ydim(){ return nYdim;}
     uint32_t    Zdim(){ return nNoZPlanes;}

     bool      ReadTiff(QString InputName, char cChannel=' ', int zstart = 0, int zend = 0, int zinc = 1);
     bool      WriteTiff(QString OutputName);
     bool      ReadHeader(QString FileName);
				      
     void      GetMinMax(uint16_t& Minval, uint16_t& Maxval, int val=0){Minval = nFmin[val]; Maxval = nFmax[val];}

     bool      setOutputType(int nConversion, char cChannel = ' ');
     void      setXPixelSize(float xps) {fXPixelSize = xps;}
     void      setYPixelSize(float yps) {fYPixelSize = yps;}
     void      setSamplesPerPixel(uint16_t spp){nSamplesPerPixel = spp;}
     void      setBitsPerSample(uint16_t bps);
     void      setDateTime(const char* _szDateTime);
     void      setDateTime(QString DateandTime);
     void      setDescription(QString _Description) {Description = _Description;}
     void      setSoftware(QString _Software) {Software = _Software;}
     void      setDescription(const char* _szDescription); 
     void      setSoftware(const char* _szSoftware);
     void      setOrigin(int nX, int nY){nX0 = nX; nY0 = nY;}
     void      setPhotometric(uint16_t _nPhotometric){nPhotometric = _nPhotometric;}
     void      setColourMap(uint16_t* r, uint16_t* g, uint16_t* b);
     bool      setDim(uint32_t _nXdim, uint32_t _nYdim, uint32_t nZdim);
     void      setRes(float fRes, uint16_t nResolutionUnit);
     void      close();

protected:

     bool    CheckDirectories();

     float   GetTIFFfloat(ttag_t tag);
     short   GetTIFFshort(ttag_t tag);
     int     GetTIFFint(ttag_t tag);
};

#endif
