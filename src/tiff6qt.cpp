/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, tiff6qt.cpp, is part of the RyR_STORM2D program and acts as an interface with the TIFF
* library, libtiff
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
**********************************************************************************************/
#include "tiff6qt.h"
#include <tiff.h>
#include <QFileInfo>
#include <QTextStream>
#include <QDateTime>
#include <new>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <fcntl.h>

using namespace std;

TiffImage::TiffImage()
{
    nX0=0;
    nY0=0;
    fXResolution = 0.0;
    fYResolution = 0.0;
    fXPixelSize = 0.0;
    fYPixelSize = 0.0;
    nResolutionUnit = 1;
    nOrientation = 1;

    nBytesPerPixel = 1;
    nSamplesPerPixel = 1;
    nBitsPerSample = 8;
    nByteMultiplier = 1;
    bByteImage = true;
    nPhotometric = 1;

    bReversedImage= false;
    bRGB = false;
    bMinMax = false;
    bReadHeaderOnly = false;

    // set arrays to zero
    ByteData = nullptr;
    tempData = nullptr;
    nFmin=nullptr;
    nFmax=nullptr;

    red = nullptr;
    green = nullptr;
    blue = nullptr;
}
TiffImage::~TiffImage()
{
    delete [] tempData;
    delete [] ByteData;
    delete [] nFmin;
    delete [] nFmax;
    delete [] red;
    delete [] green;
    delete [] blue;
}
void TiffImage::CopyParameters(const TiffImage& orgfile)
{
    nResolutionUnit=orgfile.nResolutionUnit;
    nBitsPerSample=orgfile.nBitsPerSample;
    nByteMultiplier=orgfile.nByteMultiplier;
    nCompression=orgfile.nCompression;
    nSamplesPerPixel=orgfile.nSamplesPerPixel;
    nBytesPerPixel=orgfile.nBytesPerPixel;
    nPhotometric=orgfile.nPhotometric;
    nPlanarConfiguration=orgfile.nPlanarConfiguration;
    nOrientation=orgfile.nOrientation;
    
    if(nPhotometric == 3)
    {
        uint num_entries = 1<<nBitsPerSample;
        red = new uint16_t[num_entries];
        green = new uint16_t[num_entries];
        blue = new uint16_t[num_entries];
        for (uint i = 0; i < uint(num_entries); i++)
        {
            red[i] = orgfile.red[i];
            green[i] = orgfile.green[i];
            blue[i] = orgfile.blue[i];
        }
    }

    fXResolution = orgfile.fXResolution;
    fYResolution = orgfile.fYResolution;
    fXPixelSize = orgfile.fXPixelSize;
    fYPixelSize = orgfile.fYPixelSize;

    bByteImage = orgfile.bByteImage;
    bReversedImage= orgfile.bReversedImage;
    bRGB = orgfile.bRGB;
    bMinMax = orgfile.bMinMax;

    Description = orgfile.Description;
    DateTime = orgfile.DateTime;
    Software = orgfile.Software;
}
bool TiffImage::setDim(uint32_t Xdim, uint32_t Ydim, uint32_t Zdim)
{
    nXdim = Xdim;
    nYdim = Ydim;
    nNoZPlanes = Zdim;
    nPixelsPerPlane = nXdim*nYdim;
    nPixelsPerImage = nPixelsPerPlane*nNoZPlanes;
    uint nBytesReqd = nPixelsPerImage*nBytesPerPixel;
    if(ByteData != nullptr)
    {
        if(nBytesPerImage!= nBytesReqd)
        {
            delete [] ByteData;
            try
            {
                ByteData = new uint8_t[nBytesReqd];
            }
            catch (std::bad_alloc& ba)
            {
                cerr << "\nError - Unable to allocate memory for TIFF data - Reason: " << ba.what() << "\n";
                return false;
            }
        }
    }
    else
    {
        try {
            ByteData = new uint8_t[nBytesReqd]; }
        catch  (std::bad_alloc& ba)
        {
            cerr << "\nError - Unable to allocate memory for TIFF data - Reason: " << ba.what() << "\n";
            return false;
        }
    }

    nBytesPerImage = nBytesReqd;
    nBytesPerPlane = nPixelsPerPlane*nBytesPerPixel;
    return true;
}
bool TiffImage::CopyData(const TiffImage& orgfile)
{
    if(!setDim(orgfile.nXdim, orgfile.nYdim, orgfile.nNoZPlanes))
        return false;

    nX0=orgfile.nX0;
    nY0=orgfile.nY0;
    if(bMinMax)
    {
        nFmin = new uint16_t [nSamplesPerPixel];
        nFmax = new uint16_t [nSamplesPerPixel];
        for(int j = 0; j < nSamplesPerPixel; j++)
        {
            nFmin[j]=orgfile.nFmin[j];
            nFmax[j]=orgfile.nFmax[j];
        }
    }

    uint8_t* pB = orgfile.ByteData;
    uint8_t* pBD = ByteData;
    for(uint i=0; i < nBytesPerImage; i++)
        *pBD++=*pB++;

    return true;
}
bool TiffImage::ReadHeader(QString FileName)
{
    bReadHeaderOnly = true;
    return ReadTiff(FileName);
}
bool TiffImage::ReadTiff(QString FileName, char cChannelL, int zstart, int zend, int zinc)
{
    cChannel = cChannelL;
    QString fname;

    QFileInfo fi (FileName);
    if(fi.suffix().isEmpty())
    {
        fname=fi.absoluteFilePath() + ".tif";
        fi.setFile(fname);
        if(!fi.exists())
        {
            fi.setFile(FileName);
            fname=fi.absoluteFilePath() + ".tiff";
            fi.setFile(fname);
            if(!fi.exists())
            {
                cerr << "\nError - File " << qPrintable(FileName) << " does not exist with either the .tiff or .tif extension\n\n";
                return false;
            }
        }
    }
    else
    {
        if(!fi.exists())
        {
            cerr << "\nError - File " << qPrintable(FileName) << " does not exist!\n\n";
            return false;
        }
    }
    if(!fi.isReadable())
    {
        cerr << "\nError - File " << qPrintable(fname) << " is not readable - check permissions\n";
        return false;
    }
    tif = TIFFOpen(qPrintable(fname),"r");
    if(!tif)
        return false;

    uint32_t width;
    uint32_t length;

    // read tags in first header

    TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&length);

    nXdim = width;
    nYdim = length;

    TIFFGetFieldDefaulted(tif,TIFFTAG_BITSPERSAMPLE, &nBitsPerSample);
    if(nBitsPerSample != 8 && nBitsPerSample != 16)
    {

        ErrInfo = QString("Cannot handle BitsPerSample of %1").arg(nBitsPerSample);
        TIFFError("\n ReadTiff",qPrintable(ErrInfo));
        close();
        return(false);
    }
    if(nBitsPerSample == 8)
    {
        bByteImage=true;
        nByteMultiplier=1;
    }
    else
    {
        bByteImage=false;
        nByteMultiplier=2;
    }

    TIFFGetFieldDefaulted(tif,TIFFTAG_ORIENTATION, &nOrientation);

    TIFFGetFieldDefaulted(tif,TIFFTAG_SAMPLESPERPIXEL, &nSamplesPerPixel);
    if(nSamplesPerPixel == 1)
        bRGB = false;
    else if(nSamplesPerPixel == 3)
        bRGB = true;
    else
    {
        ErrInfo = QString(" Can't handle Samples Per Pixel = %1").arg(nSamplesPerPixel);
        TIFFError("\n ReadTiff",qPrintable(ErrInfo));
        close();
        return(false);
    }

    nBytesPerPixel = nSamplesPerPixel*nByteMultiplier;
    uint16_t nSampleFormat;
    TIFFGetFieldDefaulted(tif,TIFFTAG_SAMPLEFORMAT, &nSampleFormat);
    if(nSampleFormat > 1 && nSampleFormat < 4)
    {
        ErrInfo = QString("Can't handle Sample Format = %1").arg(nSampleFormat);
        TIFFError("\n ReadTiff",qPrintable(ErrInfo));
        close();
        return(false);
    }

    TIFFGetField(tif,TIFFTAG_PHOTOMETRIC, &nPhotometric);
    if(nPhotometric != 2 && bRGB)
    {
        ErrInfo = QString("Photometic interpretation of %1 not supported").arg(nPhotometric);
        TIFFError("\n ReadTiff",qPrintable(ErrInfo));
        close();
        return(false);
    }

    if(nPhotometric == 3)
    {
        uint16_t* r;
        uint16_t* g;
        uint16_t* b;

        if(!TIFFGetField(tif, TIFFTAG_COLORMAP, &r, &g, &b))
        {
            TIFFError("\n ReadTiff","Missing color map!");
            close();
            return(false);
        }

        delete [] red;
        delete [] green;
        delete [] blue;

        uint num_entries = 1<<nBitsPerSample;
        red = new uint16_t[num_entries];
        green = new uint16_t[num_entries];
        blue = new uint16_t[num_entries];
        for (size_t i = 0; i < num_entries; i++)
        {
            red[i] = r[i];
            green[i] = g[i];
            blue[i] = b[i];
        }
    }

    if(nPhotometric==0)bReversedImage=true;
    else bReversedImage=false;

    char* text;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &text) > 0)
    {
        Description = text;
        if(Description.contains("ImageJ"))
        {
            if(Description.contains("hyperstack=true"))
            {
                TIFFError("\n ReadTiff","Cannot handle ImageJ hyperstacked images");
                close();
                return(false);
            }
        }
    }
    if(TIFFGetField(tif, TIFFTAG_SOFTWARE, &text) > 0)
        Software = text;

    if(TIFFGetField(tif, TIFFTAG_DATETIME, &text) > 0)
        DateTime = text;

    nFmin = new uint16_t[nSamplesPerPixel];
    nFmax = new uint16_t[nSamplesPerPixel];
    if(TIFFGetField(tif, TIFFTAG_MINSAMPLEVALUE, nFmin) > 0)
        bMinMax = true;
    else
        bMinMax = false;

    if(TIFFGetField(tif, TIFFTAG_MAXSAMPLEVALUE, nFmax) > 0)
        bMinMax = true;
    else
        bMinMax = false;

    if(TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &nResolutionUnit) > 0)
    {
        fXResolution=GetTIFFfloat(TIFFTAG_XRESOLUTION);
        fYResolution=GetTIFFfloat(TIFFTAG_YRESOLUTION);
        if(nResolutionUnit > 1)
        {
            if(fXResolution != 0.0f)
            {
                fXPixelSize = 1.e7f/fXResolution; // from dots per cm to nanometres
                if(nResolutionUnit == 2)fXPixelSize /= 2.54f;
            }
            if(fYResolution != 0.0f)
            {
                fYPixelSize = 1.e7f/fYResolution; // from dots per cm to nanometres
                if(nResolutionUnit == 2)fYPixelSize /= 2.54f;
            }
        }
    }

    TIFFGetFieldDefaulted(tif,TIFFTAG_PLANARCONFIG , &nPlanarConfiguration);
    if(nPlanarConfiguration > 1)
    {
        TIFFError("\n ReadTiff","Can only handle Planar Configuration = 1");
        close();
        return(false);
    }

    float X0 = GetTIFFfloat(TIFFTAG_XPOSITION);
    if(X0 < 0.0f)
        nX0 = 0;
    else
        nX0 = int(X0);
    float Y0 = GetTIFFfloat(TIFFTAG_YPOSITION);
    if(Y0 < 0.0f)
        nY0 = 0;
    else
        nY0 = int(Y0);

    if(!CheckDirectories())return false;

    if(bReadHeaderOnly)
    {
        TIFFClose(tif);
        return true;
    }

    if(nNoZPlanes == 1)
        cerr << "  There is 1 image in this file\n";
    else
        cerr << "  There are "<< nNoZPlanes << " images in this file\n";

    if(zstart == 0)
        zstart = 1;
    if(zend == 0)
        zend = nNoZPlanes;

    bool bZErr = false;
    if(zstart < 0)
    {
        cerr << " Error - Z start must be positive\n";
        bZErr = true;
    }

    if(zstart > zend)
    {
        cerr << " Error - Z end (" << zend << ") must be >= Z start (" << zstart << ")\n";
        bZErr = true;
    }

    if(uint32_t(zend) > nNoZPlanes)
    {
        cerr << " Warning: zend set to the number of planes (" << nNoZPlanes << ")\n";
        zend = nNoZPlanes;
    }

    if(bZErr)
    {
        close();
        return(false);
    }

    uint nPlanesRequired = uint((zend-zstart+1)/zinc);
    nPixelsPerPlane=nXdim*nYdim;
    nPixelsPerImage=nPixelsPerPlane*nPlanesRequired;

    nBytesPerPixel = nByteMultiplier*nSamplesPerPixel;
    nBytesPerPlane = nBytesPerPixel*nPixelsPerPlane;
    nBytesPerImage = nPlanesRequired*nBytesPerPlane;

    uint32_t nBytesPerLine = nBytesPerPixel*nXdim;

    ByteData = new uint8_t [nBytesPerImage];
    if(ByteData == NULL)
    {
        TIFFError("\n ReadTiff","Error - Out of memory for ByteData");
        close();
        return(false);
    }

    // read data

    uint8_t* buffer = 0;
    uint32_t buflen;

    int nDirNo = zstart-1;
    
    if(TIFFIsTiled(tif))
    {
        uint32_t tileWidth, tileLength;

        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
        TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
        buflen = TIFFTileSize(tif);

        int tilesAcross = (nXdim + tileWidth - 1) / tileWidth;
        int tilesDown = (nYdim + tileLength - 1) / tileLength;

        buffer = new uint8_t[buflen];

        uint32_t nBytesinTileRow = TIFFTileRowSize(tif);

        for(uint jk=0; jk < nPlanesRequired; jk++)
        {
            if(!TIFFSetDirectory(tif, nDirNo))
            {
                cerr << "\nCannot open directory for plane " << (nDirNo + 1) << "\n\n";
                close();
                return(false);
            }
            // Tiled TIFF
            uint8_t* pbplanest = ByteData + jk * nBytesPerPlane;
            for (int row = 0; row < tilesDown; row ++)
            {
                uint32_t ypos = row * tileLength;
                for (int col = 0; col < tilesAcross; col ++)
                {
                    uint32_t xpos = col * tileWidth;

                    if (TIFFReadTile(tif, buffer, xpos, ypos, 0, 0) < 0)
                    {
                        ErrInfo = QString("Error while reading tile at %1,%2").arg(xpos,ypos);
                        TIFFError("\n ReadTiff", qPrintable(ErrInfo));
                        delete [] buffer;
                        close();
                        return(false);
                    }

                    // copy from buffer into image array

                    uint8_t *pbuf = buffer;
                    uint nNoLines = tileLength;
                    uint nNoBytes =  nBytesinTileRow;

                    if (xpos + tileWidth > uint32_t(nXdim))
                        nNoBytes = (nXdim - xpos)*nBytesPerPixel;
                    if (ypos + tileLength > uint32_t(nYdim))
                        nNoLines = nYdim - ypos;

                    uint pbOffset = nBytesPerLine+nNoBytes;
                    uint pbufOffset = nBytesinTileRow-nNoBytes;

                    //flip image: copy from bottom to top
                    uint8_t* pb = pbplanest - (ypos + 1)*nBytesPerLine + xpos*nBytesPerPixel;

                    for(uint j = 0; j < nNoLines; j++)
                    {
                        for(uint i = 0; i < nNoBytes; i++)
                            *pb++ = *pbuf++;

                        pbuf += pbufOffset;
                        pb -= pbOffset;
                    }
                } // end for col..
            }  // end for row..
            nDirNo += zinc;
        }  // end for jk ..
    }
    else
    {
        buflen = TIFFScanlineSize(tif);
        buffer = new uint8_t[buflen];
        int lines=nBytesPerPlane/buflen;
        uint ipos  = 0;
        for(uint jk=1; jk <= nPlanesRequired; jk++)
        {
            if(!TIFFSetDirectory(tif, nDirNo))
            {
                cerr << "\nCannot open directory for plane " << (nDirNo + 1) << "\n\n";
                close();
                return(false);
            }
            // Line TIFF
            ipos = jk*nBytesPerPlane;
            for (int i = 0; i < lines; i++)
            {
                TIFFReadScanline(tif, buffer, i, 0);
                // flip image
                for(uint32_t k = 0; k < buflen; k++)
                    ByteData[ipos - (i+1) * buflen + k] = buffer[k];
            }
            nDirNo += zinc;
        }
    }
    delete [] buffer;

    if(bRGB && cChannel != ' ')
    {
        int nOffset = 0;
        if(cChannel == 'g')
            nOffset=1*nByteMultiplier;
        if(cChannel == 'b')
            nOffset=2*nByteMultiplier;

        tempData = new uint8_t [nBytesPerImage/3];
        uint8_t* pB = ByteData + nOffset;
        uint8_t* pBN = tempData;
        for(uint j=0; j < nPixelsPerImage; j++)
        {
            *pBN++ = *pB++;
            if(!bByteImage)
                *pBN++ = *pB++;
            pB += 2*nByteMultiplier;
        }
        nBytesPerPixel /= 3;
        nBytesPerImage /= 3;
        nBytesPerPlane /= 3;
        delete [] ByteData;
        ByteData = tempData;
        tempData = nullptr;

        bRGB = false;
        nSamplesPerPixel = 1;
        nPhotometric = 1;
    }
    nNoZPlanes = nPlanesRequired;
    return true;
}
void TiffImage::close()
{
    if(tif != nullptr)
       TIFFClose(tif);
}
bool TiffImage::CheckDirectories()
{
    int nNoDirs = 1;
    while (!TIFFLastDirectory(tif))
    {
        if(!TIFFReadDirectory(tif))
        {
            TIFFError("\n ReadTiff","Error while reading directory\n");
            return false;
        }
        nNoDirs++;
        int width=GetTIFFint(TIFFTAG_IMAGEWIDTH);
        int height=GetTIFFint(TIFFTAG_IMAGELENGTH);
        if(uint32_t(width) != nXdim || uint32_t(height) != nYdim)
        {
            cerr << " Error - width and height " << width << " x " << height;
            cerr << " in directory " << nNoDirs << " does not match that ";
            cerr << " in initial directory (" << nXdim << " x " << nYdim << ")\n";
            return false;
        }

        int bitspersample=GetTIFFshort(TIFFTAG_BITSPERSAMPLE);
        if(nBitsPerSample != bitspersample)
        {
            ErrInfo =QString("Bits per sample different in image %1").arg(nNoDirs);
            TIFFError("\n ReadTiff",qPrintable(ErrInfo));
            return false;
        }
    }

    TIFFSetDirectory(tif,0);

    nNoZPlanes = nNoDirs;
    return true;
}
float TiffImage::GetTIFFfloat(ttag_t tag)
{
    float val;
    if(TIFFGetField(tif, tag, &val) == 0)
        return 0.0;
    return (val);
}	
short TiffImage::GetTIFFshort(ttag_t tag)
{
    uint16_t val;
    if(TIFFGetField(tif, tag, &val) == 0)
        return short(0);
    return short(val);
}
int TiffImage::GetTIFFint(ttag_t tag)
{
    uint32_t val;
    if(TIFFGetField(tif, tag, &val) == 0)
        return 0;
    return int(val);
}
void TiffImage::setDescription(const char* szDescription)
{
    QString t = szDescription;
    if(t.length() > 0)
        Description = t;
}
void TiffImage::addDescription(const char* szLine)
{
    Description.append(szLine);
}
void TiffImage::setSoftware(const char* szSoftware)
{
    QString t = szSoftware;
    if(t.length() > 0)
        Software = t;
}
void TiffImage::setDateTime(const char* szDateTime)
{
    QString t = szDateTime;
    if(t.length() == 19 )
        DateTime = t;
}
void TiffImage::setDateTime(QString DateandTime)
{
    DateTime = DateandTime;
}
void TiffImage::setRes(float fPix, uint16_t nResUnit)
{
    nResolutionUnit = nResUnit;
    if(fPix > 0)
    {
        fXResolution = fPix;
        fYResolution = fPix;
    }
}
void TiffImage::setBitsPerSample(uint16_t bps)
{
    nBitsPerSample = bps;
    if(bps == 8)
        bByteImage=true;
    else if(bps == 16)
    {
        bByteImage=false;
        nByteMultiplier = 2;
    }
    else
        cerr << " Error - can\'t set bits per sample to " << bps << "\n";
}     
void TiffImage::setColourMap(uint16_t* r, uint16_t* g, uint16_t* b)
{
    delete [] red;
    delete [] green;
    delete [] blue;

    red = r;
    green = g;
    blue = b;
}     
bool TiffImage::setOutputType(int nOutput, char cChannelL)
{
    int nOutputType = nOutput;
    bByteImage = true;
    nBitsPerSample=8;
    nByteMultiplier = 1;
    nSamplesPerPixel = 1;
    nPhotometric = 1;
    bRGB = false;
    if(nOutputType >  TIFF_RGBA_8bit)
    {
        bByteImage=false;
        nBitsPerSample = 16;
        nByteMultiplier = 2;
    }
    if(nOutputType == TIFF_RGB_8bit || nOutputType == TIFF_RGB_16bit)
    {
        bRGB=true;
        nPhotometric = 2;
        nSamplesPerPixel = 3;
        cChannel = cChannelL;
    }
    if(nOutputType == TIFF_RGBA_8bit || nOutputType == TIFF_RGBA_16bit)
    {
        bRGB=true;
        nPhotometric = 2;
        nSamplesPerPixel = 4;
        cChannel = cChannelL;
    }
    if(nOutputType == TIFF_Map_8bit || nOutputType == TIFF_Map_16bit)
        nPhotometric = 3;

    if(nOutputType < TIFF_Mono_8bit || nOutputType >  TIFF_Map_16bit)
    {
        TIFFError("\n ReadTiff","Invalid output type chosen");
        return false;
    }
    nBytesPerPixel = nSamplesPerPixel*nByteMultiplier;
    return true;
}
bool TiffImage::Segment(int xstart, int xend, int ystart, int yend)
{
    if(xend == 0 || uint32_t(xend) > nXdim)xend=nXdim;
    if(yend == 0 || uint32_t(yend) > nYdim)yend=nYdim;

    bool bXlim=(xstart < 1 ||  xstart  > xend);
    bool bYlim=(ystart < 1 ||  ystart  > yend);
    if(bXlim || bYlim)
    {
        ErrInfo =QString("Invalid segment values: X = %1 %2; Y = %3 %4").arg(xstart).arg(xend).arg(ystart).arg(yend);
        return false;
    }

    // Check if segmentation is necessary

    bool XFull=(xstart == 1 && uint32_t(xend) == nXdim);
    bool YFull=(ystart == 1 && uint32_t(yend) == nYdim);
    if(XFull && YFull)
    {
        CalculateMinMax();
        return true;
    }
    uint32_t  nOldXdim = nXdim;
    uint  nOldBytesPerPlane=nBytesPerPlane;

    tempData = ByteData;
    ByteData = nullptr;

    uint32_t newXdim = uint32_t(xend - xstart + 1);
    uint32_t newYdim = uint32_t(yend - ystart + 1);
    int newX0  = int(nX0 + xstart - 1);
    int newY0  = int(nY0 + ystart - 1);

    setDim(newXdim, newYdim, nNoZPlanes);
    setOrigin(newX0, newY0);

    uint8_t* pBNew = ByteData;
    uint8_t* pBStart = tempData + ((ystart-1) * nOldXdim + xstart - 1)* nBytesPerPixel;

    for(uint32_t nPlaneNo = 0; nPlaneNo < nNoZPlanes; nPlaneNo++)
    {
        uint8_t* pBOld = pBStart + (nPlaneNo * nOldBytesPerPlane);
        for(int iy = ystart-1; iy < yend; iy++)
        {
            for (int ix = xstart-1; ix < xend; ix++)
            {
                for(uint ib = 0; ib < nBytesPerPixel; ib++)
                    *pBNew++=*pBOld++;
            }
            pBOld += (nOldXdim - xend + xstart - 1)*nBytesPerPixel;
        }
    }

    CalculateMinMax();
    delete [] tempData;
    tempData = nullptr;
    return true;
}
void TiffImage::CalculateMinMax()
{
    if(nFmin == nullptr)
        nFmin = new uint16_t[nSamplesPerPixel];
    if(nFmax == nullptr)
        nFmax = new uint16_t[nSamplesPerPixel];
    for(int j = 0; j < nSamplesPerPixel; j++)
    {
        nFmin[j] = 65535;
        nFmax[j] = 0;
    }

    if(nBitsPerSample == 8)
    {
        uint8_t* pBy = ByteData;
        for(uint j = 0; j < nPixelsPerImage; j++)
        {
            for(int k = 0; k < nSamplesPerPixel; k++)
            {
                nFmin[k] = (nFmin[k] > *pBy ? *pBy : nFmin[k]);
                nFmax[k] = (nFmax[k] < *pBy ? *pBy : nFmax[k]);
                pBy++;
            }
        }
    }
    else
    {
        uint16_t* pUS = reinterpret_cast<uint16_t*> (ByteData);
        for(uint j = 0; j < nPixelsPerImage; j++)
        {
            for(int k = 0; k < nSamplesPerPixel; k++)
            {
                nFmin[k] = (nFmin[k] > *pUS ? *pUS : nFmin[k]);
                nFmax[k] = (nFmax[k] < *pUS ? *pUS : nFmax[k]);
                pUS++;
            }
        }
    }
    for(int j = 0; j < nSamplesPerPixel; j++)
        std::cerr << "Channel " << (j+1) << " Min = " << nFmin[j] << "; Max = " << nFmax[j] << endl;
}
bool TiffImage::WriteTiff(QString FileName)
{
    if(ByteData == nullptr)
    {
        TIFFError("WriteTiff","No data to write!");
        return false;
    }

        tif = TIFFOpen(qPrintable(FileName),"w");
        if(!tif)
            return false;

    if(!nFmin)
    {
        CalculateMinMax();
    }
    switch(nSamplesPerPixel)
    {
      case 1:
        TIFFSetField(tif,TIFFTAG_MINSAMPLEVALUE,  nFmin[0]);
        TIFFSetField(tif,TIFFTAG_MAXSAMPLEVALUE,  nFmax[0]);
        break;
      case 3:
        TIFFSetField(tif,TIFFTAG_MINSAMPLEVALUE, nFmin[0], nFmin[1], nFmin[2]);
        TIFFSetField(tif,TIFFTAG_MAXSAMPLEVALUE, nFmax[0], nFmax[1], nFmax[2]);
        break;
      case 4:
        TIFFSetField(tif,TIFFTAG_MINSAMPLEVALUE, nFmin[0], nFmin[1], nFmin[2], nFmin[3]);
        TIFFSetField(tif,TIFFTAG_MAXSAMPLEVALUE, nFmax[0], nFmax[1], nFmax[2], nFmin[3]);
        break;
    }

    if(fXResolution > 0.)
    {
        TIFFSetField(tif, TIFFTAG_XRESOLUTION, fXResolution,1);
        TIFFSetField(tif, TIFFTAG_YRESOLUTION, fYResolution,1);
    }
    if(Description.length() != 0)
        TIFFSetField(tif,TIFFTAG_IMAGEDESCRIPTION,qPrintable(Description));
    if(Software.length() != 0)
        TIFFSetField(tif,TIFFTAG_SOFTWARE,qPrintable(Software));

    uint nBytesPerLine = nXdim * nBytesPerPixel;

    for(uint32_t z=0; z < nNoZPlanes; z++)
    {
        TIFFSetField(tif,TIFFTAG_IMAGEWIDTH,nXdim);
        TIFFSetField(tif,TIFFTAG_IMAGELENGTH,nYdim);

        TIFFSetField(tif,TIFFTAG_SAMPLESPERPIXEL, nSamplesPerPixel);
        if(nSamplesPerPixel == 4)
        {
            const uint16_t extras[] = {EXTRASAMPLE_ASSOCALPHA};
            TIFFSetField(tif,TIFFTAG_EXTRASAMPLES, uint16_t(1), extras);
        }
        TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, nBitsPerSample);
        TIFFSetField(tif,TIFFTAG_PLANARCONFIG, uint16_t(1));
        TIFFSetField(tif,TIFFTAG_COMPRESSION, uint16_t(1));  // no compression
        TIFFSetField(tif,TIFFTAG_PHOTOMETRIC, nPhotometric);
        if(nPhotometric == 3)
            TIFFSetField(tif,TIFFTAG_COLORMAP, red, green, blue);

        TIFFSetField(tif,TIFFTAG_ORIENTATION, nOrientation);

        if(nResolutionUnit > 1)
        {
            TIFFSetField(tif,TIFFTAG_RESOLUTIONUNIT, nResolutionUnit);
            TIFFSetField(tif,TIFFTAG_XRESOLUTION, fXResolution);
            TIFFSetField(tif,TIFFTAG_YRESOLUTION, fYResolution);
        }
        TIFFSetField(tif,TIFFTAG_XPOSITION, float(nX0));
        TIFFSetField(tif,TIFFTAG_YPOSITION, float(nY0));

        if(DateTime.length() != 19)
        {
            QDateTime t;
            t.currentDateTime();
            DateTime = t.toString("yyyy:MM:dd hh:mm:ss");
        }
        TIFFSetField(tif,TIFFTAG_DATETIME, qPrintable(DateTime));

        if(nNoZPlanes > 1)
        {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER , z+1, nNoZPlanes);
        }

        uint8_t* pbData = ByteData + z * nBytesPerPlane;
        for(uint32_t k = 0; k < nYdim; k++)
        {
            // flip image
            uint8_t* pbL = pbData + (nYdim - 1 - k)*nBytesPerLine;
            TIFFWriteScanline(tif, pbL, k, 0);
        }
        if(nNoZPlanes > 1)TIFFWriteDirectory(tif);
    }
    TIFFClose(tif);
    return true;
}
