/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, setDPI.cpp, is part of the RyR_STORM2D program and is used to  determine the magnification
* and DPI of a generated hi-res TIFF image.
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
#include "setDPI.h"
#include "ui_setDPI.h"
#include <QString>
#include <QMessageBox>

TiffMagDPI::TiffMagDPI(QOpenGLWindow *parent)
{
    setupUi(this);
    connect(this, SIGNAL(AlertMsg(QString, char)), parent, SIGNAL(AlertMsg(QString,char)));
    connect(this, SIGNAL(FBOTiff(double, int)), parent, SLOT(FBOTiff(double,int)));
    connect(pbCreateTIFF, SIGNAL(clicked()), this, SLOT(OnCreateTiff()));
    connect(pbCancel, SIGNAL(clicked()), this, SLOT(OnCancel()));
    connect(leMagnify, SIGNAL(editingFinished()),this, SLOT(CalculateImageSize()));
    connect(leDPI, SIGNAL(editingFinished()),this, SLOT(CalculateImageSize()));
    connect(cmbUnit, SIGNAL(currentIndexChanged(int)),this, SLOT(CalculateImageSize()));
    ImageDPI = 300;
    FBOMag = 1.0;
}
TiffMagDPI::~TiffMagDPI()
{
}
void TiffMagDPI::keyPressEvent( QKeyEvent *e )
{
    switch( e->key() )
    {
      case Qt::Key_Escape:
        OnCancel();
        break;
      default:
        break;
    }
}
void TiffMagDPI::displayMenu(int XInitPixelSize, int YInitPixelSize)
{
    InitXPixels = XInitPixelSize;
    InitYPixels = YInitPixelSize;
    leMagnify->setText(QString::number(FBOMag));
    leDPI->setText(QString::number(ImageDPI));
    leWidth->clear();
    leWidth->setReadOnly(true);
    leHeight->clear();
    leHeight->setReadOnly(true);
    cmbUnit->setCurrentIndex(0);
    move(600,100);
    show();
}
void TiffMagDPI::CalculateImageSize()
{
    FBOMag = leMagnify->text().toDouble();
    if(FBOMag < 1.0 || FBOMag > 8.0)
    {
        QMessageBox::warning(this,"Magnification Error","Magnification outside range of 1 to 8",QMessageBox::Ok);
        return;
    }
    ImageDPI = leDPI->text().toInt();
    if(ImageDPI < 150 || ImageDPI > 1000)
    {
        QMessageBox::warning(this,"DPI Error","DPI outside range of 150 to 1000",QMessageBox::Ok);
        return;
    }
    double scale=FBOMag/double(ImageDPI);
    double physX = scale*double(InitXPixels);
    double physY = scale*double(InitYPixels);
    if(cmbUnit->currentIndex() == 1)
    {
      physX *= 2.54;
      physY *= 2.54;
    }
    leWidth->setText(QString::number(physX,'f',1));
    leHeight->setText(QString::number(physY,'f',1));
}
void TiffMagDPI::OnCreateTiff()
{
    emit FBOTiff(FBOMag, ImageDPI);
    close();
}
void TiffMagDPI::OnCancel()
{
    close();
}
