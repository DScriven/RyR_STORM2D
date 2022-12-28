/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, blinkroi.cpp, is part of the RyR_STORM2D program and generates an interface for determining
* an ROI.
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
#include "blinkroi.h"
#include "ui_blinkroi.h"
#include "stormdensity.h"
#include <QFile>
#include <QString>

BlinkROI::BlinkROI(QOpenGLWindow *parent)
{
    setupUi(this);
    connect(this,SIGNAL(writeROI(QString, QRectF)), parent, SIGNAL(writeROI(QString, QRectF)));
    connect(this,SIGNAL(hideROI()), parent, SLOT(hideROI()));
    connect(this, SIGNAL(AlertMsg(QString, char)), parent, SIGNAL(AlertMsg(QString,char)));
    connect(pbWriteROI, SIGNAL(clicked()), this, SLOT(OnWriteROI()));
    connect(pbCancel, SIGNAL(clicked()), this, SLOT(OnCancel()));
}
BlinkROI::~BlinkROI()
{
}
void BlinkROI::keyPressEvent( QKeyEvent *e )
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
void BlinkROI::displayMenu(QRectF roinm, QString FileName)
{
    float xmin = roinm.left();
    float xmax = roinm.right();
    if(xmin > xmax)
    {
        float temp = xmin;
        xmin = xmax;
        xmax = temp;
    }
    float ymin = roinm.top();
    float ymax = roinm.bottom();
    if(ymin > ymax)
    {
        float temp = ymin;
        ymin = ymax;
        ymax = temp;
    }
    xmin = (xmin < 0 ? 0 : xmin);
    ymin = (ymin < 0 ? 0 : ymin);
    leMinX->setText(QString::number(xmin,'f',1));
    leMinX->setReadOnly(false);
    leMaxX->setText(QString::number(xmax,'f',1));
    leMaxX->setReadOnly(false);
    leMinY->setText(QString::number(ymin,'f',1));
    leMinY->setReadOnly(false);
    leMaxY->setText(QString::number(ymax,'f',1));
    leMaxY->setReadOnly(false);
    leFileName->setText(CheckOutputFileName(FileName));
    move(10,100);
    show();
}
QString BlinkROI::CheckOutputFileName(QString NewFileName)
{
    QString OutputFile;
    QString FormatString;
    QFile qdata;

    NewFileName.remove(".3dSTM");
    int segno = 1;
    int spos = NewFileName.indexOf("-seg-");
    if(spos > -1)
    {
        QString Fnseg = NewFileName.mid(spos+5);
        int spos2=Fnseg.indexOf("-");
        if(spos2 > 0)
        {
            NewFileName=NewFileName.mid(0, spos+5+spos2);
            Fnseg = Fnseg.mid(0,spos2);
            segno = Fnseg.toInt();
        }

        FormatString = "%1-%2.3dSTM";
    }
    else
        FormatString = "%1-seg-%2.3dSTM";

    bool bMismatch = false;
    do{
        bMismatch = false;
        forever
        {
            OutputFile = QString(FormatString).arg(NewFileName).arg(segno);
            qdata.setFileName(OutputFile);
            if(!qdata.exists())
                break;
            segno++;
        }
    } while(bMismatch);

    return OutputFile;
}
void BlinkROI::OnWriteROI()
{
    QRectF ROIXY;
    ROIXY.setLeft(leMinX->text().toFloat());
    ROIXY.setRight(leMaxX->text().toFloat());
    ROIXY.setBottom(leMinY->text().toFloat());
    ROIXY.setTop(leMaxY->text().toFloat());
    FileName = leFileName->text();


#ifdef _DEBUG_
    emit AlertMsg(QString("Xst = %1; Yst = %2; Xend = %3; Yend = %4").arg(ROIXY.Left()).arg(ROIXY.Bottom()).arg(ROIXY).arg(ROIXY.Top()),'g');
#endif
    emit writeROI(FileName, ROIXY);
    emit hideROI();
    close();
}
void BlinkROI::OnCancel()
{
    emit hideROI();
    close();
}
