/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, blinkroi.h, is a header file for blinkroi.cpp as is part of the RyR_STORM2D program.
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
#ifndef BLINKROI_H
#define BLINKROI_H

#include <QDialog>
#include <QOpenGLWindow>
#include <QRectF>
#include <QKeyEvent>
#include "ui_blinkroi.h"

class BlinkROI : public QDialog, private Ui::BlinkRegion
{
    Q_OBJECT

public:
    BlinkROI(QOpenGLWindow* parent);
    ~BlinkROI();
    void displayMenu(QRectF roinm,  QString Filename);

public slots:
    void OnCancel();
private slots:
    void OnWriteROI();

private:
    bool bNotShown;

    QString InputFileName;
    QString FileName;
    QString CheckOutputFileName(QString NewFileName);
    void keyPressEvent( QKeyEvent *e );
signals:
    void writeROI(QString FileName, QRectF ROIXY);
    void hideROI();
    void AlertMsg(QString, char);
    void CloseProgram();
};

#endif // BLINKROI_H
