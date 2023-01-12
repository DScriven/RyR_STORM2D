/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, param.h, is a header file for param.cpp and is part of the RyR_STORM2D program.
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
#ifndef __param_h
#define __param_h

#include "ui_param.h"
#include <QtWidgets/QDialog>
#include <qstring.h>
#include "paramvals.h"

class QScrollBar;
class QSettings;
class QThread;
class STORMdensity;
class ViewClusters;

class ParamWnd : public QDialog, private Ui::Parameters
{
   Q_OBJECT

   QThread* thread1;
   QString STORMFile;
   QString InputFile;
   QString DirectoryName;
   QScrollBar* sb;
   QSettings* imageVal;
   STORMdensity* STORMcalc;
   ViewClusters* Image;

   ParamVals ParameterValues;

   int nMinFrame;
   int nMaxFrame;
   int nFrameMin;
   int nFrameMax;

   void keyPressEvent(QKeyEvent *e);

   bool bPlotPoint;
   bool bLoaded;
   bool bLoadedFile;

   bool getFilesAndValues();
   bool getValues();

public:
   ParamWnd(QSize Screensize, qreal Scrphysdpi, qreal Scrlogdpi);

public slots:
   void showHelp();
   void setMsg(QString msg);
   void setAlert(QString msg, char color);
   void updateProgressBar(int value);
   void setFrameRange(int MinFrame, int MaxFrame);

private slots:
   void CloseDown();
   void FilePressed();
   void LoadPressed();
   void ModifyPressed();
   void changeFile();
   void OnChangeThreshold();

signals:
   void calculate_clusters(QString, ParamVals);
   void	apply_new_params(ParamVals, int, int);
};
#endif
