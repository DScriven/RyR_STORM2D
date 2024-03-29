/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022, 2023.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, main.cpp, is the entry file for the RyR_STORM2D program.
*
* The following files are part of the program:
*        blinks.h, blinkvertex.h, tetramer.h, vertex_data_traits.h
*        blinkroi.cpp & blinkroi.h
*        cluster.cpp & cluster.h
*        param.cpp & param.h
*        stormdensity.cpp & stormdensity.h
*        setDPI.cpp &setDPI.h
*        tiff6qt.cpp & tiff6.h
*        viewclusters.cpp & viewclusters.h
*
* This program can be compiled for Windows 10 and 11 or a Linux operating system.
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
#ifdef WIN32
#include <QCoreApplication>
#else
#include "unistd.h"
#endif

#include <QGuiApplication>
#include <QSurfaceFormat>
#include <QScreen>
#include "param.h"
#include <iostream>


int main(int argc, char *argv[])
{
#ifndef WIN32
    int status = fork();
    switch (status)
    {
      case -1:
      {
        std::cerr << "\nError - could not fork process\n\n";
        exit(1);
      }
      case 0: // child process
        break;
      default: /// parent process
        exit(0);
    }
    status = setsid();
    if (status == -1)
    {
        std::cerr << "\nError - could not get valid sid\n\n";
        exit(1);
    }
#else
    QCoreApplication::setApplicationName("RyR_STORM2D.exe");
    QCoreApplication::setApplicationVersion("1.3");
    QCoreApplication::setOrganizationName("Moore Lab, Life Sciences Inst., Univ. British Columbia");
    QCoreApplication::setAttribute(Qt::AA_UseDesktopOpenGL);
#endif
    QApplication app(argc, argv);
    QScreen* scrn = QApplication::primaryScreen();
    qreal Screenlogicdpi = scrn->logicalDotsPerInch();
    qreal Screenphysicaldpi = qRound(scrn->physicalDotsPerInch());
    QSize ScreenSize = scrn->availableSize();
    ParamWnd w(ScreenSize, Screenphysicaldpi,Screenlogicdpi);
    w.show();

    return app.exec();
}

