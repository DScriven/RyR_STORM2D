/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, param.cpp, is part of the RyR_STORM2D program and acts as an interface between the
* user and the main program.
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
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QFileDialog>
#include <QObject>
#include <QThread>
#include <QString>
#include <QFileInfo>
#include <QSettings>
#include <QRectF>
#include <QDesktopServices>
#include "param.h"
#include "stormdensity.h"
#include "viewclusters.h"
//#include <QLoggingCategory>

ParamWnd::ParamWnd(QSize Screensize, qreal Screendpi)
{
    qRegisterMetaType<ParamVals>("ParamVals");
    setupUi(this);
    setModal(false);
    setWindowIcon(QIcon(":/RyR_STORM2D.ico"));
    setWindowFlags(Qt::CustomizeWindowHint);

    DirectoryName ="";
    InputFile = "";
    Directory->setVisible(false);
    imageVal = new QSettings(QSettings::IniFormat, QSettings::UserScope,
                             "ChouMooreLabsUBC", "RyR_STORM2D");

    bLoadedFile = false;
    QString pName = imageVal->value("STORMFile", "").toString();
    STORMFile = pName;
    if(!pName.isEmpty())
    {
        QFileInfo q(pName);
        QString pVal = q.completeBaseName();
        Directory->setVisible(true);
        DirectoryName=q.path();
        Directory->setText(DirectoryName);
        lePosnFile->setText(pVal);
    }

    ParameterValues.MinTetramerArea = imageVal->value("MinTetramerArea",600).toDouble();
    ParameterValues.MinimumBlinksPerCluster = (imageVal->value("MinimumBlinksPerCluster", 6).toInt());
    ParameterValues.NeighbourhoodLimit  = imageVal->value("NeighbourhoodLimit",30).toDouble();
    ParameterValues.Alpha1stPass = imageVal->value("Alpha1stPass",700).toDouble();
    ParameterValues.LogDensityThreshold = imageVal->value("LogDensityThreshold",-3.0).toDouble();

    sb = teMessage->verticalScrollBar();
    leFrameStart->setEnabled(false);
    leFrameStop->setEnabled(false);
    ckDispNumbers->setChecked(false);
    ckDispOuter->setChecked(false);
    ckDispConvexHull->setChecked(false);
    ckDispMinEllipse->setChecked(false);
    ckSaveClusterData->setChecked(false);
    rbNone->setChecked(true);
    leLogDensityThresh->setText(QString().setNum(ParameterValues.LogDensityThreshold));
    leLogDensityThresh->hide();

    leAlpha1stPass->setText(QString().setNum(ParameterValues.Alpha1stPass));
    leNeighbourhoodLimit->setText(QString().setNum(ParameterValues.NeighbourhoodLimit));
    leMinTetramerArea->setText(QString().setNum(ParameterValues.MinTetramerArea));
    pbApply->hide();

    sbMinBlinksPerCluster->setValue(ParameterValues.MinimumBlinksPerCluster);

    thread1 = new QThread;
    STORMcalc = new STORMdensity;
    STORMcalc->moveToThread(thread1);
    thread1->start();
    STORMcalc->setScreenParam(Screensize, Screendpi);

    Image = new ViewClusters(STORMcalc);

    connect(STORMcalc, SIGNAL(AlertMsg(QString, char)), this, SLOT(setAlert(QString, char)));
    connect(STORMcalc, SIGNAL(Progress(int)), this, SLOT(updateProgressBar(int)));
    connect(STORMcalc, SIGNAL(setFrameRange(int, int)), this, SLOT(setFrameRange(int, int)));
    connect(STORMcalc, SIGNAL(ShowImage()), Image, SLOT(Display()));
    connect(STORMcalc, SIGNAL(RedrawImage(QString)), Image, SLOT(Redraw(QString)));
    connect(STORMcalc, SIGNAL(closeROI()), Image, SLOT(hideROI()));
    connect(Image, SIGNAL(showHelp()), this, SLOT(showHelp()));
    connect(Image, SIGNAL(writeROI(QString, QRectF)), STORMcalc, SLOT(writeROI(QString, QRectF)));

    connect(Image, SIGNAL(AlertMsg(QString, char)), this, SLOT(setAlert(QString, char)));
    connect(Image, SIGNAL(CloseProgram()), this, SLOT(CloseDown()));

    connect(rbNone, SIGNAL(toggled(bool)), this, SLOT(OnChangeThreshold()));
    connect(lePosnFile, SIGNAL(editingFinished()), this, SLOT(changeFile()));
    connect(pbLoadFiles, SIGNAL(pressed()), this, SLOT(LoadPressed()));
    connect(pbClose, SIGNAL(pressed()), this, SLOT(CloseDown()));
    connect(pbBrowseFiles, SIGNAL(pressed()), this, SLOT(FilePressed()));
    connect(pbApply, SIGNAL(pressed()), this, SLOT(ModifyPressed()));

    connect(this, SIGNAL(calculate_clusters(QString, ParamVals)), STORMcalc, SLOT(CalculateStructure(QString, ParamVals)));
    connect(this, SIGNAL(apply_new_params(ParamVals, int, int)), STORMcalc, SLOT(NewClusterValues(ParamVals, int, int)));

    sb = teMessage->verticalScrollBar();
    progressBar->setRange(0,100);
    progressBar->hide();
    pbApply->hide();
    move(0,15);
//    QLoggingCategory::setFilterRules("*.debug=false\n"
//                                      "qt.qpa.gl=true");
}
void ParamWnd::CloseDown()
{
    if(bLoadedFile)
       Image->close();
    QString pfile= QString("STORMFile");
    QString SaveFile;
    if(!InputFile.isEmpty() && InputFile != STORMFile)
        SaveFile = InputFile;
    else
        SaveFile = STORMFile;

    imageVal->setValue(pfile, SaveFile);
    imageVal->setValue("MinTetramerArea",ParameterValues.MinTetramerArea);
    imageVal->setValue("MinimumBlinksPerCluster",ParameterValues.MinimumBlinksPerCluster );
    imageVal->setValue("NeighbourhoodLimit",ParameterValues.NeighbourhoodLimit);
    imageVal->setValue("Alpha1stPass",ParameterValues.Alpha1stPass);
    imageVal->setValue("LogDensityThreshold",ParameterValues.LogDensityThreshold);
    delete Image;
    delete imageVal;
    delete STORMcalc;
    close();
}
void ParamWnd::keyPressEvent(QKeyEvent *e)
{
    switch(e->key())
    {
    case Qt::Key_Escape:
        CloseDown();
        break;
    case Qt::Key_F1:
        showHelp();
        break;
    }
}
void ParamWnd::showHelp()
{
    QString InfoFiles = "file:///" % QDir::currentPath() % "/docs/index.html";
    QDesktopServices::openUrl(QUrl(InfoFiles, QUrl::TolerantMode));
}
void ParamWnd::FilePressed()
{
    QString ImageFilename =  QFileDialog::getOpenFileName( this,"Open STORM file", DirectoryName, "STORM files (*.3dSTM)");
    QFileInfo q(ImageFilename);
    DirectoryName=q.path();
    QString DirectoryMinName=DirectoryName;
    int pos = qMax(0,DirectoryMinName.indexOf("My Documents"));
    Directory->setText(DirectoryName.mid(pos));
    lePosnFile->setText(q.completeBaseName());
    changeFile();
}
void ParamWnd::changeFile()
{
    if(bLoadedFile)
    {
       pbApply->hide();
       pbLoadFiles->show();
    }
}
void ParamWnd::OnChangeThreshold()
{
    if(rbNone->isChecked())
       leLogDensityThresh->hide();
    else
       leLogDensityThresh->show();
}
void ParamWnd::setMsg(QString msg)
{
    teMessage->append(msg);
    sb->setValue(sb->maximum());
}
void ParamWnd::updateProgressBar(int value)
{
    if(value == -1)
    {
        progressBar->hide();
        pbApply->show();
    }
    else
    {
        if(!progressBar->isVisible())
        {
            progressBar->show();
            pbApply->hide();
        }
        progressBar->setValue(value);
    }
}
void ParamWnd::setAlert(QString msg, char color)
{
    if(color == 'c')
        msg = "<font color=cyan>" + msg + "</font>";
    if(color == 'g')
        msg = "<font color=green>" + msg + "</font>";
    if(color == 'o')
        msg = "<font color=orange>" + msg + "</font>";
    if(color == 'r')
        msg = "<font color=red>" + msg + "</font>";
    if(color == 'b')
        msg = "<font color=blue>" + msg + "</font>";
    if(color == 'm')
        msg = "<font color=magenta>" + msg + "</font>";

    setMsg(msg);
}
void ParamWnd::LoadPressed()
{
    bLoadedFile = true;
    if(getFilesAndValues())
    {
        setAlert("****************************************************************",'m');
        setAlert(QString(" File = %1").arg(InputFile),'m');
        emit calculate_clusters(InputFile, ParameterValues);
        pbLoadFiles->hide();
    }
}
void ParamWnd::ModifyPressed()
{
    if(getValues())
        emit apply_new_params(ParameterValues, nFrameMin, nFrameMax);
}
bool ParamWnd::getFilesAndValues()
{
    if(lePosnFile->text().isEmpty())
    {
        setAlert("Error - You must enter a 3dSTM file",'r');
        return false;
    }
    QString File1=DirectoryName+"/"+lePosnFile->text()+".3dSTM";
    QFileInfo f1(File1);
    if(!f1.exists())
    {
        setAlert(QString("File %1 not found!").arg(File1),'r');
        return false;
    }
    InputFile = File1;
    return getValues();
}
bool ParamWnd::getValues()
{
    ParameterValues.paramGeometry = frameGeometry();
    double NeighbourhoodLimit = leNeighbourhoodLimit->text().toDouble();
    if((NeighbourhoodLimit < 30.0) ||  (NeighbourhoodLimit > 300.0))
    {
        setAlert(QString("Error - Cluster Distance = %1 nm - out of range").arg(NeighbourhoodLimit),'r');
        return false;
    }
    ParameterValues.NeighbourhoodLimit = NeighbourhoodLimit;

    if(leLogDensityThresh->isEnabled())
    {
        double LogDensityThreshold = leLogDensityThresh->text().toDouble();
        if((LogDensityThreshold < -7) || (LogDensityThreshold > 2))
        {
            setAlert(QString("Error - Log Density Threshold = %1 - out of range").arg(LogDensityThreshold),'r');
            return false;
        }
        ParameterValues.LogDensityThreshold = LogDensityThreshold;
    }

    double MinTetramerArea = leMinTetramerArea->text().toDouble();
    if((MinTetramerArea < 200) || (MinTetramerArea > 2000))
    {
        if(MinTetramerArea != 0.0)
        {
            setAlert(QString("Error - Minimum Tetramer Area = %1 nm - out of range").arg(MinTetramerArea),'r');
            return false;
        }
        else
            setAlert("No tetramer area exclusion",'r');
    }
    ParameterValues.MinTetramerArea = MinTetramerArea;

    double Alpha1stPass = leAlpha1stPass->text().toDouble();
    if((Alpha1stPass < 100) || (Alpha1stPass > 10000))
    {
        setAlert(QString("Error - 1st pass Alpha = %1 nm - out of range").arg(Alpha1stPass),'r');
        return false;
    }
    ParameterValues.Alpha1stPass = Alpha1stPass;

    if(leFrameStart->isEnabled())
    {
        nFrameMin = leFrameStart->text().toInt();
        nFrameMax = leFrameStop->text().toInt();
        if(nFrameMin < nMinFrame)
        {
            nFrameMin = nMinFrame;
            leFrameStart->setText(QString().setNum(nMinFrame));
        }
        if(nFrameMax > nMaxFrame)
        {
            nFrameMax = nMaxFrame;
            leFrameStop->setText(QString().setNum(nMaxFrame));
        }
        if(nFrameMin > nFrameMax)
        {
            leFrameStart->setText(QString().setNum(nMinFrame));
            leFrameStop->setText(QString().setNum(nMaxFrame));
            return false;
        }
    }

    ParameterValues.SaveClusterData=false;
    if(ckSaveClusterData->isChecked())
        ParameterValues.SaveClusterData=true;

    ParameterValues.DispNoThreshold=false;
    if(rbNone->isChecked())
        ParameterValues.DispNoThreshold=true;

    ParameterValues.DispExcludedData=false;
    if(ckDispExcludedData->isChecked())
        ParameterValues.DispExcludedData=true;

    ParameterValues.DispConvexHull=false;
    if(ckDispConvexHull->isChecked())
        ParameterValues.DispConvexHull=true;

    ParameterValues.DispMinEllipse=false;
    if(ckDispMinEllipse->isChecked())
        ParameterValues.DispMinEllipse=true;

    ParameterValues.DispNumbers=false;
    if(ckDispNumbers->isChecked())
        ParameterValues.DispNumbers=true;

    ParameterValues.DispBoundary=true;
    if(!ckDispOuter->isChecked())
        ParameterValues.DispBoundary=false;

    ParameterValues.MinimumBlinksPerCluster = sbMinBlinksPerCluster->value();

    return true;
}
void ParamWnd::setFrameRange(int MinFrame, int MaxFrame)
{
    nMinFrame = MinFrame;
    nMaxFrame = MaxFrame;
    leFrameStart->setEnabled(true);
    leFrameStart->setText(QString().setNum(nMinFrame));
    leFrameStop->setEnabled(true);
    leFrameStop->setText(QString().setNum(nMaxFrame));
}
