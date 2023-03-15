/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from ryanodine receptors on the
* surface of ventricular cardiomyocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, viewclusters.cpp, is part of the RyR_STORM2D program and uses OpenGL to display the data
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
#include "cluster.h"
#include "viewclusters.h"
#include "blinkroi.h"
#include "tiff6qt.h"
#include <algorithm>
#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QMoveEvent>
#include <QInputDialog>
#include <QSurfaceFormat>
#include <QOffscreenSurface>
#include <QOpenGLShaderProgram>
#include <QOpenGLFramebufferObject>
#include <QOpenGLVertexArrayObject>
#include <QMatrix4x4>
#include <QPainter>
#include <QOpenGLPaintDevice>
#include <QInputDialog>
#include <QRubberBand>
#include <QDateTime>

using namespace std;

ViewClusters::ViewClusters(STORMdensity* SCalc) : QOpenGLWindow()
{
    SD = SCalc;
    basefont = QFont("Arial",8,QFont::Bold);
    imageMenu = new QMenu(nullptr);
    setFlag(Qt::CustomizeWindowHint);
    setFlag(Qt::WindowTitleHint);
    setIcon(QIcon(":/RyR_STORM2D.ico"));
    roiinfo = nullptr;
    bd = nullptr;
    boundline = nullptr;
    hullpts = nullptr;
    ellipsepts = nullptr;
    tetrameroutline = nullptr;
    rubberBand = nullptr;
    roiinfo = nullptr;
    image_fbo = nullptr;

    vs_start = nullptr;
    vs_count = nullptr;
    bn_start = nullptr;
    bn_count = nullptr;
    hl_start = nullptr;
    hl_count = nullptr;
    el_start = nullptr;
    el_count = nullptr;
    tl_start = nullptr;
    tl_count = nullptr;

    PaintDevice = nullptr;

    colclust = nullptr;
    boundclust = nullptr;
    tetplot = nullptr;

    initValues();

    degtorad=4.0*atan(1.0)/180.0;    // some values we need
    tetramerwidth=27;
    sqrt2 = sqrt(2.0);
    tetramerdiag=tetramerwidth*sqrt2/2;
    fortyfivedeg = 45.0 * degtorad;
    nPalette = 2;

    fSizeFactor = 0.97f;
}
ViewClusters::~ViewClusters()
{
    std::vector<int>().swap(TetMappedCluster);
    std::vector<NNDval>().swap(NNDtet);
    ctree.clear();

    if(colclust != nullptr)
    {
        delete colclust;
        colclust = nullptr;
    }
    if(boundclust != nullptr)
    {
        delete boundclust;
        boundclust = nullptr;
    }
    if(tetplot != nullptr)
    {
        delete tetplot;
        tetplot = nullptr;
    }

    makeCurrent();
    vsp_vbo.destroy();
    bound_vbo.destroy();
    hull_vbo.destroy();
    ellipse_vbo.destroy();
    scalebar_vbo.destroy();
    tetramers_vbo.destroy();
    singletetramer_vbo.destroy();
    doneCurrent();

    delete PaintDevice;
    delete rubberBand;
    delete roiinfo;

    delete [] bd;
    delete [] boundline;
    delete [] hullpts;
    delete [] ellipsepts;
    delete [] tetrameroutline;

    delete [] vs_start;
    delete [] vs_count;
    delete [] bn_start;
    delete [] bn_count;
    delete [] hl_start;
    delete [] hl_count;
    delete [] el_start;
    delete [] el_count;
    delete [] tl_start;
    delete [] tl_count;
}
void ViewClusters::initValues()
{
    fZoom = 1.0f;
    fZoomDelta = 0.25f;
    PointSize = 1;
    ExPointSize = 2;
    ImageOffset = Point_2(0,0);
    xshift = 0.0f;
    yshift = 0.0f;
    bLinePlot = false;
    bLeftButtonDown = false;
    bScaleBar = false;
    bScaleBarPositionChosen = false;

    bROI = false;
    bSavedTetFile = false;
    bCalcNearestNeighbours = true;
    bModifyTetramer=false;
    bAddTetramer=false;
    bClassifyTetramer=false;
    bPlacingTetramer=false;
    bDisplayNND = false;
    bShowTetramer=false;
    bPlotTetramers=false;
    bIsolateTetramers = false;
    bDisplayClusterInfo = false;
    bLeftButtonDown = false;
    bDoubleClick = false;

    nNoTetramers = 0;
    nImagedTetramers = 0;
    nBeingModified = -1;
    nBeingClassified = -1;
    NNDlinecolor = Qt::black;
    NNDtextcolor = Qt::black;
}
void ViewClusters::Display()
{
    QFileInfo fi(SD->FileName);
    fi.makeAbsolute();
    setTitle(fi.baseName());
    fXScreenMax = SD->ScreenSize.width()-550;
    fYScreenMax = SD->ScreenSize.height()-80;
    initValues();

    if(SD->displayNoThreshold)
        ThresholdColourIndex = 7;
    else
        ThresholdColourIndex = 0;

    float maxx=float(SD->Limits.width());
    float maxy=float(SD->Limits.height());
    float fMinRatio = qMin(fXScreenMax/maxx,fYScreenMax/maxy);
    InitXImage = qRound(fMinRatio*maxx);
    InitYImage = qRound(fMinRatio*maxy);
    PixelSize =  SD->Limits.width()/(InitXImage);
    if(PixelSize < 6.0)
    {
       double dReduce = PixelSize/6.0;
       InitXImage *= dReduce;
       InitYImage *= dReduce;
       PixelSize /= dReduce;
    }
    emit AlertMsg(QString("Zoom = %1 - PixelSize = %2 nm").arg(fZoom).arg(PixelSize,0,'f',2),'m');
    emit AlertMsg(QString("Screen: Xmax = %1 - Ymax = %2 pixels").arg(fXScreenMax).arg(fYScreenMax),'m');
    emit AlertMsg(QString("Initial Image size: X = %1, Y = %2 pixels").arg(InitXImage).arg(InitYImage),'m');
    Centre = Point_2(SD->Limits.width()/2.0, SD->Limits.height()/2.0);

    XImageSize = InitXImage;
    YImageSize = InitYImage;
    setupTransforms();

    NoClusters=uint(SD->Clusters.size());
    setupClusterSearch();
    QRect pG = SD->ParamGeometry;
    QPoint winposn = QPoint(qMin(pG.x()+pG.width()+3,550),pG.y());
    if(SD->reassignClusters)
    {
      emit AlertMsg("Reassigning Cluster Labels",'m');
      reassignClusterNos();
      nLabelledClusters = ClusterCount();
    }

    ImageOffset=Point_2(0.,0.);
    if(isVisible())
    {
      bLinePlot = false;
      bScaleBar = false;
      bScaleBarPositionChosen = false;
      resize(InitXImage, InitYImage);
      Redraw("Points Outer Inner");
    }
    else
    {
      setGeometry(0,0,InitXImage, InitYImage);
      createMenus();
      show();
    }
    setFramePosition(winposn);
}
void ViewClusters::CreateScreenShot()
{
    CreateTIFF(XImageSize, YImageSize, BackBuffer);
}
void ViewClusters::CreateFBOImage()
{
    bool ok;
    fFBOZoom = QInputDialog::getDouble(nullptr, tr("FBO magnification"),
                                         tr("FBO Magnification:"), 1.5, 1.05, 3.0, 2, &ok, Qt::WindowFlags(), 0.05);
    if(!ok)
      return;
    QOpenGLContext* ctx = context();
    QOffscreenSurface* offsurf = new QOffscreenSurface();
    QSurfaceFormat format;
    format.setDepthBufferSize( 4 );
    format.setSamples(24);
    format.setVersion(4,3);
    format.setRenderableType(QSurfaceFormat::OpenGL);
    offsurf->setFormat(format);
    ctx->makeCurrent(offsurf);

    int FullXdim = ceil(fFBOZoom * InitXImage);
    int FullYdim = ceil(fFBOZoom * InitYImage);
    image_fbo=new QOpenGLFramebufferObject(FullXdim, FullYdim, QOpenGLFramebufferObject::CombinedDepthStencil,GL_TEXTURE_2D,GL_RGBA);
    image_fbo->bind();
    PaintDevice = new QOpenGLPaintDevice;
    PaintDevice->setSize(QSize(FullXdim,FullYdim));
    orthoproj.setToIdentity();
    orthoproj.ortho(SD->Limits);
    translation.setToIdentity();
    scaling.setToIdentity();
    scaling.scale(fSizeFactor,-fSizeFactor,1.0);
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
    glViewport(0,0,FullXdim,FullYdim);
    paintGL();
    CreateTIFF(FullXdim, FullYdim, FrameBuffer);
    ctx->doneCurrent();
    image_fbo->release();
    delete image_fbo;
    delete PaintDevice;
    delete offsurf;

// restore settings for display
    makeCurrent();
    PaintDevice = new QOpenGLPaintDevice;
    PaintDevice->setSize(QSize(XImageSize,YImageSize));
    orthoproj.setToIdentity();
    orthoproj.ortho(DisplayLimits);
    translation.setToIdentity();
    translation.setColumn(3,QVector4D(ImageOffset.x(),ImageOffset.y(),0.0,1.0));
    scaling.setToIdentity();
    scaling.scale(fZoom*fSizeFactor,-fZoom*fSizeFactor,1.0);
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
    update();
}
void ViewClusters::setupTransforms()
{
    translation.setToIdentity();
    orthoproj.setToIdentity();
    orthoproj.ortho(SD->Limits);
    DisplayLimits = SD->Limits;
    Centre = Point_2(SD->Limits.width()/2.0, SD->Limits.height()/2.0);
    scaling.setToIdentity();
    scaling.scale(fZoom*fSizeFactor,-fZoom*fSizeFactor,1.0);
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
}
void ViewClusters::setupClusterSearch()
{
    ctree.clear();
    std::list<APoint>::iterator l;
    for(l=SD->CentroidList.begin(); l!=SD->CentroidList.end(); ++l)
       ctree.insert(*l);
}
void ViewClusters::reassignClusterNos()
{
    int nNotetramers = int(SD->tm.size());
    if (nNotetramers == 0)
        return;
    emit AlertMsg("Reassigning placed tetramer cluster numbers",'g');
    for(int k=0; k < nNotetramers; k++)
    {
        Tetramer itm = SD->tm[k];
        int cno= findCluster(Point_2(itm.xcentre, itm.ycentre));
        if(cno == -1)
        {
            emit AlertMsg(QString("Cannot reassign cluster number for tetramer at %1,%2").arg(itm.xcentre).arg(itm.ycentre),'r');
            continue;
        }
        itm.ClusterNo = cno;
        SD->tm[k] = itm;
    }
}
uint ViewClusters::ClusterCount()
{
    TetMappedCluster.clear();
    Tetcount.clear();
    nNoTetramers = uint(SD->tm.size());
    if(nNoTetramers == 0)
      return uint(0);

    for(uint i=0; i < nNoTetramers; i++)
       TetMappedCluster.push_back(SD->tm[i].ClusterNo);

    uint j = 1;
    uint numTet = 1;
    std::sort(TetMappedCluster.begin(),TetMappedCluster.end());
    while (j < nNoTetramers)
    {

        if(TetMappedCluster[j-1] != TetMappedCluster[j])
        {
            std::pair<uint,uint> tc = std::make_pair(TetMappedCluster[j-1], numTet);
            Tetcount.insert(tc);
            numTet = 0;
        }
        j++;
        numTet++;
    }
    if(numTet > 0)
    {
        std::pair<uint,uint> tc = std::make_pair(TetMappedCluster[j-1], numTet);
        Tetcount.insert(tc);
    }
    TetMappedCluster.erase(std::unique(TetMappedCluster.begin(),TetMappedCluster.end()), TetMappedCluster.end());
    return(uint(TetMappedCluster.size()));
}
void ViewClusters::displayClusterInfo(uint cno)
{
    Cluster c = SD->Clusters[cno];
    std::tuple<double,double> el = c.getEllipseData();
//    double Eccentricity = std::get<0>(el)/std::get<1>(el);
    double cArea = c.getArea();
    double BlinkDensity = c.getNoBlinks()/cArea;
    double PointDensity = c.getClustSize()/cArea;
    QString cinfo[9];
    int nString = 6;
    int nSubClustOffset =0;
    int NoSubClusters=c.getNoSubClusters();
    if(NoSubClusters > 1)
    {
       nString++;
       nSubClustOffset=1;
    }
    cinfo[0] = QString("Cluster No. %1 - Neighbourhood Limit = %2 nm").arg(c.getClustNo()).arg(SD->NeighbourhoodLimit);
    cinfo[1] = QString("Area = %1 nm<sup>2</sup> -- Alpha = %2 -- Convex Hull Area = %3 nm<sup>2</sup>").arg(cArea).arg(SD->Alpha1stPass).arg(c.getConvexHullArea());
    cinfo[2] = QString("# Blinks = %1 - Density = %2").arg(c.getNoBlinks()).arg(BlinkDensity,0,'f',3);
    cinfo[3] = QString("# Points = %1 - Density = %2").arg(c.getClustSize()).arg(PointDensity,0,'f',3);
    if(nSubClustOffset > 0)
    {
       std::vector<std::pair<uint, double>> sc = c.getBaseClusterInfo();
       QString SubClust = QString("No. of SubClusters = %1 - Areas: ").arg(NoSubClusters);
       QString punc=", ";
       for(int k = 0; k < NoSubClusters; k++)
       {
          if(k == NoSubClusters-1)
              punc = " nm<sup>2</sup>";
          SubClust += QString("%1 %2").arg(sc[k].second).arg(punc);
       }
       cinfo[4]=SubClust;
    }
    cinfo[4+nSubClustOffset] = QString("NND = %1 nm - NearestArea = %2 nm<sup>2</sup>").arg(c.getNND(),0,'f',1).arg(c.getNearestArea(),0,'f',1);
    cinfo[5+nSubClustOffset] = QString("Ellipse : Min = %1 nm - Max = %2 nm").arg(std::get<1>(el),0,'f',1).arg(std::get<0>(el),0,'f',1);
    if(nLabelledClusters  > 0)
    {
        auto search = Tetcount.find(cno);
        if (search != Tetcount.end()) {
          nString += 2;
          int nTetinCluster = search->second;
          double Occupancy = nTetinCluster * tetramerwidth*tetramerwidth /cArea;
          double OccupancyCom = nTetinCluster * 900.0/cArea;
          cinfo[6+nSubClustOffset] = QString("Tetramers = %1 - Occupancy = %2 (%3 nm) = %4 (30 nm)").arg(nTetinCluster).arg(Occupancy,0,'f',3).arg(tetramerwidth,0,'f',0).arg(OccupancyCom,0,'f',3);
          cinfo[7+nSubClustOffset] = QString("Blinks/Tetramer = %1").arg(c.getNoBlinks()/nTetinCluster);
        }
    }
    emit AlertMsg("---------------------------------",'b');
    for(int k=0; k < nString; k++)
        emit AlertMsg(cinfo[k],'b');
    emit AlertMsg("---------------------------------",'b');
}
void ViewClusters::imageLine()
{
    bLinePlot = !bLinePlot;
    if(bLinePlot)
       StartPt = EndPt = Point_2(0.0f,0.0f);
    update();
}
void ViewClusters::increasePointSize()
{
    PointSize += 1;
    emit AlertMsg(QString("PointSize = %1").arg(PointSize),'b');
    update();
}
void ViewClusters::increaseExPointSize()
{
    ExPointSize += 1;
    emit AlertMsg(QString("Excluded Points Size = %1").arg(ExPointSize),'b');
    update();
}
void ViewClusters::decreasePointSize()
{
    if(PointSize == 1)
      return;
    PointSize -= 1;
    emit AlertMsg(QString("PointSize = %1").arg(PointSize),'b');
    update();
}
void ViewClusters::decreaseExPointSize()
{
    if(ExPointSize == 1)
      return;
    ExPointSize -= 1;
    emit AlertMsg(QString("Excluded Points Size = %1").arg(ExPointSize),'b');
    update();
}
void ViewClusters::keyPressEvent( QKeyEvent *e )
{
    bool bAlterTetramer = bModifyTetramer || bAddTetramer;

    if(bAlterTetramer)
    {
        validAddModifyKeys(e);
        return;
    }
    if(bClassifyTetramer)
    {
        validClassifyKeys(e);
        return;
    }
    if(bDisplayClusterInfo)
    {
        if(e->key() == Qt::Key_F4)
          clusterInfo();
        else
          emit AlertMsg("Invalid key pressed in Cluster Info mode",'r');
        return;
    }
    switch( e->key() )
    {
    case Qt::Key_Escape:
        if(bROI)
        {
          toggleROI();
          return;
        }
        CloseWindow();
        break;
    case Qt::Key_A:
        addTetramer();
        break;
    case Qt::Key_C:
        classifyTetramer();
        break;
    case Qt::Key_L: //LinePlot
        imageLine();
        break;
    case Qt::Key_M:
        modifyTetramer();
        break;
    case Qt::Key_N:
        showNearestNeighbours();
        break;
    case Qt::Key_P:
        if(nPalette == 1)
           nPalette = 2;
        else
           nPalette = 1;
        update();
        break;
    case Qt::Key_F1:
        emit showHelp();
        break;
    case Qt::Key_F2:
        resetDisplay();
        break;
    case Qt::Key_F3:
        toggleROI();
        break;
    case Qt::Key_F4:
        clusterInfo();
        break;
    case Qt::Key_F5:
        getScaleBarSize();
        break;
    case Qt::Key_F6:
        isolateTetramers();
        break;
    case Qt::Key_F7:
        displayTetramers();
        break;
    case Qt::Key_F8:
        changeTetramerWidth();
        break;
    case Qt::Key_F9:
        saveTetramerFile();
        break;
    case Qt::Key_F10:      //Increase Image Pointsize
        if (e->modifiers() & Qt::AltModifier)
            increaseExPointSize();
        else
            increasePointSize();
        break;
    case Qt::Key_F11:     //Decrease Image Pointsize
        if (e->modifiers() & Qt::AltModifier)
            decreaseExPointSize();
        else
            decreasePointSize();
        break;
    case Qt::Key_F12:    //Create Tiff Image
        if (e->modifiers() & Qt::AltModifier)
          CreateFBOImage();
        else
          CreateScreenShot();
        break;
    case Qt::Key_Plus:  // Zoom In
        fZoomDelta = 0.25;
        if (e->modifiers() & Qt::AltModifier)
            fZoomDelta = 0.5;
        if (e->modifiers() & Qt::ControlModifier)
            fZoomDelta = 1.0;
        PosnWanted=Point_2(Centre.x()-ImageOffset.x(),Centre.y()-ImageOffset.y());
        imageZoomIn();
        break;
    case Qt::Key_Minus: // Zoom Out
        fZoomDelta = 0.25;
        if (e->modifiers() & Qt::AltModifier)
            fZoomDelta = 0.5;
        if (e->modifiers() & Qt::ControlModifier)
            fZoomDelta = 1.0;
        PosnWanted=Point_2(Centre.x()-ImageOffset.x(),Centre.y()-ImageOffset.y());
        imageZoomOut();
        break;
    default:
//        emit AlertMsg(QString("Key %1 - not assigned a function - ignored").arg(e->key()),'r');
        break;
    }
}
void ViewClusters::validAddModifyKeys(QKeyEvent *e)
{
    switch( e->key() )
    {
    case Qt::Key_Escape:
        setCursor(Qt::ArrowCursor);
        if(bAddTetramer)
        {
            bAddTetramer = false;
            if(bPlacingTetramer)
            {
                emit AlertMsg(QString("Tetramer not added"),'r');
                bPlacingTetramer = false;
                loadTetramerData();
                update();
            }
            return;
        }
        if(bModifyTetramer)
        {
            bModifyTetramer = false;
            if(nBeingModified > -1)
            {
                emit AlertMsg(QString("Tetramer %1 NOT modified").arg(nBeingModified + 1),'r');
                nBeingModified = -1;
                loadTetramerData();
                update();
            }
            return;
        }
        break;
    case Qt::Key_Right:
        ntm.xcentre += PixelSize;
        loadSingleTetramerData();
        break;
    case Qt::Key_Left:
        ntm.xcentre -= PixelSize;
        loadSingleTetramerData();
        break;
    case Qt::Key_Up:
        ntm.ycentre += PixelSize;
        loadSingleTetramerData();
        break;
    case Qt::Key_Down:
        ntm.ycentre -= PixelSize;
        loadSingleTetramerData();
        break;
    case Qt::Key_D:
        deleteTetramer();
        break;
    case Qt::Key_S:
        saveTetramer();
        break;
    default:
        emit AlertMsg(QString("Invalid key %1 during add/modify").arg(e->text()),'r');
        break;
    }
}
void ViewClusters::validClassifyKeys(QKeyEvent *e)
{
    switch( e->key() )
    {
    case Qt::Key_Escape:
        setCursor(Qt::ArrowCursor);
        bClassifyTetramer = false;
        if(nBeingClassified > -1)
        {
           emit AlertMsg(QString("Tetramer %1 NOT classified").arg(nBeingClassified + 1),'r');
           nBeingClassified = -1;
           loadTetramerData();
           update();
        }
        return;
        break;
    case Qt::Key_B:
        tetramerClassification('b');
        break;
    case Qt::Key_C:
        tetramerClassification('c');
        break;
    case Qt::Key_E:
        bClassifyTetramer = false;
        nBeingClassified = -1;
        loadTetramerData();
        update();
        break;
    case Qt::Key_I:
        tetramerClassification('i');
        break;
    case Qt::Key_S:
        tetramerClassification('s');
        break;
    case Qt::Key_U:
        tetramerClassification('u');
        break;
    default:
        emit AlertMsg(QString("Invalid key %1 for classification").arg(e->text()),'r');
        break;
    }
}
void ViewClusters::CloseWindow()
{
    if(bPlacingTetramer || bModifyTetramer || bClassifyTetramer)
    {
        emit AlertMsg("Cannot close program while adding or modifying tetramers",'r');
    }
    if(nNoTetramers > 0 && (!bSavedTetFile || !bCalcNearestNeighbours))
      saveTetramerFile();
    emit CloseProgram();
}
void ViewClusters::toggleROI()
{
    bROI = !bROI;
    if(bROI)
    {
      if(!roiinfo)
        roiinfo = new BlinkROI(this);
      setCursor(Qt::CrossCursor);
      ROIAct->setText("Hide ROI {F3}");
    }
    else
    {
      if(rubberBand)
         roiinfo->OnCancel();
       setCursor(Qt::ArrowCursor);
       ROIAct->setText("Create ROI {F3}");
    }
}
void ViewClusters::hideROI()
{
    bROI = false;
    rubberBand->hide();
    rubberBand->setGeometry(0,0,0,0);
    ROIAct->setText("Create ROI {F3}");
}
void ViewClusters::changeTetramerWidth()
{
    bool ok;
    double tsize = QInputDialog::getDouble(nullptr, tr("Change Tetramer width"),
                                         tr("Tetramer width (nm):"), 27.0, 24.0, 50.0, 1, &ok);
    if(!ok)
      return;
    tetramerwidth = tsize;
    tetramerdiag = tetramerwidth/sqrt2;
    if(bPlacingTetramer||bModifyTetramer||bClassifyTetramer)
        loadSingleTetramerData();
    loadTetramerData();
    update();
}
void ViewClusters::clusterInfo()
{
    bDisplayClusterInfo = !bDisplayClusterInfo;
    if(bDisplayClusterInfo)
    {
       setCursor(Qt::PointingHandCursor);
       nLabelledClusters = ClusterCount();
    }
    else
       setCursor(Qt::ArrowCursor);
    update();
}
void ViewClusters::addTetramer()
{
    setCursor(Qt::CrossCursor);
    bPlotTetramers=true;
    bAddTetramer = true;
    bPlacingTetramer=false;
    bDisplayNND = false;
    emit AlertMsg(" Adding a tetramer",'b');
}
void ViewClusters::modifyTetramer()
{
    if(nNoTetramers == 0)
    {
        emit AlertMsg(" No tetramers to modify",'r');
        return;
    }
    setCursor(Qt::CrossCursor);
    bModifyTetramer = true;
    bDisplayNND = false;
    emit AlertMsg(" Modifying a tetramer",'b');
}
void ViewClusters::classifyTetramer()
{
    if(nNoTetramers == 0)
    {
        emit AlertMsg(" No tetramers to classify",'r');
        return;
    }
    setCursor(Qt::CrossCursor);
    bClassifyTetramer = true;
    bDisplayNND = false;
    emit AlertMsg(" Classifying a tetramer",'b');
}
void ViewClusters::deleteTetramer()
{
    if(bAddTetramer)
    {
        bAddTetramer = false;
        bPlacingTetramer = false;
    }
    if(bModifyTetramer)
    {
        if(nBeingModified > -1)
        {
            SD->tm.erase(SD->tm.begin()+nBeingModified);
            bModifyTetramer = false;
            loadTetramerData();
        }
    }
    bCalcNearestNeighbours = true;
    update();
}
void ViewClusters::saveTetramer()
{
    // save a tetramer

    if(bAddTetramer)
    {
        ntm.classification = 'u';
        bAddTetramer = false;
        bPlacingTetramer = false;
        SD->tm.push_back(ntm);
    }
    if(bModifyTetramer)
    {
        SD->tm[nBeingModified] = ntm;
        bModifyTetramer = false;
        nBeingModified = -1;
    }
    loadTetramerData();
    bCalcNearestNeighbours = true;
    update();
}
void ViewClusters::tetramerClassification(char c)
{
/*************** classification *************
     c = checkerboard
     s = side by side
     b = both checkerboard and side by side
     i = isolated
     u = unclassified
*********************************************/
    char ttype[] = {'c', 's', 'b', 'i', 'u'};
    int k = -1;
    for(int i = 0; i < 5; i++)
    {
        if(c == ttype[i])
        {
            k = i;
            break;
        }
    }
    if(k == -1)
        c = 'u';

    emit AlertMsg(QString("Tetramer %1 classied as %2").arg(nBeingClassified).arg(c),'b');
    SD->tm[uint(nBeingClassified)].classification = c;
    nBeingClassified = -1;
    bClassifyTetramer = false;
    loadTetramerData();
    update();
}
int ViewClusters::findCluster(Point_2 posn)
{
    // first find nearest centroids to point

    std::vector<uint> NearClusters;
//    std::cerr << "The point (" << posn.x() << "," << posn.y() << ")\n";
    APoint clusterquery = APoint(posn.x(),posn.y());
    Neighbor_search2 searchc(ctree,clusterquery, 4);
    Neighbor_search2::iterator ittc = searchc.begin();
    int maxx = int(SD->Limits.right());
    for(int iz=0; iz < 4; iz++)
    {
        APoint NearC = (*ittc).first;
        ittc++;
        int cOffset = int(NearC.y()) * maxx + int(NearC.x());
        int cno = SD->CentroidMap.value(cOffset, -1);
        if(cno == -1)
        {
            emit AlertMsg("findCluster: Cannot find cluster in map",'r');
            return -1;
        }
//        emit AlertMsg(QString("Found cluster %1").arg(cno),'b');
        NearClusters.push_back(uint(cno));
    }
    int matchedCluster = -1;
    for(uint ic=0; ic < NearClusters.size(); ic++)
    {
        uint icn=NearClusters[ic];
        Polygon_2 pgn = SD->Clusters[icn].getConvexHull();

        CGAL::Bounded_side bs = pgn.bounded_side(clusterquery);

        if(bs == CGAL::ON_BOUNDED_SIDE)
        {
             matchedCluster = int(icn);
//             emit AlertMsg(QString("findCluster: Matched point to cluster %1").arg(icn),'b');
             break;
        }
    }
    if(matchedCluster == -1)
        emit AlertMsg(QString("findCluster: Could not match point to a cluster"),'r');
    return matchedCluster;
}
bool ViewClusters::findTetramer(Point_2 posn)
{
    int clusterNo = findCluster(posn);
    if (clusterNo == -1)
    {
        bClassifyTetramer = bModifyTetramer = false;
        return false;
    }
    else
    {
        emit AlertMsg(QString("Identified cluster %1").arg(clusterNo),'b');
        ntm.ClusterNo = clusterNo;
    }
    double tmeas= tetramerwidth*tetramerwidth/4.0;

    if(bAddTetramer)
    {
        double xdiff = (ntm.xcentre - posn.x());
        double ydiff = (ntm.ycentre - posn.y());
        double dist = xdiff*xdiff + ydiff*ydiff;
        if (dist < tmeas)
        {
            loadSingleTetramerData();
        }
    }
    else
    {
        for (uint k = 0; k < nNoTetramers; k++)
        {
            if(SD->tm[k].ClusterNo != clusterNo)
                continue;

            double xdiff = (SD->tm[k].xcentre - posn.x());
            double ydiff = (SD->tm[k].ycentre - posn.y());
            double dist = xdiff*xdiff + ydiff*ydiff;
            if (dist < tmeas)
            {
                if(bModifyTetramer)
                    nBeingModified = int(k);

                if(bClassifyTetramer)
                    nBeingClassified = int(k);

                ntm=SD->tm[k];
                loadSingleTetramerData();
                loadTetramerData();
                update();
                break;
            }
        }
        if(nBeingModified == -1 && nBeingClassified == -1)
        {
            bClassifyTetramer = bModifyTetramer = false;
            emit AlertMsg("Cannot find tetramer to be modified or classified",'r');
            return false;
        }
    }
    return true;
}
void ViewClusters::mouseDoubleClickEvent( QMouseEvent *e )
{
    bDoubleClick = true;
    Point_2 whereClicked = LimitMouse(e->x(),e->y());
    PosnWanted = to2DObjectPoint(whereClicked);
    fZoomDelta = 0.25f;
    imageZoomIn();
}
void ViewClusters::mousePressEvent( QMouseEvent *e )
{
    if(bPlacingTetramer || bDoubleClick)
        return;
    Point_2 mousePosn= LimitMouse(e->x(),e->y());
    Point_2 ImagePt = to2DObjectPoint(mousePosn);
    StartImagePt = ImagePt;
//    emit AlertMsg(QString("Mouse Posn at Start = (%1, %2)").arg(mousePosn.x()).arg(mousePosn.y()),'r');
    switch (e->button())
    {
     case Qt::LeftButton:
        StartPt = EndPt = mousePosn;
        bLeftButtonDown = true;
        if(bROI)
        {

            winoffset = mapToGlobal(QPoint(0,0));
            int xoffset = mousePosn.x()+winoffset.x();
            int yoffset = mousePosn.y()+winoffset.y();
//            QPoint convert = mapToGlobal(QPoint(e->x(),e->y()));
//            emit AlertMsg(QString("Posn using winoffset = %1 %2; ToGlobal = %3 %4").arg(xoffset).arg(yoffset).arg(convert.x()).arg(convert.y()),'r');
            if (!rubberBand)
                rubberBand = new QRubberBand(QRubberBand::Rectangle);
            rubberBand->setGeometry(QRect(xoffset, yoffset, 0, 0));
            rubberBand->show();
            return;
        }
        if(bModifyTetramer || bClassifyTetramer)
        {
            setCursor(Qt::ArrowCursor);
            if(!findTetramer(ImagePt))
                return;
        }
        if(bDisplayClusterInfo)
        {
            int ClusterNo = findCluster(ImagePt);
            if(ClusterNo == -1)
                return;

            CurrentCluster = ClusterNo;
            displayClusterInfo(ClusterNo);
            return;
        }
        if(bAddTetramer)
        {
            setCursor(Qt::ArrowCursor);
            ntm.xcentre = ImagePt.x();
            ntm.ycentre = ImagePt.y();
            ntm.ClusterNo = findCluster(ImagePt);
            if(ntm.ClusterNo == -1)
            {
                emit AlertMsg("Adding tetramer - Could not identify cluster", 'r');
                return;
            }
            ntm.angle = 0.0;
//           emit AlertMsg(QString("Adding tetramer at %1,%2").arg(ntm.xcentre).arg(ntm.ycentre), 'b');
            bPlacingTetramer = true;
            loadSingleTetramerData();
            update();
            return;
        }
        if(bLinePlot)
            setCursor(Qt::CrossCursor);
        else { if (bDisplayClusterInfo)
            setCursor(Qt::PointingHandCursor);
        else
            setCursor(Qt::OpenHandCursor);}

        break;
     case Qt::MiddleButton:
        if(bScaleBar && !bScaleBarPositionChosen)
        {

          ScaleBarPosition = to2DObjectPoint(mousePosn);
//         emit AlertMsg(QString("Mouse Posn = %1 %2 - Image posn = %3 %4").arg(mousePosn.x()).arg(mousePosn.y()).arg(ScaleBarPosition.x()).arg(ScaleBarPosition.y()),'r');
          bScaleBarPositionChosen = true;
          CalculateScaleBar();
        }
        break;
    case Qt::RightButton:
        imageMenu->exec(e->globalPos());
        break;
    default:
        QWindow::mousePressEvent(e);
        break;
    }
}
void ViewClusters::mouseReleaseEvent( QMouseEvent *e )
{
    if(bDoubleClick)
       bDoubleClick = false;

    EndPt = LimitMouse(e->x(), e->y());
    EndImagePt = to2DObjectPoint(EndPt);
//    emit AlertMsg(QString("Mouse Posn on Release = (%1, %2)").arg(EndImagePt.x()).arg(EndImagePt.y()),'r');

    switch (e->button())
    {
     case Qt::LeftButton:
        bLeftButtonDown=false;
        if(bROI)
        {
           QRectF ROIimage = rubberBand->geometry();
           QPointF topLeftnm = to2DObjectPoint(ROIimage.topLeft()-winoffset);
           QPointF bottomRightnm = to2DObjectPoint(ROIimage.bottomRight()-winoffset);
           QRectF ROInm = QRectF(topLeftnm,bottomRightnm).normalized();
           setCursor(Qt::ArrowCursor);
           roiinfo->displayMenu(ROInm, SD->FileName);
           return;
        }
        if(bLinePlot)
        {
            if(StartPt != EndPt)
                update();
        }
        break;
    default:
        QWindow::mouseReleaseEvent(e);
        break;

     }
     if(!bDisplayClusterInfo)
        unsetCursor();
}
void ViewClusters::mouseMoveEvent( QMouseEvent *e )
{
    Point_2 newMousePosn = LimitMouse(e->x(),e->y());
    Point_2 newImagePosn = to2DObjectPoint(newMousePosn);
    if(bDoubleClick)
        return;
    if(bLeftButtonDown)
    {
        if(bROI)
        {
           EndPt = newMousePosn;
           rubberBand->setGeometry(QRect(QPoint(StartPt.x()+winoffset.x(), StartPt.y()+winoffset.y()), QPoint(EndPt.x()+winoffset.x(),EndPt.y()+winoffset.y())).normalized());
           return;
        }
        if(bLinePlot)
        {
            EndPt = newMousePosn;
//            emit AlertMsg(QString("New Mouse Posn = (%1, %2)").arg(newMousePosn.x()).arg(newMousePosn.y()),'r');
            update();
        }
        else
        {
            Point_2 newOffset = Point_2(newImagePosn.x() - StartImagePt.x(), newImagePosn.y() - StartImagePt.y() );
            ImageOffset = Point_2(ImageOffset.x() + newOffset.x(), ImageOffset.y() + newOffset.y());
            changeTranslation(ImageOffset.x(), ImageOffset.y());
            update();
        }
    }
}
void ViewClusters::wheelEvent(QWheelEvent *event)
{
    int angle = event->angleDelta().y();
    if(bPlacingTetramer || bModifyTetramer)
    {
        ntm.angle += int(angle/40);
        ntm.angle %= 90;
        loadSingleTetramerData();
        return;
    }
/*    PosnWanted=ImageOffset;
    fZoomDelta = +0.1f;
    if(angle > 0)
       imageZoomIn();
    else
       imageZoomOut(); */
}
Point_2 ViewClusters::LimitMouse(int x, int y)
{
    // Transform y coordinate
    if(x < 0) x = 0;
    else if(x > XImageSize-1)x = XImageSize-1;

    if(y < 0)y = 0;
    else if(y > YImageSize-1)y = YImageSize-1;

    return Point_2(x, y);
}
void ViewClusters::imageZoomIn()
{
    fOldZoom = fZoom;
    if(fZoom >= 1.0f)
       fZoom += fZoomDelta;
    else
       fZoom += 0.1f;

    if (fZoom <= 0.0f)
    {
      emit AlertMsg(QString("Zoom is zero - reset to 1.0"),'r');
      fZoom = 1.0f;
    }
    changeSize();
}
void ViewClusters::imageZoomOut()
{
    fOldZoom = fZoom;
    if(fZoom >= 1.0f)
       fZoom -= fZoomDelta;
    else
       fZoom -= 0.25f;


    if (fZoom <= 0.0f)
    {
      emit AlertMsg(QString("Zoom is zero - reset to 1.0"),'r');
      fZoom = 1.0f;
    }
    changeSize();
}
void ViewClusters::changeSize()
{
    if(std::abs(fZoom - 1.00f) < FLT_EPSILON)
    {
        XImageSize = InitXImage;
        YImageSize = InitYImage;
        fZoom = 1.0f;
        changeProjection(1.0f,1.0f);
        resize(XImageSize, YImageSize);
    }
    else //if(fZoom != 1.0)
    {
        float XnewImageSize = qMin(InitXImage*fZoom, fXScreenMax);
        float YnewImageSize = qMin(InitYImage*fZoom, fYScreenMax);
        changeProjection(XnewImageSize/InitXImage, YnewImageSize/InitYImage);
        if (std::abs(XnewImageSize - XImageSize) > 0 || std::abs(YnewImageSize - YImageSize) > 0 )
        {
            XImageSize = XnewImageSize;
            YImageSize = YnewImageSize;
            resize(XImageSize, YImageSize);
        }
    }
    changeScalingandTranslation(ImageOffset.x(), ImageOffset.y());
    calculatePixelSize();
    emit AlertMsg(QString("Zoom = %1 - PixelSize = %2 nm").arg(fZoom).arg(PixelSize,0,'f',1),'m');
    if(bLinePlot)
       bLinePlot = !bLinePlot;

    if(bScaleBar)
       CalculateScaleBar();
    update();
}
void ViewClusters::changeProjection(float xm, float ym)
{
    double x, y, width, height;
    SD->Limits.getRect(&x, &y, &width, &height);
    orthoproj.setToIdentity();
    DisplayLimits = QRectF(float(x), float(y), float(xm*width), float(ym*height));
    orthoproj.ortho(DisplayLimits);
    float NewXCentre = DisplayLimits.width()/2.0;
    float NewYCentre = DisplayLimits.height()/2.0;
//    emit AlertMsg(QString("Old Centre = %1, %2").arg(Centre.x()).arg(Centre.y()),'r');
//    emit AlertMsg(QString("New Centre = %1, %2").arg(NewXCentre).arg(NewYCentre),'r');
    float NewXoffset = (NewXCentre - PosnWanted.x());
    float NewYoffset = (NewYCentre - PosnWanted.y());
    Centre = Point_2(NewXCentre, NewYCentre);
    ImageOffset = Point_2(NewXoffset, NewYoffset);
}
void ViewClusters::changeTranslation(float xs, float ys)
{
//    emit AlertMsg(QString("Moving image by %1, %2").arg(xs).arg(ys),'r');
    translation.setToIdentity();
    translation.setColumn(3,QVector4D(xs,ys,0.0,1.0));
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
}
void ViewClusters::changeScalingandTranslation(float xs, float ys)
{
    scaling.setToIdentity();
    scaling.scale(fZoom*fSizeFactor,-fZoom*fSizeFactor,1.0);
    translation.setToIdentity();
    translation.setColumn(3,QVector4D(xs,ys,0.0,1.0));
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
}
void ViewClusters::changeScaling()
{
    scaling.setToIdentity();
    scaling.scale(fZoom*fSizeFactor,-fZoom*fSizeFactor,1.0);
    transform=scaling*orthoproj*translation;
    invtransform = transform.inverted();
}
void ViewClusters::resetDisplay()
{
    fZoom=1.0f;
    xshift = 0.0f;
    yshift = 0.0f;
    ImageOffset = Point_2(0,0);
    Centre = Point_2(SD->Limits.width()/2.0, SD->Limits.height()/2.0);
    resize(InitXImage, InitYImage);
    XImageSize = InitXImage;
    YImageSize = InitYImage;
    PaintDevice->setSize(QSize(InitXImage,InitYImage));
    setupTransforms();
    calculatePixelSize();
    emit AlertMsg(QString("Reset - Zoom = %1 - PixelSize = %2 nm").arg(fZoom).arg(PixelSize,0,'f',1),'m');
    update();
}
void ViewClusters::calculatePixelSize()
{
    double XPixelSize =  SD->Limits.width()/(fZoom*InitXImage);
    double YPixelSize = SD->Limits.height()/(fZoom*InitYImage);
    if(abs(XPixelSize - YPixelSize) > 1.e-2)
    {
         emit AlertMsg(QString("X & Y Pixel sizes differ: X = %1 nm; Y = %2 nm ").arg(XPixelSize,1).arg(YPixelSize,1),'m');
         PixelSize = (XPixelSize + YPixelSize)/2.0;
    }
    else
    {
         PixelSize = XPixelSize;
    }
}
void ViewClusters::getScaleBarSize()
{
    bScaleBar=!bScaleBar;
    bScaleBarPositionChosen = false;
    if(bScaleBar)
    {
      int minsize=20;
      nScaleBarLength=QInputDialog::getInt(nullptr,"ScaleBar Value","Enter Scalebar length (nm)", nScaleBarLength, minsize, 4000, 10);
    }
    else
      update();
}
void ViewClusters::calculateNearestNeighbours()
{
    NNDtet.clear();
    nNoTetramers = uint(SD->tm.size());
    if(nNoTetramers < 2)
       return;

    uint nNoClusters = ClusterCount();
    emit AlertMsg(QString("No of populated clusters = %1").arg(nNoClusters),'b');
    if(nNoClusters == 1)
       ClusterNearestNeighbours(SD->tm);
    else
    {
       uint jstart=0;
       uint jend = 0;
       std::sort(SD->tm.begin(), SD->tm.end());
       for(uint k=0; k < nNoClusters; k++)
       {
         int ClusterNo = TetMappedCluster[k];
         uint nNoTetramersinCluster = Tetcount[ClusterNo];
         jend += nNoTetramersinCluster;
         std::vector<Tetramer> tmj;

         for(uint j=jstart; j < jend; j++)
         {
            if(SD->tm[j].ClusterNo != ClusterNo)
            {
               emit AlertMsg(QString("Mapping Error: Tetramer %1 Cluster %2 != %3 (tm value)<br>Unable to calculate nearest neighbours").arg(j).arg(ClusterNo).arg(SD->tm[k].ClusterNo),'r');
               bDisplayNND = false;
               return;
            }
            tmj.push_back(SD->tm[j]);
         }
         jstart = jend;
         if(tmj.size() > 1)
           ClusterNearestNeighbours(tmj);
         tmj.clear();
       }
    }

    bDisplayNND = true;
    bCalcNearestNeighbours = false;
    update();
}
void ViewClusters::ClusterNearestNeighbours(std::vector<Tetramer> tetm)
{
    uint nTetramers = uint(tetm.size());
    // Because there are usually less than 25 tetramers in a dyad
    // do a simple search rather than setting up a tree
    for(uint j = 0; j < nTetramers; j++)
    {
        double xpos = tetm[j].xcentre;
        double ypos = tetm[j].ycentre;

        double nndmin = 10000000;
        NNDval a;
        a.ClustNo = tetm[j].ClusterNo;
        a.NNDStart = Point_2(xpos,ypos);

        for(uint i = 0; i < nTetramers; i++)
        {
            if(i == j) continue;

            double xd = xpos-tetm[i].xcentre;
            double yd = ypos-tetm[i].ycentre;
            double dist=xd*xd + yd*yd;
            nndmin = qMin(nndmin,dist);
            if(std::abs(dist - nndmin) < DBL_EPSILON)
            {
              a.NND = sqrt(nndmin);
              a.NNDEnd = Point_2(tetm[i].xcentre, tetm[i].ycentre);
            }
        }
        NNDtet.push_back(a);
    }
}
void ViewClusters::displayTetramers()
{
     bPlotTetramers = !bPlotTetramers;
     if(bPlotTetramers)
     {
         displayTetramersAct->setText("Hide Tetramers {F7}");
         loadTetramerData();
     }
     else
     {
         displayTetramersAct->setText("Show Tetramers {F7}");
         if(bDisplayNND)
         {
             showNearestNeighboursAct->setText("Show Nearest Neighbours {n}");
             bDisplayNND =false;
         }
     }
     update();
}
void ViewClusters::showNearestNeighbours()
{
     if(nNoTetramers == 0)
        return;
     bDisplayNND = !bDisplayNND;
     if(bDisplayNND)
     {
         bShowTetramer = true;
         if(bCalcNearestNeighbours)
            calculateNearestNeighbours();
         showNearestNeighboursAct->setText("Hide Nearest Neighbours {n}");
     }
     else
         showNearestNeighboursAct->setText("Show Nearest Neighbours {n}");

     update();
}
void ViewClusters::isolateTetramers()
{
    if(nNoTetramers == 0)
       return;
    bIsolateTetramers = !bIsolateTetramers;
    if(bIsolateTetramers)
    {
        nPalette = 1;
        IsolateTetramersAct->setText("Unisolate Tetramers (Show Blinks) {F6}");
    }
    else
    {
        nPalette = 2;
        IsolateTetramersAct->setText("Isolate Tetramers (Hide Blinks) {F6}");
    }
    loadTetramerData();
    update();
}
void ViewClusters::initializeGL()
{
    QSurfaceFormat format;
    format.setDepthBufferSize( 4 );
    format.setSamples(24);
    format.setVersion(4,3);
    format.setRenderableType(QSurfaceFormat::OpenGL);
    setFormat(format);
    PaintDevice = new QOpenGLPaintDevice;
    PaintDevice->setSize(QSize(InitXImage,InitYImage));
    makeCurrent();

    initializeOpenGLFunctions();
    const char* Renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));
    emit AlertMsg(QString("Renderer : %1").arg(Renderer),' ');
    std::string Vendor = reinterpret_cast<const char *>(glGetString(GL_VENDOR));
    emit AlertMsg(QString("Vendor : %1").arg(Vendor.c_str()),' ');
    if(Vendor.find("NVIDIA") != std::string::npos)
    {
        int GPU_Memory;
        glGetIntegerv(GL_GPU_MEMORY_INFO_TOTAL_AVAILABLE_MEMORY_NVX, &GPU_Memory);
        QString unit = " kB";
        if (GPU_Memory > 100000)
        {
           GPU_Memory /= 1024;
           unit = " MB";
        }
        emit AlertMsg(QString("GPU memory available is %1 %2").arg(GPU_Memory).arg(unit),' ');
    }
    const char* Version = reinterpret_cast<const char *>(glGetString(GL_VERSION));
    emit AlertMsg(QString("Open GL version = %1").arg(Version),' ');

    glEnable(GL_PROGRAM_POINT_SIZE);
    initShaders();
}
template <class T>
T ViewClusters::to2DViewportPoint(T posn)
{
    QVector4D rp = QVector4D(posn.x(),posn.y(),0.0f,1.0f);
    QVector4D rt = transform*rp;
    int newx = int(XImageSize*(rt.x()+1.0f)/2.0f);
    int newy = int(YImageSize*(1.0f-(rt.y()+1.0f)/2.0f));
    return (T(newx, newy));
}
template <class T>
T ViewClusters::to2DObjectPoint(T posn)
{
    float newx = (2.0f*(float(posn.x())/float(XImageSize)) - 1.0f);
    float newy = (2.0f*(YImageSize-float(posn.y()))/YImageSize) - 1.0f;
    QVector4D rp = QVector4D(newx,newy,0.0,1.0);
    QVector4D rt = invtransform*rp;
    return (T(rt.x(), rt.y()));
}
void ViewClusters::CalculateScaleBar()
{
    Point_2 sy1 = to2DObjectPoint(Point_2(0.0,0.0));
    Point_2 sy2 = to2DObjectPoint(Point_2(0.0,8.0));
    float sbht=float(sy1.y()-sy2.y());
    Point_2 st = ScaleBarPosition;
    ScaleBar[0] = QVector2D(st.x(), st.y());
    ScaleBar[1] = QVector2D(st.x(), st.y() + sbht);
    ScaleBar[2] = QVector2D(st.x() + nScaleBarLength, st.y());
    ScaleBar[3] = QVector2D(st.x() + nScaleBarLength, st.y() + sbht);
    scalebar_vbo.bind();
    scalebar_vbo.allocate(ScaleBar, int(4 * sizeof(QVector2D)));
    scalebar_vbo.release();
    update();
}
void ViewClusters::initShaders()
{
    colclust = new QOpenGLShaderProgram();

    if(!colclust->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/denpnts.vert"))
    {
       emit AlertMsg("Error when compiling the point vertex shader",'r');
       emit AlertMsg(QString(colclust->log()),'r');
       return;
    }

    if(!colclust->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/clustcolor.frag"))
    {
       emit AlertMsg("Error when compiling the point fragment shader",'r');
       emit AlertMsg(QString(colclust->log()),'r');
       return;
    }
    colclust->bindAttributeLocation("blinks", 0);

    if(!colclust->link())
    {
       emit AlertMsg("Error when linking the point shader program",'r');
       emit AlertMsg(QString(colclust->log()),'r');
       return;
    }

    colclust->bind();

    vao_blinks.create();
    vao_blinks.bind();

    vsp_vbo.create();
    vsp_vbo.bind();
    vsp_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    loadBlinkData();
    std::size_t doff = offsetof(blinkdata,logDensity);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(blinkdata), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(blinkdata), reinterpret_cast<GLvoid*>(doff));
    vsp_vbo.release();

/*
    //  Cluster colours
    const QVector4D ClustColor[7] = {
        {1.0f,0.0f,0.0f,1.0f}, // red
        {0.0f,0.0f,1.0f,1.0f}, // blue
        {0.0f,1.0f,1.0f,1.0f}, // cyan
        {1.0f,0.0f,1.0f,1.0f}, // magenta
        {0.0f,1.0f,0.0f,1.0f}, // green
        {1.0f,1.0f,0.0f,1.0f}, // yellow
        {0.0f,0.0f,0.0f,1.0f}  // black
    };
    int CCol=colclust->uniformLocation("ClustColor");
    colclust->setUniformValueArray(CCol,ClustColor,7);
*/
    const QVector4D DensityColor[8] = {
        {1.0f, 1.0f, 1.0f, 1.0f},    //white
        {1.0f, 1.0f, 0.0f, 1.0f},    //yellow
        {0.0f, 1.0f, 1.0f, 1.0f},    //cyan
        {0.0f, 1.0f, 0.0f, 1.0f},    //green
        {1.0f, 0.65f, 0.0f, 1.0f},   //orange
        {1.0f, 0.0f, 0.0f, 1.0f},    //red
        {1.0f, 0.0f, 1.0f, 1.0f},    //magenta
        {0.0f, 0.0f, 0.0f, 1.0f}     //black
     };

    int DCol=colclust->uniformLocation("DensityColor");
    colclust->setUniformValueArray(DCol,DensityColor,8);

    PointSizeLoc=colclust->uniformLocation("PointSize");
    TransformColClustLoc=colclust->uniformLocation("transform");

    DenThresholdLoc=colclust->uniformLocation("denThreshold");
    ThresholdColourIndexLoc=colclust->uniformLocation("ThresholdColourIndex");

    vao_blinks.release();

    boundclust = new QOpenGLShaderProgram();

    if(!boundclust->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/boundary.vert"))
    {
       emit AlertMsg("Error when compiling the segment vertex shader",'r');
       emit AlertMsg(QString(boundclust->log()),'r');
       return;
    }

    if(!boundclust->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/clustcolor.frag"))
    {
        emit AlertMsg("Error when compiling the segment fragment shader",'r');
        emit AlertMsg(QString(boundclust->log()),'r');
        return;
    }
    boundclust->bindAttributeLocation("boundline", 0);
    boundclust->bind();

    vao_bounds.create();
    vao_bounds.bind();

    bound_vbo.create();
    bound_vbo.bind();
    bound_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    loadOuterSegmentData();

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), nullptr);
    bound_vbo.release();
    //  Scale & shift

    TransformBoundClustLoc=boundclust->uniformLocation("transform");
    BoundaryColourLoc=boundclust->uniformLocation("BoundaryColour");
    ExPointSizeLoc=boundclust->uniformLocation("ExPointSize");

    vao_bounds.release();

// Excluded blinks

    vao_exblinks.create();
    vao_exblinks.bind();

    exblinks_vbo.create();
    exblinks_vbo.bind();
    exblinks_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    loadExcludedBlinks();

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), nullptr);
    exblinks_vbo.release();
    vao_exblinks.release();

// Convex Hull

    vao_hull.create();
    vao_hull.bind();

    hull_vbo.create();
    hull_vbo.bind();
    hull_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    loadHullData();

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), nullptr);
    hull_vbo.release();
    vao_hull.release();

// Min Ellipse

    loadEllipseData();

    vao_ellipse.create();
    vao_ellipse.bind();

    ellipse_vbo.create();
    ellipse_vbo.bind();
    ellipse_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    ellipse_vbo.allocate(ellipsepts, int(totalEllipsePnts * sizeof(QVector2D)));

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), nullptr);
    ellipse_vbo.release();
    vao_ellipse.release();

// scalebar

    vao_scalebar.create();
    vao_scalebar.bind();

    scalebar_vbo.create();
    scalebar_vbo.bind();
    scalebar_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(QVector2D), nullptr);
    scalebar_vbo.release();
    vao_scalebar.release();

// tetramer plot

    tetplot = new QOpenGLShaderProgram();

    if(!tetplot->addShaderFromSourceFile(QOpenGLShader::Vertex, ":/tetramers.vert"))
    {
       emit AlertMsg("Error when compiling the tetramer vertex shader",'r');
       emit AlertMsg(QString(tetplot->log()),'r');
       return;
    }

    if(!tetplot->addShaderFromSourceFile(QOpenGLShader::Fragment, ":/clustcolor.frag"))
    {
       emit AlertMsg("Error when compiling the tetramer fragment shader",'r');
       emit AlertMsg(QString(tetplot->log()),'r');
       return;
    }
    tetplot->bindAttributeLocation("tetramers", 0);

    if(!tetplot->link())
    {
       emit AlertMsg("Error when linking the tetramer shader program",'r');
       emit AlertMsg(QString(tetplot->log()),'r');
       return;
    }

    tetplot->bind();

    vao_tetramers.create();
    vao_tetramers.bind();

    tetramers_vbo.create();
    tetramers_vbo.bind();
    tetramers_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    std::size_t toff = offsetof(tetramerimage,linecolour);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(tetramerimage), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(tetramerimage), reinterpret_cast<GLvoid*>(toff));

    TransformTetramerLoc=tetplot->uniformLocation("transform");
    TransparencyLoc=tetplot->uniformLocation("transparency");

    tetramers_vbo.release();
    vao_tetramers.release();

// single tetramer

    vao_singletetramer.create();
    vao_singletetramer.bind();

    singletetramer_vbo.create();
    singletetramer_vbo.bind();
    singletetramer_vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(tetramerimage), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(tetramerimage), reinterpret_cast<GLvoid*>(toff));

    singletetramer_vbo.release();
    vao_singletetramer.release();

/*
    Colours[0] = QVector4D(0.,0.,0.,1.); // black
    Colours[1] = QVector4D(0.902, 0.624, 0.0, 1.0); // orange
    Colours[2] = QVector4D(0.337, 0.706, 0.913, 1.0); // skyblue
    Colours[3] = QVector4D(0., 0.62, 0.451,1.0); // bluegreen
    Colours[4] = QVector4D(0.941, 0.894, 0.259, 1.0); //yellow
    Colours[5] = QVector4D(0., 0.447, 0.698, 1.0); // blue
    Colours[6] = QVector4D(0.835, 0.369, 0.,1.0); // vermillion
    Colours[7] = QVector4D(0.8, 0.475, 0.655, 1.0); // redpurple
*/
    CUD_Colours[0] = {0.f,0.f,0.f,1.f}; //black
    CUD_Colours[1] = {0.902f, 0.624f, 0.0f, 1.0f}; // orange
    CUD_Colours[2] = {0.337f, 0.706f, 0.913f, 1.0f}; // skyblue
    CUD_Colours[3] = {0.f, 0.62f, 0.451f,1.0f}; // bluegreen
    CUD_Colours[4] = {0.941f, 0.894f, 0.259f, 1.0f}; //yellow
    CUD_Colours[5] = {0.f, 0.447f, 0.698f, 1.0f}; // blue
    CUD_Colours[6] = {0.835f, 0.369f, 0.f,1.0f}; // vermillion
    CUD_Colours[7] = {0.8f, 0.475f, 0.655f, 1.0f}; // redpurple

}
void ViewClusters::paintGL()
{
    GLint tetline = GL_QUADS;
    float transparency=1.0f;
    QPainter p(PaintDevice);
    p.beginNativePainting();
    if(bIsolateTetramers)
    {
        glClearColor(0, 0, 0, 1);
        tetline = GL_LINE_LOOP;
    }
    else
    {
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glClearColor(1, 1, 1, 1);
        transparency = 0.5f;
    }
    glClear(GL_COLOR_BUFFER_BIT);
    glLineWidth(2.0f);
    if(!bIsolateTetramers)
    {
        glUseProgram(colclust->programId());
        vao_blinks.bind();
        if(SD->displayNoThreshold)
        {
            denThreshold = 3.0f;
            ThresholdColourIndex = 7;
        }
        else
        {
            denThreshold=float(SD->LogDensityThreshold);
            ThresholdColourIndex = 3;
        }
        colclust->setUniformValue(PointSizeLoc, float(PointSize));
        colclust->setUniformValue(DenThresholdLoc, denThreshold);
        colclust->setUniformValue(ThresholdColourIndexLoc, ThresholdColourIndex);
        colclust->setUniformValueArray(TransformColClustLoc,&transform, 1);

        glMultiDrawArrays(GL_POINTS,vs_start,vs_count,GLsizei(nValidClusters));
        vao_blinks.release();

        if(SD->displayBoundary || SD->displayConvexHull || SD->displayMinEllipse || bScaleBar || SD->displayExcludedBlinks)
        {
            glUseProgram(boundclust->programId());
            boundclust->setUniformValueArray(TransformBoundClustLoc,&transform, 1);
        }

        if(SD->displayExcludedBlinks)
        {
            vao_exblinks.bind();
            BoundaryColour=QVector4D(1.0,0.0,1.0,1.0);
            boundclust->setUniformValue(ExPointSizeLoc, float(ExPointSize));
            boundclust->setUniformValue(BoundaryColourLoc, BoundaryColour);
            glDrawArrays(GL_POINTS,0,nExBlinks);
            vao_exblinks.release();
        }


        if(SD->displayBoundary)
        {
            vao_bounds.bind();
            if(SD->displayNoThreshold) // || SD->displayConvexHull)
                BoundaryColour=QVector4D(1.0,0.0,0.0,1.0);
            else
                BoundaryColour=QVector4D(0.0,0.0,1.0,1.0);
            boundclust->setUniformValue(BoundaryColourLoc, BoundaryColour);
            glMultiDrawArrays(GL_LINES,bn_start,bn_count,GLsizei(nValidBoundarySegs));
            vao_bounds.release();
        }


        if(SD->displayConvexHull)
        {
            vao_hull.bind();
            BoundaryColour=QVector4D(0.0,0.0,1.0,1.0);
            boundclust->setUniformValue(BoundaryColourLoc, BoundaryColour);
            glMultiDrawArrays(GL_LINE_LOOP,hl_start,hl_count,GLsizei(nValidHulls));
            vao_hull.release();
        }

        if(SD->displayMinEllipse)
        {
            vao_ellipse.bind();
            BoundaryColour=QVector4D(0.0,1.0,0.0,1.0);
            boundclust->setUniformValue(BoundaryColourLoc, BoundaryColour);
            glMultiDrawArrays(GL_LINE_LOOP,el_start,el_count,GLsizei(nValidEllipses));
            vao_ellipse.release();
        }

        if(bScaleBar)
        {
            vao_scalebar.bind();
            ScaleBarColour = QVector4D(0.0,0.0,0.0,1.0);
            boundclust->setUniformValue(BoundaryColourLoc, ScaleBarColour);
            glDrawArrays(GL_TRIANGLE_STRIP,0,4);
            vao_scalebar.release();
        }
    }
    if(bPlotTetramers)
    {
        glUseProgram(tetplot->programId());
        tetplot->setUniformValueArray(TransformTetramerLoc,&transform, 1);
        tetplot->setUniformValue(TransparencyLoc, 1.0f);
        if(bPlacingTetramer || bModifyTetramer || bClassifyTetramer)
        {
            vao_singletetramer.bind();
            glDrawArrays(GL_LINE_LOOP,0,4);
            vao_singletetramer.release();
        }
        if(nImagedTetramers > 0)
        {
            tetplot->setUniformValue(TransparencyLoc, transparency);
            if(!bIsolateTetramers)
                glEnable(GL_BLEND);
            vao_tetramers.bind();
            glMultiDrawArrays(tetline,tl_start,tl_count,GLsizei(nImagedTetramers));
            vao_tetramers.release();
            if(!bIsolateTetramers)
                glDisable(GL_BLEND);
        }
    }

    p.endNativePainting();

    if(bIsolateTetramers)
    {
        NNDlinecolor = Qt::white;
        NNDtextcolor = Qt::white;
    }
    else
    {
        NNDlinecolor = Qt::black;
        NNDtextcolor = Qt::black;
    }

    if(bDisplayNND)
    {
        p.setPen(NNDlinecolor);
        p.setFont(basefont);
        for (uint i = 0; i < NNDtet.size(); i++)
        {
            Point_2 p1 = to2DViewportPoint(NNDtet[i].NNDStart);
            Point_2 p2 = to2DViewportPoint(NNDtet[i].NNDEnd);
            p.setPen(NNDlinecolor);
            p.drawLine(QPointF(p1.x(),p1.y()),QPointF(p2.x(),p2.y()));
            float xpos = ((p1.x()+p2.x())/2);
            float ypos = ((p1.y()+p2.y())/2);
            p.setPen(NNDtextcolor);
            p.drawText(xpos, ypos, QString::number(NNDtet[i].NND,'f',1));
        }
    }

    if(bLinePlot && EndPt != StartPt)
    {
        p.setPen(Qt::red);
        p.setFont(basefont);
        QString LineVal = "";

        p.drawLine(QPointF(StartPt.x(),StartPt.y()),QPointF(EndPt.x(),EndPt.y()));
        Point_2 StartL = to2DObjectPoint(StartPt);
        Point_2 EndL = to2DObjectPoint(EndPt);
        char mu = char(0xb5);
        double dXdist= double(EndL.x() - StartL.x());
        double dYdist= double(EndL.y() - StartL.y());
        double dDist = sqrt(dXdist*dXdist + dYdist*dYdist);

        if(dDist > 1000.0)
            LineVal = QString("%1 %2m").arg(QString::number(dDist/1000.0, 'g',3)).arg(mu);
        else
            LineVal = QString("%1 nm").arg(int(dDist));

        p.drawText(int(StartPt.x() + (EndPt.x() - StartPt.x())/2) , int(StartPt.y() + (EndPt.y() - StartPt.y())/2 + 20), LineVal);
    }

    if(SD->displayNumbers)
    {
        p.setPen(Qt::red);
        p.setFont(basefont);
        for(uint i=0; i < SD->Clusters.size(); i++)
        {
            Cluster a = SD->Clusters[i];
            uint cno = a.getClustNo();
            Point_2 posn = to2DViewportPoint(a.getLowestPoint());
            if(!SD->displayNoThreshold)
            {
              if(a.isWrongSize())
                  p.setPen(Qt::gray);
              else{
                  if(a.isBoundaryValid())
                     p.setPen(Qt::red);
                  else
                     p.setPen(Qt::blue);
                  }
            }
            if(a.isBoundaryValid() || SD->displayNoThreshold)
                p.drawText(posn.x(), posn.y() + 20, QString().setNum(cno));
        }
    }

    p.end();
}
void ViewClusters::loadBlinkData()
{
    if(vs_start != nullptr)
    {
       delete[] vs_start;
       delete[] vs_count;
       delete[] bd;
    }

    std::vector<bool> validcluster;
    validcluster.reserve(NoClusters);

    for(uint k=0; k < NoClusters; k++)
        validcluster.push_back(!SD->Clusters[k].isWrongSize() || SD->displayNoThreshold);

    uint totalBlinks = 0;
    nValidClusters=0;
    for(uint k=0; k < NoClusters; k++)
    {
        if(validcluster[k])
        {
           totalBlinks += SD->Clusters[k].getClustSize();
           nValidClusters++;
        }
    }
    vs_start = new GLint[nValidClusters];
    vs_count = new GLsizei[nValidClusters];

    bd = new blinkdata[totalBlinks];
    uint iv = 0;
    uint ic = 0;
    for(uint k=0; k < NoClusters; k++)
    {
        if (!(validcluster[k]))
            continue;
        std::vector<ClusterPoint> a = SD->Clusters[k].getClustPoints();
        uint nNoBlinks=uint(a.size());
        vs_start[ic] = GLint(iv);
        for(uint i=0; i < nNoBlinks; i++)
        {
          blinkdata q;
          q.logDensity = float(a[i].logdensity);
          q.vpos=QVector2D(a[i].posn.x(),a[i].posn.y());
          bd[iv++] = q;
        }
        vs_count[ic] = GLsizei(nNoBlinks);
        ic++;
    }

    emit AlertMsg(QString("Number of blinks entered = %1 - Total in set = %2").arg(iv).arg(totalBlinks),'b');
//    emit AlertMsg(QString("Max X = %1 - Max Y = %2").arg(maxx).arg(maxy), ' ');
    vsp_vbo.allocate(bd, int(totalBlinks * sizeof(blinkdata)));
}
void ViewClusters::loadExcludedBlinks()
{
    nExBlinks = uint(SD->ExcludedBlinks.size());
    QVector2D* bdex = new QVector2D[nExBlinks];
    for(int j=0; j < nExBlinks; j++)
    {
        bdex[j]= SD->ExcludedBlinks.at(j);
    }
    exblinks_vbo.allocate(bdex, int(nExBlinks*sizeof(QVector2D)));
    delete [] bdex;
}
void ViewClusters::loadOuterSegmentData()
{
    if(bn_start != nullptr)
    {
       delete[] bn_start;
       delete[] bn_count;
       delete[] boundline;
    }
    nValidBoundarySegs = 0;
    uint totalBoundaryPnts = 0;
    uint NoBaseClusters = uint(SD->BaseClusters.size());

//    emit AlertMsg(QString("No. of Clusters = %1 ").arg(NoClusters),'b');

    for(size_t k=0; k < NoBaseClusters; k++)
    {
        if (!SD->BaseClusters[k].isValid())
            continue;
        int nsegpts = SD->BaseClusters[k].getNumBoundarySegmentPts();
        totalBoundaryPnts += nsegpts;
        nValidBoundarySegs++;
    }
    bn_start = new GLint [size_t(nValidBoundarySegs)];
    bn_count = new GLsizei [size_t(nValidBoundarySegs)];
//    emit AlertMsg(QString("Total Points from boundary segments = %1").arg(totalBoundaryPnts),'b');
    boundline = new QVector2D[size_t(totalBoundaryPnts)];
    int ib=0;
    int j=0;
    for(uint k=0; k < NoBaseClusters; k++)
    {
        if (!(SD->BaseClusters[k].isValid()))
            continue;
        std::vector<QVector2D> obs = SD->BaseClusters[k].getBoundarySegmentPts();
        bn_start[j] = ib;
        bn_count[j] = GLsizei(obs.size());
        std::copy(obs.begin(),obs.end(), &boundline[ib]);
        ib += bn_count[j];
        j++;
    }

    bound_vbo.allocate(boundline, int(totalBoundaryPnts * sizeof(QVector2D)));
}
void ViewClusters::loadHullData()
{
    if(hl_start != nullptr)
    {
       delete[] hl_start;
       delete[] hl_count;
       delete[] hullpts;
    }
    uint totalHullPnts = 0;
    nValidHulls = 0;
    for(uint k=0; k < NoClusters; k++)
    {
        if(SD->Clusters[k].isWrongSize())
            continue;
        uint nHullPnts = SD->Clusters[k].getNumHullPoints();
        totalHullPnts += nHullPnts;
        nValidHulls++;
    }
    hl_start = new GLint [nValidHulls];
    hl_count = new GLsizei [nValidHulls];
    hullpts = new QVector2D[totalHullPnts];
    int ich=0;
    int j3=0;
    for(uint k=0; k < NoClusters; k++)
    {
       if (SD->Clusters[k].isWrongSize())
           continue;
       std::vector<Point_2> ihs = SD->Clusters[k].getConvexHullPoints();
       hl_start[j3] = ich;
       hl_count[j3] = GLsizei(ihs.size());
       for(uint j4=0; j4 < ihs.size(); j4++)
       {
          hullpts[ich] = QVector2D(ihs[j4].x(),ihs[j4].y());
          ich++;
       }
       j3++;
    }
    hull_vbo.allocate(hullpts, int(totalHullPnts * sizeof(QVector2D)));
}
void ViewClusters::loadEllipseData()
{
    if(el_start != nullptr)
    {
        delete[] el_start;
        delete[] el_count;
        delete[] ellipsepts;
    }
    totalEllipsePnts = 0;
    nValidEllipses = 0;
    for(uint k=0; k < NoClusters; k++)
    {
        if(SD->Clusters[k].isWrongSize())
            continue;
        int nEllipsePnts = SD->Clusters[k].getNumEllipsePoints();
        totalEllipsePnts += nEllipsePnts;
        nValidEllipses++;
    }
    el_start = new GLint [nValidEllipses];
    el_count = new GLsizei [nValidEllipses];
    ellipsepts = new QVector2D[totalEllipsePnts];
    uint ich=0;
    uint j3=0;
    for(uint k=0; k < NoClusters; k++)
    {
       if (SD->Clusters[k].isWrongSize())
           continue;
       std::vector<Point_2> ies = SD->Clusters[k].getEllipsePoints();
       el_start[j3] = ich;
       el_count[j3] = GLsizei(ies.size());
       for(uint j4=0; j4 < ies.size(); j4++)
       {
          ellipsepts[ich] = QVector2D(ies[j4].x(),ies[j4].y());
          ich++;
       }

       j3++;
    }
}
void ViewClusters::loadTetramerData()
{
    if(tl_start != nullptr)
    {
        delete[] tl_start;
        delete[] tl_count;
        delete[] tetrameroutline;
        tl_start = nullptr;
        tl_count = nullptr;
        tetrameroutline = nullptr;
    }
    nNoTetramers = uint(SD->tm.size());
    if(nNoTetramers == 0)
        return;
    nImagedTetramers = nNoTetramers;
    if(nBeingModified > -1 || nBeingClassified > -1)
        nImagedTetramers--;
    if(nImagedTetramers == 0)
        return;

    imageindex = 0;

    tl_start = new GLint[size_t(nImagedTetramers)];
    tl_count = new GLsizei[size_t(nImagedTetramers)];
    tetrameroutline = new tetramerimage[nImagedTetramers*4];
    QVector3D PlacingColor;
    QVector3D UnlabeledColor;
    QVector3D ChkrboardColor;
    QVector3D tcolor;
    if (nPalette == 1)
    {
        PlacingColor = QVector3D(1.0f,1.0f,1.0f);
        UnlabeledColor = QVector3D(0.7f,0.0f,0.0f);
        ChkrboardColor = QVector3D(0.0f,1.0f,0.0f);
    }
    else
    {
        PlacingColor = QVector3D(1.0f,0.0f,0.0f);
        UnlabeledColor = QVector3D(0.0f,0.7f,0.0f);
        ChkrboardColor = QVector3D(0.0f,1.0f,0.0f);
    }
    int j = 0;
    for(int i = 0; i < int(nNoTetramers); i++)
    {
        if(i == nBeingModified || i == nBeingClassified)
            continue;
        else
        {
            switch(SD->tm[i].classification)
            {
            case 'u':
                tcolor=UnlabeledColor;
                break;
            case 'c':
                tcolor=ChkrboardColor;
                break;
            case 's':
                tcolor=QVector3D(1.0,1.0,0.0);
                break;
            case 'i':
                tcolor=QVector3D(0.0,1.0,1.0);
                break;
            case 'b':
                tcolor=QVector3D(1.0,0.0,1.0);
                break;
            default:
                tcolor=PlacingColor;
            }
            calcTetramer(SD->tm[i],tcolor, tetrameroutline);
        }
        tl_start[j]=GLint(j*4);
        tl_count[j]=GLsizei(4);
        j++;
    }
    tetramers_vbo.bind();
    tetramers_vbo.allocate(tetrameroutline, int(4*j*sizeof(tetramerimage)));
    tetramers_vbo.release();
}
void ViewClusters::calcTetramer(Tetramer tmer, QVector3D tcolor, tetramerimage *tmi)
{
    double angleinradians = tmer.angle*degtorad+fortyfivedeg;

    float x1 = float(tetramerdiag * cos (angleinradians));
    float y1 = float(tetramerdiag * sin (angleinradians));
    float xpos = tmer.xcentre;
    float ypos = tmer.ycentre;
    tetramerimage td;
    td.vpos = QVector2D(xpos + x1, ypos + y1); td.linecolour = tcolor;
    tmi[imageindex++] = td;
    td.vpos = QVector2D(xpos - y1, ypos + x1);
    tmi[imageindex++] = td;
    td.vpos = QVector2D(xpos - x1, ypos - y1);
    tmi[imageindex++] = td;
    td.vpos = QVector2D(xpos + y1, ypos - x1);
    tmi[imageindex++] = td;
    return;
}
void ViewClusters::loadSingleTetramerData()
{
    QVector3D Placing;
    if (nPalette == 1)
       Placing = QVector3D(1.0,1.0,1.0);
    else
       Placing = QVector3D(0.0,0.0,0.0);

    imageindex = 0;
    calcTetramer(ntm,Placing, SingleTetramer);

    singletetramer_vbo.bind();
    singletetramer_vbo.allocate(SingleTetramer,4*sizeof(tetramerimage));
    singletetramer_vbo.release();
    update();
}
void ViewClusters::resizeGL(int width, int height)
{
    PaintDevice->setSize(QSize(width,height));
}
void ViewClusters::saveTetramerFile()
{
    if(SD->tm.size() == 0)
      return;
    if(nBeingModified > -1 ||  nBeingClassified > -1 || bPlacingTetramer)
    {
       emit AlertMsg("Cannot save Tetramer Data while adding or changing a tetramer",'r');
       return;
    }
    int nNoUnclassified = 0;
    int nNoCheckerboard = 0;
    int nNoSidebySide = 0;
    int nNoIsolated = 0;
    int nNoBoth = 0;
    emit AlertMsg(QString("Saving tetramer file: %1 ").arg(SD->TetramerFile),'b');
    QFile file(SD->TetramerFile);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        emit AlertMsg(QString("Tetramer save error: ").arg(file.errorString()),'r');
        return;
    }

    QTextStream out(&file);
    out << "Neighbourhood Limit = " << SD->NeighbourhoodLimit << "\n";

    for (uint jj = 0; jj < nNoTetramers; jj++)
    {
        Tetramer tout = SD->tm[jj];
        out << tout.xcentre << ' '<< tout.ycentre << ' ' << tout.angle << ' ' << tout.classification << ' ' << tout.ClusterNo <<"\n";
    }
    file.close();


    for (uint jj = 0; jj < nNoTetramers; jj++)
    {
        switch(SD->tm[jj].classification)
        {
        case 'u':
            nNoUnclassified++;
            break;
        case 'c':
            nNoCheckerboard++;
            break;
        case 's':
            nNoSidebySide++;
            break;
        case 'i':
            nNoIsolated++;
            break;
        case 'b':
            nNoBoth++;
            break;
        }
    }

    if(bCalcNearestNeighbours)
        calculateNearestNeighbours();

    QFile nnfile(SD->NNDFile);  // nearest neighbour distance file
    if(!nnfile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        emit AlertMsg(QString("NND file open error: %1").arg(nnfile.errorString()),'r');
        return;
    }
    QTextStream nnout(&nnfile);

    // output to both file and screen
    QString cs;
    uint nLabeledClusters = uint(Tetcount.size());
    for(uint jk = 0; jk < nLabeledClusters; jk++)
    {
       uint clustNo = TetMappedCluster[jk];
       double area = SD->Clusters[clustNo].getArea();
       uint noTets = Tetcount[clustNo];
       cs=QString("Cluster %1 has area %2 and contains %3 tetramers").arg(clustNo).arg(area).arg(noTets);
       nnout << qPrintable(cs) << "\n";
       emit AlertMsg(cs,'b');
    }
    cs = QString("Total tetramers = %1").arg(nNoTetramers);
    nnout << qPrintable(cs) << "\n";
    emit AlertMsg(cs,'b');
    if(nNoCheckerboard > 0)
    {
        cs=QString("No. checkerboard = %1").arg(nNoCheckerboard);
        nnout << qPrintable(cs) << "\n";
        emit AlertMsg(cs,'b');
    }
    if(nNoSidebySide > 0)
    {
        cs = QString("No. side by side = %1").arg(nNoSidebySide);
        nnout << qPrintable(cs) << "\n";
        emit AlertMsg(cs,'b');
    }
    if(nNoIsolated > 0)
    {
        cs = QString("No. Isolated = %1").arg(nNoIsolated);
        nnout << qPrintable(cs) << "\n";
        emit AlertMsg(cs,'b');
    }
    if(nNoBoth > 0)
    {
        cs = QString("No. Both = %1").arg(nNoBoth);
        nnout << qPrintable(cs) << "\n";
        emit AlertMsg(cs,'b');
    }
    if(nNoUnclassified > 0)
    {
        cs = QString("No. Unclassified = %1").arg(nNoUnclassified);
        nnout << qPrintable(cs) << "\n";
        emit AlertMsg(cs,'b');
    }
    if (NNDtet.size() > 0)
    {
        for(uint k= 0; k < NNDtet.size(); k++)
        {
            cs = QString("%1").arg(QString::number(NNDtet[k].NND, 'f',1));
            nnout << qPrintable(cs) << "\n";
        }

        nnfile.close();

        emit AlertMsg(QString("Saving nearest neighbour file : %1").arg(SD->NNDFile),'b');
    }
    bSavedTetFile = true;
}

void ViewClusters::Redraw(QString elements)
{
    NoClusters=uint(SD->Clusters.size());

    ctree.clear();
    setupClusterSearch();
    if(SD->reassignClusters)
    {
      reassignClusterNos();
      nLabelledClusters = ClusterCount();
    }

    if(elements.contains("Points"))
    {
        vsp_vbo.bind();
        loadBlinkData();
        vsp_vbo.release();

        exblinks_vbo.bind();
        loadExcludedBlinks();
        exblinks_vbo.release();
    }

    if(elements.contains("Outer"))
    {
        bound_vbo.bind();
        loadOuterSegmentData();
        bound_vbo.release();
        hull_vbo.bind();
        loadHullData();
        hull_vbo.release();
        loadEllipseData();
        ellipse_vbo.bind();
        ellipse_vbo.allocate(ellipsepts, int(totalEllipsePnts * sizeof(QVector2D)));
        ellipse_vbo.release();
    }

    update();
}
void ViewClusters::createMenus()
{
    showHelpAct = new QAction(tr("Help {F1}"),this);
    connect(showHelpAct, SIGNAL(triggered()), this, SIGNAL(showHelp()));

    resetDisplayAct = new QAction(tr("reset Display {F2}"),this);
    connect(resetDisplayAct, SIGNAL(triggered()), this, SLOT(resetDisplay()));

    ROIAct = new QAction("Create ROI {F3}",this);
    connect(ROIAct, SIGNAL(triggered()), this, SLOT(toggleROI()));

    clusterInfoAct = new QAction(tr("Cluster Info {F4}"),this);
    connect(clusterInfoAct, SIGNAL(triggered()), this, SLOT(clusterInfo()));

    addScaleBarAct = new QAction(tr("Add Scale Bar {F5}"),this);
    connect(addScaleBarAct, SIGNAL(triggered()), this, SLOT(getScaleBarSize()));

    IsolateTetramersAct = new QAction(tr("Isolate Tetramers (Hide Blinks) {F6}"),this);
    connect(IsolateTetramersAct, SIGNAL(triggered()), this, SLOT(isolateTetramers()));

    displayTetramersAct = new QAction(tr("Show Tetramers {F7}"),this);
    connect(displayTetramersAct, SIGNAL(triggered()), this, SLOT(displayTetramers()));

    changeTetramerWidthAct = new QAction(tr("Change Tetramer Width {F8}"), this);
    connect(changeTetramerWidthAct, SIGNAL(triggered()), this, SLOT(changeTetramerWidth()));

    saveTetramerFileAct = new QAction(tr("Save Tetramer File {F9}"),this);
    connect(saveTetramerFileAct, SIGNAL(triggered()), this, SLOT(saveTetramerFile()));

    increasePointSizeAct = new QAction(tr("Increase Point Size {F10}"),this);
    connect(increasePointSizeAct, SIGNAL(triggered()), this, SLOT(increasePointSize()));

    decreasePointSizeAct = new QAction(tr("Decrease Point Size {F11}"),this);
    connect(decreasePointSizeAct, SIGNAL(triggered()), this, SLOT(decreasePointSize()));

    createScreenShotAct = new QAction(tr("Save TIFF screenshot {F12}"),this);
    connect(createScreenShotAct, SIGNAL(triggered()), this, SLOT(CreateScreenShot()));

    createFBOImageAct = new QAction(tr("Save Hi-Res Zoomed image {Alt + F12}"),this);
    connect(createFBOImageAct, SIGNAL(triggered()), this, SLOT(CreateFBOImage()));

    LineAct = new QAction(tr("Measure Distance {l}"),this);
    connect(LineAct, SIGNAL(triggered()), this, SLOT(imageLine()));

    ZoomInAct = new QAction(tr("Zoom In {+}"),this);
    connect(ZoomInAct, SIGNAL(triggered()), this, SLOT(imageZoomIn()));

    ZoomOutAct = new QAction(tr("Zoom Out {-}"),this);
    connect(ZoomOutAct, SIGNAL(triggered()), this, SLOT(imageZoomOut()));

    QuitAct = new QAction(tr("Quit {esc}"),this);
    connect(QuitAct, SIGNAL(triggered()), this, SLOT(close()));

    addTetramerAct = new QAction(tr("Add {a}"),this);
    connect(addTetramerAct, SIGNAL(triggered()), this, SLOT(addTetramer()));

    modifyTetramerAct = new QAction(tr("Modify {m}"),this);
    connect(modifyTetramerAct, SIGNAL(triggered()), this, SLOT(modifyTetramer()));

    deleteTetramerAct = new QAction(tr("Delete {d}"),this);
    connect(deleteTetramerAct, SIGNAL(triggered()), this, SLOT(deleteTetramer()));

    saveTetramerAct = new QAction(tr("Save {s}"),this);
    connect(saveTetramerAct, SIGNAL(triggered()), this, SLOT(saveTetramer()));

    showNearestNeighboursAct = new QAction(tr("Show Nearest Neighbours {n}"),this);
    connect(showNearestNeighboursAct, SIGNAL(triggered()), this, SLOT(showNearestNeighbours()));

    imageMenu->addAction(showHelpAct);
    imageMenu->addSeparator();
    imageMenu->addAction(resetDisplayAct);
    imageMenu->addAction(addScaleBarAct);
    imageMenu->addSeparator();
    imageMenu->addAction(createScreenShotAct);
    imageMenu->addAction(createFBOImageAct);
    imageMenu->addSeparator();
    imageMenu->addAction(increasePointSizeAct);
    imageMenu->addAction(decreasePointSizeAct);
    imageMenu->addSeparator();
    imageMenu->addAction(ROIAct);
    imageMenu->addSeparator();
    imageMenu->addAction(clusterInfoAct);
    imageMenu->addSeparator();
    imageMenu->addAction(ZoomInAct);
    imageMenu->addAction(ZoomOutAct);
    imageMenu->addSeparator();

    tetramerMenu=imageMenu->addMenu(tr("Tetramers"));
    tetramerMenu->addAction(changeTetramerWidthAct);
    tetramerMenu->addAction(displayTetramersAct);
    tetramerMenu->addSeparator();
    tetramerMenu->addAction(addTetramerAct);
    tetramerMenu->addAction(deleteTetramerAct);
    tetramerMenu->addAction(modifyTetramerAct);
    tetramerMenu->addAction(saveTetramerAct);
    tetramerMenu->addSeparator();
    tetramerMenu->addAction(showNearestNeighboursAct);
    tetramerMenu->addAction(IsolateTetramersAct);
    tetramerMenu->addAction(saveTetramerFileAct);
}
void ViewClusters:: CreateTIFF(int XSize, int YSize, int imagetype)
{
   double ActualPixelSize = PixelSize;
   if(imagetype == BackBuffer)  // Screenshot
       emit AlertMsg(QString("Creating Screenshot image (%1 x %2) as TIFF ....").arg(XSize).arg(YSize),' ');
   else
   {
       ActualPixelSize = PixelSize/fFBOZoom;
       emit AlertMsg(QString("Creating FBO image (%1 x %2) as TIFF ....").arg(XSize).arg(YSize),' ');
   }
   int nNoRGBABytes = XSize*YSize*4;

   GLubyte* ScreenData=nullptr;
   try{
       ScreenData = new GLubyte[size_t(nNoRGBABytes)];
   }
   catch (std::bad_alloc& ba)
   {
       emit AlertMsg(QString("Unable to allocate memory for TIFF creation<br>Reason: %1").arg(ba.what()),'r');
       return;
   }
#ifdef WIN32
   QOpenGLContext* ctx;
   QSurface* ctsurf;
   if(imagetype == BackBuffer)
   {
     ctx = context();
     ctsurf = ctx->surface();
     ctx->swapBuffers(ctsurf);
   }
#endif

   glReadPixels(0, 0, XSize, YSize, GL_RGBA, GL_UNSIGNED_BYTE, ScreenData);
   GLenum err = glGetError();
   if (err != GL_NO_ERROR)
   {
      const char* errString = "something"; //reinterpret_cast<const char*> (gluErrorString (err));
      emit AlertMsg(QString("OpenGL Error  - %1 - %2 at glReadPixels").arg(err).arg(errString),'r');
      return;
   }
#ifdef WIN32
   if(imagetype == BackBuffer)
   {
       ctx->swapBuffers(ctsurf);
       makeCurrent();
   }
#endif
   QFileInfo fi(SD->FileName);
   fi.makeAbsolute();
   QDateTime whenTIFFcreated(QDateTime::currentDateTime());
   QString DateTimeInfo=whenTIFFcreated.toString("yyyy:MM:dd hh:mm:ss");

   QString fname=fi.filePath();
   fname.remove("."+fi.suffix());

   QFile qdata;
   QString OutputFile;
   int picno = 1;

   forever
   {
       OutputFile = QString("%1_p%2.tiff").arg(fname).arg(picno);
       qdata.setFileName(OutputFile);
       if(!qdata.exists())
           break;
       picno++;
   }

   QString Description = QString("RyR_STORM2D image : Pixel size = %1 nm").arg(ActualPixelSize,0,'f',1);
   if(SD->displayConvexHull || SD->displayExcludedBlinks || SD->displayMinEllipse)
     Description += QString("\nNeighbourhood Limit = %1 nm  - Minimum Blinks Per Cluster = %2 - Minimum Tetramer Size = %3 nm2").arg(SD->NeighbourhoodLimit).arg(SD->MinimumBlinksPerCluster).arg(SD->MinTetramerArea);
   if(!SD->displayNoThreshold)
     Description += QString("\nLog Density Threshold = %1").arg(SD->LogDensityThreshold);
   TiffImage Tiff_Image;
   Tiff_Image.setOutputType(TIFF_RGBA_8bit);
   Tiff_Image.setDim(uint32(XSize), uint32(YSize), 1);
   if(imagetype==BackBuffer)
       Tiff_Image.setRes(SD->ScreenDPI,2);
   else
       Tiff_Image.setRes(300.0f,2);
   Tiff_Image.setBytePointer(ScreenData);
   Tiff_Image.setDateTime(qPrintable(DateTimeInfo));
   Tiff_Image.setOrigin(0,0);
   Tiff_Image.setSoftware("RyR_STORM2D, Moore Lab, UBC, Vancouver, Canada");
   Tiff_Image.CalculateMinMax();
   Tiff_Image.setDescription(Description);
   if(Tiff_Image.WriteTiff(OutputFile))
       emit AlertMsg(QString("Image File: %1 saved").arg(OutputFile),'b');
   else
       emit AlertMsg(QString("Unable to save image file"),'r');

   return;
}
