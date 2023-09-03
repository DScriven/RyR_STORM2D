/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, viewclusters.h, is a header file for viewcluster.cpp, part of the RyR_STORM2D program
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
#ifndef VIEWCLUSTERS_H
#define VIEWCLUSTERS_H

#include <QtGui>
#include <QOpenGLWindow>
#include <QOpenGLFunctions_4_3_Compatibility>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QSurfaceFormat>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QMatrix4x4>
#include <QRectF>
#include <QColor>
#include <QString>
#include <QFont>
#include <vector>
#include <iostream>
#include <fstream>
#include "cluster.h"
#include "stormdensity.h"

//#pragma pack(push,1)
struct blinkdata{
  QVector2D vpos;
  float logDensity;
};
struct tetramerimage{
  QVector2D vpos;
  QVector3D linecolour;
};

//#pragma pack(pop)

class QOpenGLShaderProgram;
class QOpenGLFramebufferObject;
class QOpenGLPaintDevice;
class QAction;
class QMenu;
class QRubberBand;
class BlinkROI;

struct NNDval
{
    int ClustNo;
    double NND;
    Point_2 NNDStart;
    Point_2 NNDEnd;
};

class ViewClusters: public QOpenGLWindow, protected QOpenGLFunctions_4_3_Compatibility
{
    Q_OBJECT
public:
    ViewClusters(STORMdensity* SCalc);
    ~ViewClusters() override;
    void initializeGL () override;
    void paintGL() override;
    void resizeGL(int w, int h) override;
public slots:
    void Redraw(QString elements);
    void Display();
    void FBOTiff(double FBOMag, int dpi);
private:
    void keyPressEvent( QKeyEvent *e ) override;
    void mousePressEvent( QMouseEvent *e ) override;
    void mouseDoubleClickEvent( QMouseEvent *e ) override;
    void mouseReleaseEvent( QMouseEvent *e ) override;
    void mouseMoveEvent( QMouseEvent *e ) override;
    void wheelEvent(QWheelEvent *event) override;
    void validAddModifyKeys(QKeyEvent *e);
    void validClassifyKeys(QKeyEvent *e);
    void createMenus();
    void initShaders();
    void initValues();
    void loadBlinkData();
    void loadExcludedBlinks();
    void loadOuterSegmentData();
    void loadHullData();
    void loadEllipseData();
    void loadTetramerData();
    void loadSingleTetramerData();
    void CalculateScaleBar();
    void setupTransforms();
    void changeProjection(float xm, float ym);
    void changeTranslation(float xshift, float yshift);
    void changeScalingandTranslation(float xs, float ys);
    void changeScaling();
    void changeSize();
    void calculatePixelSize();
    void readTetramerFile();
    void reassignClusterNos();

    void CreateTIFF(int Xdim, int Ydim, int imagetype, int ImageDPI, double Magnification);

    Point_2 LimitMouse(int x, int y);
    template <class T> T to2DViewportPoint(T posn);
    template <class T> T to2DObjectPoint(T posn);

    void setupClusterSearch();
    void calcTetramer(Tetramer tmer, QVector3D tcolor, tetramerimage *tm);
    bool findTetramer(Point_2 posn);
    int  findCluster(Point_2 posn);
    void tetramerClassification(char c);
    void calculateArea();
    void calculateNearestNeighbours();
    void ClusterNearestNeighbours(std::vector<Tetramer> tmk);
    void displayClusterInfo(uint cno);
    uint ClusterCount();

    enum {BackBuffer, FrameBuffer};
    Tree2 ctree;

    QMenu* tetramerMenu;
    QMenu* imageMenu;

    QAction* showHelpAct;
    QAction* resetDisplayAct;
    QAction* addScaleBarAct;
    QAction* createScreenShotAct;
    QAction* createFBOImageAct;
    QAction* BlackBreakAct;
    QAction* increaseBlinkSizeAct;
    QAction* decreaseBlinkSizeAct;
    QAction* LineAct;
    QAction* ZoomInAct;
    QAction* ZoomOutAct;
    QAction* QuitAct;
    QAction* ROIAct;

    QAction* addTetramerAct;
    QAction* deleteTetramerAct;
    QAction* modifyTetramerAct;
    QAction* showTetramersAct;
    QAction* saveTetramerAct;
    QAction* displayTetramersAct;
    QAction* showNearestNeighboursAct;
    QAction* IsolateTetramersAct;
    QAction* changeTetramerWidthAct;
    QAction* saveTetramerFileAct;
    QAction* clusterInfoAct;

    float fXtrans;
    float fYtrans;
    float fSizeFactor;
    float fZoom;
    float fOldZoom;
    float fZoomDelta;
    Point_2 PosnWanted;
    int nScaleBarLength;

    uint nNoTetramers;
    uint nImagedTetramers;
    uint nValidClusters;
    uint nValidHulls;
    uint nValidEllipses;
    uint nValidBoundarySegs;
    uint nLabelledClusters;
    uint nExcludedPoints;

    uint NoClusters;
    int InitXImage;
    int InitYImage;
    int XImageSize;
    int YImageSize;

    bool bROI;
    bool bLinePlot;
    bool bLeftButtonDown;
    bool bDoubleClick;
    bool bScaleBar;
    bool bScaleBarPositionChosen;

    int nBeingModified;
    int nBeingClassified;
    uint imageindex;
    uint CurrentCluster;

    double degtorad;
    double tetramerdiag;
    double tetramerwidth;
    double baseangleinradians;
    double centredistance;
    double sqrt2;
    double fortyfivedeg;

    double fFBOZoom;
    int ImageDPI;

    std::vector<int> TetMappedCluster;
    std::map<uint,uint> Tetcount;
    std::vector<NNDval> NNDtet;
    Tetramer ntm;
    double PixelSize;

    BlinkROI* roiinfo;

    bool bModifyTetramer;
    bool bAddTetramer;
    bool bClassifyTetramer;
    bool bPlacingTetramer;
    bool bPlotTetramers;
    bool bDisplayNND;
    bool bShowTetramer;
    bool bPickedCentre;
    bool bCalcNearestNeighbours;
    bool bSavedTetFile;
    bool bIsolateTetramers;
    bool bDisplayClusterInfo;

    QColor NNDlinecolor;
    QColor NNDtextcolor;

    Point_2 StartPt;
    Point_2 EndPt;
    Point_2 StartImagePt;
    Point_2 EndImagePt;
    Point_2 Centre;
    Point_2 ImageOffset;

    float xshift;
    float yshift;

    float fXScreenMax;
    float fYScreenMax;

    float fXScale;
    float fYScale;

    QFont basefont;
    QVector4D BoundaryColour;

    GLint viewportparam[4];
    GLsizei nExBlinks;

    GLint* vs_start;
    GLsizei* vs_count;
    blinkdata* bd;

    GLint* bn_start;
    GLsizei* bn_count;
    QVector2D* boundline;

    GLint* hl_start;
    GLsizei* hl_count;
    QVector2D* hullpts;

    GLint* el_start;
    GLsizei* el_count;
    QVector2D* ellipsepts;
    uint totalEllipsePnts;

    GLint* tl_start;
    GLsizei* tl_count;
    tetramerimage* tetrameroutline;
    tetramerimage SingleTetramer[4];

    Point_2 ScaleBarPosition;
    QVector4D ScaleBarColour;
    QVector2D ScaleBar[4];

    int BlinkSize;
    int ExBlinkSize;
    float denThreshold;
    int DenThresholdLoc;
    int BlinkSizeLoc;
    int ExBlinkSizeLoc;

    int TransparencyLoc;
    int TransformColClustLoc;
    int TransformBoundClustLoc;
    int TransformTetramerLoc;
    int BelowThresholdColourIndexLoc;
    int BelowThresholdColourIndex;
    int AboveThresholdColourIndexLoc;
    int AboveThresholdColourIndex;
    int BoundaryColourLoc;
    int nPalette;
    int imagetype;

    STORMdensity* SD;

    QSurfaceFormat sformat;
    QOpenGLPaintDevice* PaintDevice;
    QRectF DisplayLimits;

    enum CUD{black, orange, skyblue, bluegreen, yellow, blue, vermillion, redpurple};
    QVector4D CUD_Colours[8];

    QMatrix4x4 orthoproj;
    QMatrix4x4 translation;
    QMatrix4x4 scaling;
    QMatrix4x4 transform;
    QMatrix4x4 invtransform;

    QOpenGLVertexArrayObject vao_blinks;
    QOpenGLVertexArrayObject vao_exblinks;
    QOpenGLVertexArrayObject vao_bounds;
    QOpenGLVertexArrayObject vao_hull;
    QOpenGLVertexArrayObject vao_ellipse;
    QOpenGLVertexArrayObject vao_support;
    QOpenGLVertexArrayObject vao_scalebar;
    QOpenGLVertexArrayObject vao_tetramers;
    QOpenGLVertexArrayObject vao_singletetramer;

    QOpenGLBuffer vsp_vbo;
    QOpenGLBuffer exblinks_vbo;
    QOpenGLBuffer bound_vbo;
    QOpenGLBuffer hull_vbo;
    QOpenGLBuffer ellipse_vbo;
    QOpenGLBuffer support_vbo;
    QOpenGLBuffer scalebar_vbo;
    QOpenGLBuffer tetramers_vbo;
    QOpenGLBuffer singletetramer_vbo;

    QOpenGLFramebufferObject* image_fbo;

    QOpenGLShaderProgram* boundclust;
    QOpenGLShaderProgram* colclust;
    QOpenGLShaderProgram* tetplot;

    QPoint winoffset;
    QRubberBand* rubberBand;

private slots:
    void imageLine();
    void resetDisplay();
    void getScaleBarSize();
    void CreateScreenShot();
    void CreateFBOImage();
    void imageZoomIn();
    void imageZoomOut();
    void increaseBlinkSize();
    void decreaseBlinkSize();
    void increaseExBlinkSize();
    void decreaseExBlinkSize();

    void addTetramer();
    void saveTetramer();
    void deleteTetramer();
    void modifyTetramer();
    void classifyTetramer();
    void clusterInfo();
    void saveTetramerFile();
    void changeTetramerWidth();
    void displayTetramers();
    void showNearestNeighbours();
    void isolateTetramers();
    void CloseWindow();
    void toggleROI();
    void hideROI();

public:

signals:
    void AlertMsg(QString msg, char colour);
    void CloseProgram();
    void showHelp();
    void writeROI(QString FileName, QRectF ROI);
};

#endif // ViewClusters_H
