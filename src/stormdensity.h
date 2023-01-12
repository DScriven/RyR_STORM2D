/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, stormdensity.h, is a header file for stormdensity.cpp and is part of the RyR_STORM2D program.
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
#ifndef STORMDENSITY_H
#define STORMDENSITY_H

#include "paramvals.h"
#include "param.h"
#include <QObject>
#include <QThread>
#include <QString>
#include <QElapsedTimer>
#include <list>
#include <string>
#include <vector>
#include "blinks.h"
#include "blinkvertex.h"
#include "tetramer.h"
#include <CGAL/Cartesian.h>
#include "cluster.h"
#include <utility>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef CGAL::Search_traits_2<K> TreeTraits2;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits2> Neighbor_search2;
typedef Neighbor_search2::Tree Tree2;

struct connections
{
    bool visitflag;
    int stpnt;
    std::vector<int> endpnts;
    connections(bool bv, int st, const std::vector<int>& s) : visitflag(bv), stpnt(st), endpnts(s) {}

    bool operator < (const connections& str) const
    {
        return (stpnt < str.stpnt);
    }
};

struct distances
{
    int stpnt;
    std::vector<std::pair<int,double>> edgelength;
    distances(int st, const std::vector<std::pair<int,double>>& s) : stpnt(st), edgelength(s) {}
};


class STORMdensity: public QThread
{
    Q_OBJECT

public:
    explicit STORMdensity();
    ~STORMdensity() override;
public slots:
    void CalculateStructure(QString File, ParamVals ParameterValues);
    void NewClusterValues(ParamVals ParameterValues, int FrameMin, int FrameMax);
    void writeROI(QString OutputFile, QRectF ROIXY);
    void setScreenParam(QSize ScrnSize, qreal Screendpi) {ScreenSize = ScrnSize; ScreenDPI = Screendpi;}
protected:
    bool read_blinks();
    void clearVectors();
    void CopyBlinkArray();
    bool FrameFilter();
    bool CalculateDensity();
    void CalculateClusters(std::vector<int>& v2, std::vector<Cluster>& Clust);
    void CalculateClusterProperties_1stpass();
    Cluster CheckConvexDisparity(Cluster a);
    void CalculateClusterProperties_2ndpass();
    void CreateClusterMap(std::vector<Cluster>& Clust, QHash <uint,int>& CMap, std::list<APoint>& CList);
    bool FindNearestCluster();
    void GetPointsinClusters(double Limit, std::vector<Cluster>& Clust, std::vector<distances>& edist);
    void secondPass();
    void secondPass_Clusters();
    void secondPass_NL30();
    bool CalculateGroupAreas();
    bool writeClusterData();
    bool readTetramers();

    std::vector<int> tracepoints(int i);

    std::vector<InputBlinks> dSTORMInput;
    std::vector<Blinks> dSTORMres;
    std::vector<Blinks> dSTORM;
    std::vector<blinkVertex> bv;
    std::vector<QVector2D> ExcludedBlinks;

    std::vector<connections> edgedata;
    std::vector<distances> edgedistance;
    std::vector<distances> edgedistance_2ndpass;
    QElapsedTimer* tic;

    std::vector<Cluster> BaseClusters;
    std::vector<Cluster> Clusters;
    QHash<uint, int> BaseCentroidMap;
    std::list<APoint> BaseCentroidList;
    QHash<uint, int> CentroidMap;
    std::list<APoint> CentroidList;
    std::vector<Tetramer> tm;

    QSize ScreenSize;
    qreal ScreenDPI;
    QRect ParamGeometry;
    QRectF Limits;
    QString FileName;
    QString TetramerFile;
    QString NNDFile;

    bool    bFrameFilter;
    bool    bUnMatched;

    int     nFrameMax;
    int     nFrameMin;
    int     nMaxFrameNo;
    int     nMinFrameNo;
    uint    BlinkCoordinates;

    double dXmin;
    double dXmax;
    double dYmin;
    double dYmax;
    double dZmin;
    double dZmax;
    double TotalClusterArea;

    double XInputOffset;
    double YInputOffset;

    double dMaxAmplitude;
    double dMinAmplitude;

    double Alpha1stPass;
    double ConvexAlphaRatio;

    bool bHaveProcessedFile;
    bool bSaveClusterData;
    bool reassignClusters;
    bool displayNumbers;
    bool displayBoundary;
    bool displayNoThreshold;
    bool displayExcludedBlinks;
    bool displayConvexHull;
    bool displayMinEllipse;
    bool Tetramersloaded;
    uint MinimumBlinksPerCluster;
    double NeighbourhoodLimit;
    double LogDensityThreshold;
    double MinTetramerArea;
    int PassStage;

    friend class ViewClusters;

public:
signals:
    void ShowImage();
    void RedrawImage(QString instruction);
    void AlertMsg(QString msg, char Color);
    void Progress(int value);
    void setFrameRange(int nMinFrameNo, int nMaxFrameNo);
    void closeROI();
};

#endif // STORMDENSITY_H
