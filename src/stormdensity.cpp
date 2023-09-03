/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022, 2023.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, stormdensity.cpp, is part of the RyR_STORM2D program and groups blink data into clusters.
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
#include "stormdensity.h"
#include "blinks.h"
#include "cluster.h"
#include "vertex_data_traits.h"
#include <QTextStream>
#include <QFileInfo>
#include <QStringList>
#include <QRegularExpression>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <cmath>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Bbox_2.h>
#include <QCoreApplication>
#include <QDir>
#include <tuple>

struct BlinkInfo
{
    Point_2 posn;
    int frame;
    double amp;
};

using namespace std;

typedef CGAL::Triangulation_vertex_base_with_info_2<blinkVertex, K>  Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                     Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                       DT;

typedef DT::Point                   Point;
typedef DT::Locate_type             Locate_type;
typedef DT::Face_handle             Cell_handle;
typedef DT::Face_iterator           Cell_iterator;
typedef DT::Finite_faces_iterator   Finite_cells_iterator;
typedef DT::Vertex_handle           Vertex_handle;
typedef CGAL::Timer                 Timer;

// Nearest Neighbour
typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits2> NN_incremental_search2;
typedef NN_incremental_search2::iterator NN_iterator2;
typedef NN_incremental_search2::Tree ITree2;

// A functor that returns true, iff the squared distance of a dD point is > 810000 (i.e 900 squared);
struct dist2_gt_val {
    bool operator()(const NN_iterator2& it) { return ((*it).second > 810000);  }
};
// An iterator that only enumerates dD points with distance < 900 (squared)
typedef CGAL::Filter_iterator<NN_iterator2, dist2_gt_val> NN_distance_iterator2;

// A functor that returns true, iff the squared distance of a dD point is > 1440000 (i.e 1200 squared);
struct dist2_from_centre {
    bool operator()(const NN_iterator2& it) { return ((*it).second > 1440000);  }
};
// An iterator that only enumerates dD points with distance < 550 (squared)
typedef CGAL::Filter_iterator<NN_iterator2, dist2_from_centre> NN_centroid_distance_iter;


STORMdensity::STORMdensity()
{
    bFrameFilter = false;
    bHaveProcessedFile = false;
    reassignClusters = false;
    tic = new QElapsedTimer;
}
STORMdensity::~STORMdensity()
{
    clearVectors();
}
void STORMdensity::clearVectors()
{
    std::vector<InputBlinks>().swap(dSTORMInput);
    std::vector<Blinks>().swap(dSTORMres);
    std::vector<Blinks>().swap(dSTORM);
    std::vector<blinkVertex>().swap(bv);
    std::vector<connections>().swap(edgedata);
    std::vector<distances>().swap(edgedistance);
    std::vector<Cluster>().swap(Clusters);
    std::vector<Tetramer>().swap(tm);
    std::vector<QVector2D>().swap(ExcludedBlinks);

    CentroidMap.clear();
    CentroidList.clear();
}
void STORMdensity::CalculateStructure(QString File, ParamVals ParameterValues)
{
    if(bHaveProcessedFile)
       clearVectors();
    ParamGeometry = ParameterValues.paramGeometry;
    MinimumBlinksPerCluster = ParameterValues.MinimumBlinksPerCluster;
    MinTetramerArea = ParameterValues.MinTetramerArea;
    NeighbourhoodLimit = ParameterValues.NeighbourhoodLimit;
    Alpha1stPass = ParameterValues.Alpha1stPass;
    LogDensityThreshold = ParameterValues.LogDensityThreshold;
    displayNoThreshold = ParameterValues.DispNoThreshold;
    displayExcludedBlinks = ParameterValues.DispExcludedData;
    displayBoundary = ParameterValues.DispBoundary;
    displayNumbers = ParameterValues.DispNumbers;
    displayConvexHull = ParameterValues.DispConvexHull;
    displayMinEllipse = ParameterValues.DispMinEllipse;
    bSaveClusterData = ParameterValues.SaveClusterData;

 // set cluster statics
    Cluster::MinBlinksPerCluster = MinimumBlinksPerCluster;
    Cluster::MinTetramerArea = MinTetramerArea;

    if(displayExcludedBlinks) // display all blinks
        emit AlertMsg("Showing all blinks including excluded ones",'b');

    FileName = File;
    if(!read_blinks())
       return;
    FrameFilter();
    CalculateDensity();
    emit AlertMsg(" Beginning 1st pass",'m');
    CalculateClusterProperties_1stpass();
    CalculateClusterProperties_2ndpass();

    TetramerFile = FileName;
    QFileInfo q(TetramerFile);
    TetramerFile.remove("." + q.suffix());
    NNDFile = TetramerFile;
    TetramerFile += QString("_tetramer.dat");
    NNDFile += QString("_nnd.dat");

    readTetramers();

    bHaveProcessedFile = true;
    emit ShowImage();
}
void STORMdensity::NewClusterValues(ParamVals ParameterValues, int FrameMin, int FrameMax)
{
    ParamGeometry = ParameterValues.paramGeometry;
    int nOldFrameMax = nFrameMax;
    int nOldFrameMin = nFrameMin;

    nFrameMax = FrameMax;
    nFrameMin = FrameMin;

    bool reDisplay=false;
    QString displayChanges = "";
    if(displayBoundary != ParameterValues.DispBoundary)
    {
       displayBoundary = ParameterValues.DispBoundary;
       displayChanges += QString("Outer Boundaries ")  + QString((displayBoundary) ? "ARE " : "ARE NOT ") + QString("being displayed<br>");
       reDisplay = true;
    }
    if(displayNumbers != ParameterValues.DispNumbers)
    {
       displayNumbers = ParameterValues.DispNumbers;
       displayChanges += QString("Cluster Numbers ")  + QString((displayNumbers) ? "ARE " : "ARE NOT ") + QString("being displayed<br>");
       reDisplay = true;
    }
    if(displayConvexHull != ParameterValues.DispConvexHull)
    {
       displayConvexHull = ParameterValues.DispConvexHull;
       displayChanges += QString("Convex Hulls ")  + QString((displayConvexHull) ? "ARE " : "ARE NOT ") + QString("being displayed<br>");
       reDisplay = true;
    }
    if(displayMinEllipse != ParameterValues.DispMinEllipse)
    {
       displayMinEllipse = ParameterValues.DispMinEllipse;
       displayChanges += QString("Min. Ellipses ")  + QString((displayMinEllipse) ? "ARE " : "ARE NOT ") + QString("being displayed<br>");
       reDisplay = true;
    }
    if(displayNoThreshold != ParameterValues.DispNoThreshold)
    {
        QString exstr;
        displayNoThreshold = ParameterValues.DispNoThreshold;
        if (displayNoThreshold)
            exstr  = QString("Showing blinks without threshold<br>");
        else
        {
            LogDensityThreshold = ParameterValues.LogDensityThreshold;
            exstr  = QString("Blink log density threshold = %1<br>").arg(LogDensityThreshold);
        }
        displayChanges += exstr;
        reDisplay = true;
    }
    else
    {
        if(abs(LogDensityThreshold - ParameterValues.LogDensityThreshold) > DBL_EPSILON)
        {
            LogDensityThreshold = ParameterValues.LogDensityThreshold;
            QString exstr  = QString("Blink log density threshold = %1<br>").arg(LogDensityThreshold);
            displayChanges += exstr;
            reDisplay = true;
        }
    }
    if(displayExcludedBlinks != ParameterValues.DispExcludedData)
    {
        QString exstr;
        displayExcludedBlinks = ParameterValues.DispExcludedData;
        if(displayExcludedBlinks)
          exstr=QString("Showing excluded blinks<br>");
        else
          exstr=QString("Omitting excluded blinks<br>");

        displayChanges += exstr;
        reDisplay = true;
    }
    if(reDisplay)
    {
        QString ChangeMsg = QString("-------------------------------------------<br>Display changes:<br>") + displayChanges + QString("<br>-------------------------------------------");
        emit AlertMsg(ChangeMsg,'m');
    }

    bool bChangeFilter = false;
    if(nOldFrameMax != nFrameMax || nOldFrameMin != nFrameMin)
    {
        if((nFrameMin == nMinFrameNo) && (nFrameMax == nMaxFrameNo))
        {
            emit AlertMsg(QString("Removing Frame filter"),'b');
            bChangeFilter = true;
            bFrameFilter = false;
        }
        else
        {
            emit AlertMsg(QString("Changing Frame filter value range: %1 to %2").arg(nFrameMin).arg(nFrameMax),'b');
            bFrameFilter = true;
        }
    }
    else
    {
        if((nFrameMin == nMinFrameNo) && (nFrameMax == nMaxFrameNo))
            bFrameFilter = false;
        else
            bFrameFilter = true;
    }

    bool bRecalcBaseClusters = false;
    bool bRecalcClusters = false;

    reassignClusters = false;

    QString Changes = "";
    if(bSaveClusterData != ParameterValues.SaveClusterData)
    {
        bSaveClusterData = ParameterValues.SaveClusterData;
        Changes += (QString("Cluster Data ") + QString((bSaveClusterData) ? "IS " : "IS NOT ") + QString("being saved<br>"));
    }
    if(std::abs(NeighbourhoodLimit - ParameterValues.NeighbourhoodLimit) > DBL_EPSILON)
    {
        NeighbourhoodLimit = ParameterValues.NeighbourhoodLimit;
        Changes += QString("Cluster Distance = %1 nm<br>").arg(NeighbourhoodLimit);
        bRecalcClusters = true;
        if(tm.size() > 0)
        {
          Changes += QString("Clusters will be reassigned");
          reassignClusters = true;
        }
    }
    if(std::abs(Alpha1stPass - ParameterValues.Alpha1stPass) > DBL_EPSILON)
    {
        Alpha1stPass = ParameterValues.Alpha1stPass;
        Changes += QString("1st pass Alpha Value = %1 nm<br>").arg(Alpha1stPass);
        bRecalcBaseClusters = true;
    }

    if(std::abs(MinTetramerArea - ParameterValues.MinTetramerArea) > DBL_EPSILON)
    {
        MinTetramerArea = ParameterValues.MinTetramerArea;
        if(std::abs(MinTetramerArea) > DBL_EPSILON)
        {
            Changes += QString("Minimum Tetramer Area = %1 nm2<br>").arg(MinTetramerArea);
            // set cluster statics
            Cluster::MinTetramerArea = MinTetramerArea;
            bRecalcBaseClusters = true;
        }
    }

    if(MinimumBlinksPerCluster != ParameterValues.MinimumBlinksPerCluster)
    {
        MinimumBlinksPerCluster = ParameterValues.MinimumBlinksPerCluster;
        // set cluster statics
        Cluster::MinBlinksPerCluster = MinimumBlinksPerCluster;
        Changes += QString("Minimum Blinks Per Cluster = %1<br>").arg(MinimumBlinksPerCluster);
        bRecalcBaseClusters = true;
    }


    if(std::abs(LogDensityThreshold - ParameterValues.LogDensityThreshold) > DBL_EPSILON)
    {
        LogDensityThreshold = ParameterValues.LogDensityThreshold;
        Changes += QString("Log Density Threshold = %1<br>").arg(LogDensityThreshold);
    }

    if(bFrameFilter || bChangeFilter)
    {
       emit AlertMsg(QString("Recalculating everything because input blinks have changed:"),'g');
       if(!FrameFilter())
       {
          emit AlertMsg("Failure in FrameFilter - aborting",'r');
          return;
       }
       if(!Changes.isEmpty())
           emit AlertMsg(Changes,'g');
       CalculateDensity();
       CalculateClusterProperties_1stpass();
       CalculateClusterProperties_2ndpass();
       emit RedrawImage("Points,Outer");
       return;
    }

    if(Changes.isEmpty())
    {
        if(reDisplay)
            emit RedrawImage("");
        return;
    }
    else
    {
        QString ChangeMsg = QString("-------------------------------------------<br>Parameter changes:<br>") + Changes + QString("<br>-------------------------------------------");
        emit AlertMsg(ChangeMsg,'m');
    }

    if(bRecalcBaseClusters)
    {
        CalculateClusterProperties_1stpass();
        if(tm.size() > 0)
          reassignClusters = true;
        bRecalcClusters = true;
    }
    if(bRecalcClusters)
    {
        std::vector<Cluster>().swap(Clusters);
        CalculateClusterProperties_2ndpass();
        emit RedrawImage("Points Outer");
        return;
    }
    QString RedrawType;
    if(displayExcludedBlinks)
      RedrawType = "Points";

    if(bSaveClusterData)
       writeClusterData();

    if(!RedrawType.isEmpty() || reDisplay)
        emit RedrawImage(RedrawType);
}
void STORMdensity:: CalculateClusterProperties_1stpass()
{
    PassStage = Passtype::First;
    ConvexAlphaRatio = 3.0;
    uint nNoClusters = uint(BaseClusters.size());
    emit AlertMsg(QString("Base Clusters omitted:"),'m');

    uint badBlinks = 0;
    for(uint k=0; k < nNoClusters; k++)
    {
        Cluster a=BaseClusters[k];
        if(a.hasTooFewBlinks())
        {
           emit AlertMsg(QString("Base Cluster # %1 - %2").arg(k).arg(a.getStatus()),'r');
           badBlinks += a.getClustSize();
           continue;
        }
// need to delete fiduciary beads that might be hiding in a cluster
// cluster is deleted if #Blinks > 5 x #Points
        double blinkratio = double(a.getNoBlinks())/double(a.getClustSize());
        if(blinkratio > 5)
        {
            a.setStatus("Blink ratio too high - bead?");
            emit AlertMsg(QString("Base Cluster # %1 - %2").arg(k).arg(a.getStatus()),'r');
            continue;
        }
        a.CalculateCentroid();
        a.CalculateConvexHull();
        a.CalculateBoundaries(Alpha1stPass);
        double oarea = a.getArea();
        double carea = a.getConvexHullArea();
        if(a.isBoundaryValid() && (ConvexAlphaRatio*oarea  < carea))
        {
            emit AlertMsg(QString("Convex Disparity - BaseCluster %1 - Con Area = %2 Alpha Area = %3 ").arg(a.getClustNo()).arg(carea).arg(oarea),'m');
            a = CheckConvexDisparity(a);
        }
        if(a.isBoundaryValid() && oarea > 0.0)
        {
            if (oarea < MinTetramerArea)
            {
               a.setTooSmall(true);
               badBlinks += a.getClustSize();
               emit AlertMsg(QString("Base Cluster # %1 - %2").arg(k).arg(a.getStatus()),'r');
            }
            else if (oarea < 2000)
            {
               a.CalculateMinEllipse();
               if(a.hasBadShape())
               {
                  badBlinks += a.getClustSize();
                  emit AlertMsg(QString("Base Cluster # %1 - %2").arg(k).arg(a.getStatus()),'r');
               }
            }
        }
        else
            emit AlertMsg(QString("Base Cluster # %1 - %2").arg(k).arg(a.getStatus()),'r');

        BaseClusters[k]=a;
    }

    uint CountExcluded = 0;
    std::vector<QVector2D>().swap(ExcludedBlinks);
    for(uint k=0; k < nNoClusters; k++)
    {
        Cluster a=BaseClusters[k];
        if(!a.isValid())
        {
           std::vector<ClusterPoint> cp = a.getClustPoints();
           uint nNoBlinks=uint(cp.size());
           CountExcluded += nNoBlinks;
           for(uint i=0; i < nNoBlinks; i++)
           {
             QVector2D ex = QVector2D(cp[i].posn.x(), cp[i].posn.y());
             ExcludedBlinks.push_back(ex);
           }
           continue;
        }
    }
    CreateClusterMap(BaseClusters, BaseCentroidMap, BaseCentroidList);
    emit AlertMsg(QString("Bad Blinks = %1; Excluded = %2").arg(badBlinks).arg(CountExcluded),'r');
    emit AlertMsg(QString("Exclude Blinks - vector size = %1").arg(ExcludedBlinks.size()),'r');
}
Cluster STORMdensity::CheckConvexDisparity(Cluster b)
{
// Corrects a disparity between convex hull and alpha shape area. This is usually due to  a too
// low alpha value when you have a thin neck connecting the small part of the cluster to the
// larger part - the small part gets snipped off and thats the area you get, instead of the full
// cluster area - corrected by increasing the alpha until the disparity is resolved

    double oarea = ConvexAlphaRatio*b.getArea();
    double alpha = Alpha1stPass;
    double carea = b.getConvexHullArea();
    bool bCorrected = false;
    for(int k=2; k < 11 ; k++)
    {
        alpha *= k;
        b.CalculateBoundaries(alpha);
        oarea = b.getArea();
        if ((ConvexAlphaRatio*oarea) > carea)
        {
           emit AlertMsg(QString("Convex Disparity corrected - Con Area = %1 Alpha Area = %2").arg(carea).arg(oarea),'m');
           bCorrected = true;
           break;
        }
     }
     if(!bCorrected)
     {
        b.setStatus(QString("Uncorrected disparity: Alpha area = %1 - Convex Hull area = %2 - alpha = %3").arg(oarea).arg(carea).arg(alpha));
        b.setHasBadShape(true);
        b.setBoundaryValidity(false);
     }
     return b;
}
void STORMdensity::CreateClusterMap(std::vector<Cluster>& Clust, QHash <uint,int>& CMap, std::list<APoint>& CList)
{
    CList.clear();
    CMap.clear();
    uint nNoClusters = uint(Clust.size());
    uint maxx = uint(Limits.right());
    for (uint k =0; k < nNoClusters; k++)
    {
        Cluster a = Clust[k];
        if(a.isWrongSize())
           continue;
        uint clno=a.getClustNo();
        Point_2 cp2 = a.getCentroid();
//        QString cst = QString("%1 %2 - Centroid = (%3, %4)").arg(CType).arg(clno).arg(cp2.x()).arg(cp2.y());
//        emit AlertMsg (cst,'g');
//        std::cerr << qPrintable(cst) << "\n";
        APoint cp = APoint(cp2.x(),cp2.y());
        CList.push_back(cp);
        uint CvalOffset = uint(cp.y() * maxx + cp.x());
        CMap.insert(CvalOffset, clno);
    }
}
bool STORMdensity:: FindNearestCluster()
{
    int maxx = int(Limits.right());
    uint nNoClusters = uint(Clusters.size());
    if(nNoClusters < 2)
      return true;

    //   emit AlertMsg(QString("Calculating Nearest Neighbour of Centroids...."),' ');
    //   emit AlertMsg(QString("Found %1 Centroids").arg(Centroid.size()),'g');

    // now examine the centroid set for nearest neighbours
    std::vector<int> NearClusters;

    Tree2 ctree(CentroidList.begin(), CentroidList.end());
    ITree2 itree2(CentroidList.begin(), CentroidList.end());
    std::list<APoint>::iterator l;
    int centroidno = 0;

    emit AlertMsg("Cluster Nearest neighbour Edge to Edge distances",'m');
    for(l=CentroidList.begin(); l!=CentroidList.end(); ++l)
    {
        APoint centroidquery(*l);
        centroidno++;
        if(nNoClusters > 4)
        {
            // Find centroids within 900 nm of this one
            NN_incremental_search2 NN(itree2, centroidquery);
            NN_distance_iterator2 it2(NN.end(), dist2_gt_val(), NN.begin()), end(NN.end(), dist2_gt_val());
            it2++;
            NearClusters.clear();
            while(it2 != end)
            {
                APoint NearbyCentroid = (*it2).first;
                uint cOffset = uint((NearbyCentroid.y()) * maxx + NearbyCentroid.x());
                int cno = CentroidMap.value(cOffset, -1);
                if(cno == -1)
                {
                    emit AlertMsg("Cannot find centroid in map",'r');
                    return false;
                }
                NearClusters.push_back(cno);
                it2++;
            }
            // if two or less - look further  get nearest 3
            if(NearClusters.size() < 3)
            {
                NearClusters.clear();
//                emit AlertMsg(QString("Centroid at (%1, %2)- cannot find nearby clusters < 900 nm away<br>Trying nearest neighbours").arg(centroidquery.x()).arg(centroidquery.y()),'r');
                // nearest cluster centroid is greater 800 nm then look for nearest neighbours
                Neighbor_search2 searchc(ctree,centroidquery, 4);
                Neighbor_search2::iterator ittc = searchc.begin();
                for(int iz=0; iz < 3; iz++)
                {
                    ittc++;
                    APoint NearC = (*ittc).first;
                    uint cOffset = uint((NearC.y()) * maxx + NearC.x());
                    int cno = CentroidMap.value(cOffset, -1);
                    if(cno == -1)
                    {
                        emit AlertMsg("Cannot find centroid in map",'r');
                        return false;
                    }
                    NearClusters.push_back(cno);
                }
            }
        }
        else
        {
            NearClusters.clear();
            std::list<APoint>::iterator l2;
            for(l2=CentroidList.begin(); l2!=CentroidList.end(); ++l2)
            {
                if(l2 != l)
                {
                    APoint newC(*l2);
                    uint cOffset = uint((newC.y()) * maxx + newC.x());
                    int cno = CentroidMap.value(cOffset, -1);
                    if(cno == -1)
                    {
                        emit AlertMsg("Cannot find centroid in map",'r');
                        return false;
                    }
                    NearClusters.push_back(cno);
                }
            }
        }
        // create query set from the initial the cluster you are looking at
        int cqno = int(centroidquery.y()) * maxx + int(centroidquery.x());
        int ThisCluster = CentroidMap.value(cqno, -1);
        if(ThisCluster == -1)
        {
            emit AlertMsg("Cannot find query cluster in centroid map",'r');
            return false;
        }

        Cluster CurrentCluster = Clusters[size_t(ThisCluster)];
        uint nNoComPoints = CurrentCluster.getClustSize();
        std::vector<ClusterPoint> ct = CurrentCluster.getClustPoints();
        uint nNearClusters = uint(NearClusters.size());
        if(nNearClusters == 0)
        {
            emit AlertMsg(QString("No near clusters for cluster %1").arg(ThisCluster),'r');
            CurrentCluster.setEdgetoEdgeDist(1.e30);
            CurrentCluster.setNearestClusterArea(0.0);
            Clusters[ThisCluster] = CurrentCluster;
            continue;
        }
        // create nearest cluster tree for each cluster
        std::vector <tuple<double, double, int>> etoe;
        for(uint k=0; k < nNearClusters; k++)
        {
            double etoe_nndsq = 1.0e30;
            uint cno = NearClusters[k];
            uint nNoPoints=Clusters[cno].getClustSize();

            Tree2 treec;
            std::vector<ClusterPoint> cp = Clusters[cno].getClustPoints();
            for(uint j = 0; j < nNoPoints; j++)
            {
                APoint val = APoint(cp[j].posn.x(),cp[j].posn.y());
                treec.insert(val);
            }
            for(uint jk = 0; jk < nNoComPoints; jk++)
            {
                APoint query2 = APoint(ct[jk].posn.x(),ct[jk].posn.y());
                Neighbor_search2 search2(treec, query2, 1);
                Neighbor_search2::iterator itt = search2.begin();
                etoe_nndsq = qMin(itt->second, etoe_nndsq);
            }
            treec.clear();
            cp.clear();
            etoe.push_back({sqrt(etoe_nndsq), Clusters[cno].getArea(), cno});
        }
        double etoe_no = 1.e30;
        double nearArea = 0.0;
        int ncno;
        for(const auto &i : etoe)
        {
            double etoe_val = get<0>(i);
            if(etoe_val < etoe_no)
            {
                etoe_no = etoe_val;
                nearArea = get<1>(i);
                ncno = get<2>(i);
            }
        }
        CurrentCluster.setEdgetoEdgeDist(etoe_no);
        CurrentCluster.setNearestClusterArea(nearArea);
        QString etoe_str = QString::number(etoe_no,'f',0);
        QString nearArea_str = QString::number(nearArea,'f',0);
        emit AlertMsg(QString("Clust %1 : %2 nm to Clust %3 (area = %4 nm<sup>2</sup>)").arg(ThisCluster).arg(etoe_str).arg(ncno).arg(nearArea_str),'b');
        Clusters[ThisCluster] = CurrentCluster;
        ct.clear();
    }

    if(bSaveClusterData)
        return writeClusterData();
    return true;
}
bool STORMdensity::readTetramers()
{
    QFileInfo dfi(TetramerFile);
    if(!dfi.exists())
        return false;

    QFile file(TetramerFile);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        emit AlertMsg(QString("Error - %1").arg(file.errorString()),'r');
        return false;
    }
    int loadedLimit = 0;
    QTextStream input(&file);
    forever
    {
        QString tline=input.readLine();
        if (tline.isEmpty())
            break;

        if(tline.contains("Neighbourhood Limit ="))
        {
            loadedLimit=(tline.remove("Neighbourhood Limit =")).toInt();
//            emit AlertMsg(QString("Tetramer file loaded NLimit : %1").arg(loadedLimit),'b');
        }
        else
        {
            QStringList tl=tline.split(QRegularExpression("\\s+"));
            Tetramer itm;
            itm.xcentre = tl.at(0).toFloat();
            itm.ycentre = tl.at(1).toFloat();
            itm.angle = tl.at(2).toInt();
            itm.classification = tl.at(3).at(0).toLatin1();
            itm.ClusterNo = tl.at(4).toInt();
            tm.push_back(itm);
        }
    }
    file.close();
    uint nT = uint(tm.size());
    if(nT == 0)
    {
        emit AlertMsg("No valid tetramers in input",'r');
        return false;
    }
    else
    {
        emit AlertMsg(QString("%1 mapped tetramers in input").arg(nT),'g');
        if(loadedLimit != int(NeighbourhoodLimit))
            reassignClusters = true;
        else
            reassignClusters = false;

        return true;
    }
}
bool STORMdensity::read_blinks()
{
    QString STORMPath = FileName;
    if(!STORMPath.isEmpty())
    {
        QFileInfo q(STORMPath);
        QString exten = q.suffix();
        if(exten != "3dSTM")
        {
            emit AlertMsg(QString("%1 - not a valid file extension").arg(exten),'r');
            return false;
        }
    }
    else
    {
        return false;
    }
    tic->start();
    emit AlertMsg(QString("Reading STORM data %1").arg(STORMPath),' ');
    char* posn_buffer = nullptr;
    uint bufsize = 0;
    ifstream inFile(qPrintable(STORMPath), std::ifstream::binary);
    if(inFile)
    {
        // get length of file:
        inFile.seekg (0, inFile.end);
        bufsize = inFile.tellg();
        inFile.seekg (0, inFile.beg);
        posn_buffer = new char[bufsize+10];
        inFile.read(posn_buffer, bufsize);
        inFile.close();
    }
    else
    {
        emit AlertMsg(QString("Error - Can't open STORM data file %1").arg(STORMPath),'r');
        return false;
    }
//    cerr << "Reading in file took " << tic->elapsed() << " ms\n";

    char* pb = posn_buffer;
    int nNoBlinks = 0;
    for(uint i=0; i < bufsize; i++)
    {
        if( *pb == '\n')
            nNoBlinks++;
        pb++;
    }
    emit AlertMsg(QString("File contains %1 blinks").arg(nNoBlinks),' ');

    int nInBlinks = 0;
    double x=0.0;
    double y=0.0;
    double z=0.0;
    double dX=0.0;
    double dY=0.0;
    double dZ=0.0;
    double A=0.0;
    int picno=0;
    tic->restart();
    QString storminfo;
    QByteArray pdata = QByteArray::fromRawData(posn_buffer, bufsize);
    QTextStream pin(&pdata, QIODevice::ReadOnly);

    double maxdXe = 0;
    double maxdYe = 0;
    double maxdZe = 0;
    double maxAmpe = 0;
    double minAmpe = 1.e30;

    int MinFrameNo = 100000000;
    int MaxFrameNo = 0;
    while (nInBlinks < nNoBlinks)
    {
        storminfo = pin.readLine();
        nInBlinks++;
        QTextStream bin(&storminfo, QIODevice::ReadOnly);
        bin >> x >> y >> z >> dX >> dY >> dZ >> A >> picno;
        int perRead = (nInBlinks*100)/nNoBlinks;
        if(perRead % 5 == 0)
            emit Progress(perRead);

        InputBlinks a;
        a.x=x;
        a.y=y;
        a.z=z;
        a.dX = dX;
        a.dY = dY;
        a.dZ = dZ;
        maxdXe = (dX > maxdXe) ? dX : maxdXe;
        maxdYe = (dY > maxdYe) ? dY : maxdYe;
        maxdZe = (dZ > maxdZe) ? dZ : maxdZe;
        a.Amplitude = A;
        maxAmpe = (A > maxAmpe) ? A : maxAmpe;
        minAmpe = (A < minAmpe) ? A : minAmpe;
        a.pic_no=picno;
        MinFrameNo = (picno < MinFrameNo) ? picno : MinFrameNo;
        MaxFrameNo = (picno > MaxFrameNo) ? picno : MaxFrameNo;
        dSTORMInput.push_back(a);
    }
    emit Progress(-1);
    emit setFrameRange(MinFrameNo, MaxFrameNo);
    nMaxFrameNo = MaxFrameNo;
    nMinFrameNo = MinFrameNo;
    nFrameMin = MinFrameNo;
    nFrameMax = MaxFrameNo;

    delete [] posn_buffer;
    emit AlertMsg(QString("STORM data was read in %1 s").arg(double(tic->elapsed()) /1000.0),' ');

    emit AlertMsg(QString("Blinks input : %1").arg(nInBlinks),'b');

    if(dSTORMInput.size() < 50)
    {
        emit AlertMsg("Too few points - abandoned",'r');
        return false;
    }

    CopyBlinkArray();

    dSTORM = dSTORMres;

    return true;
}
void STORMdensity::CopyBlinkArray()
{
    uint nNoBlinks= uint(dSTORMInput.size());
    emit AlertMsg(QString("No of Input blinks = %1").arg(nNoBlinks),'b');

    dSTORMres.clear();
    for (uint j = 0; j < nNoBlinks; j++)
    {
        InputBlinks a = dSTORMInput.at(j);
        Blinks b;
        b.x = a.x;
        b.y = a.y;
        b.z = a.z;
        b.dX = a.dX;
        b.dY = a.dY;
        b.dZ = a.dZ;
        b.pic_no = a.pic_no;
        b.Amplitude = a.Amplitude;
        dSTORMres.push_back(b);
    }
}
bool STORMdensity::FrameFilter()
{
    std::vector<blinkVertex>().swap(bv);
    std::vector<connections>().swap(edgedata);
    std::vector<distances>().swap(edgedistance);
    std::vector<Cluster>().swap(Clusters);
    std::vector<QVector2D>().swap(ExcludedBlinks);
    if(!bFrameFilter)
    {
        dSTORM = dSTORMres;
    }
    else
    {
       std::vector<Blinks>().swap(dSTORM);
       int nFrameExcluded = 0;
       for (unsigned int j = 0; j < dSTORMres.size(); j++)
       {
          Blinks a = dSTORMres.at(j);
          if(bFrameFilter)
          {
            if(a.pic_no < nFrameMin || a.pic_no > nFrameMax)
            {
                nFrameExcluded++;
                continue;
            }
          }
          dSTORM.push_back(a);
       }
       if (nFrameExcluded > 0)
       {
         if(bFrameFilter && nFrameExcluded > 0)
            emit AlertMsg(QString("Frame excluded blinks = <br> %1 - frame numbers &lt. %2 or &gt. %3").arg(nFrameExcluded).arg(nFrameMin).arg(nFrameMax),'r');
       }
    }
    if(dSTORM.size() < 50)
    {
        emit AlertMsg(" Too few blinks to generate a usable image!",'r');
        return false;
    }
    return true;
}
bool STORMdensity::CalculateDensity()
{
    emit AlertMsg("Calculating point density",'m');
    uint nBlinks = uint(dSTORM.size());

    double minX = 1.e30;
    double minY = 1.e30;
    for (uint j=0; j < nBlinks; j++)
    {
        Blinks a = dSTORM.at(j);
        minX = (a.x < minX) ? a.x : minX;
        minY = (a.y < minY) ? a.y : minY;
    }

    XInputOffset = minX-10;
    YInputOffset = minY-10;
//    emit AlertMsg(QString("Min X = %1 - Min Y = %2").arg(minX).arg(minY),' ');

    std::vector<BlinkInfo> blinkposn;

    double minnewX = 1.e30;
    double minnewY = 1.e30;
    double maxnewX = 0;
    double maxnewY = 0;
    for (uint j=0; j < nBlinks; j++)
    {
        Blinks a = dSTORM.at(j);
        double newX =(a.x-XInputOffset);
        double newY =(a.y-YInputOffset);
        minnewX = (newX < minnewX) ? newX : minnewX;
        minnewY = (newY < minnewY) ? newY : minnewY;
        maxnewX = (newX > maxnewX) ? newX : maxnewX;
        maxnewY = (newY > maxnewY) ? newY : maxnewY;
        BlinkInfo b;
        b.posn = Point_2(newX,newY);
        b.frame = a.pic_no;
        b.amp = a.Amplitude;
        blinkposn.push_back(b);
    }

    emit AlertMsg(QString("New Min X = %1 - New Min Y = %2").arg(minnewX).arg(minnewY),' ');
    emit AlertMsg(QString("New Max X = %1 - New Max Y = %2").arg(maxnewX).arg(maxnewY),' ');
    Limits= QRectF(minnewX,minnewY,maxnewX,maxnewY);

    std::stable_sort(blinkposn.begin(),blinkposn.end(), [&](auto &a, auto &b){
        return  (a.posn.x() < b.posn.x()) ? true : (a.posn.x() == b.posn.x()) ?  a.posn.y() < b.posn.y() : false;});

    decltype(blinkposn.size()) ixy=0;

    std::vector<blinkVertex>().swap(bv);

    ixy = 0;
    int ndup = 0;
    uint j = 0;
    uint xylimit=uint(blinkposn.size())-1;
    std::vector<double> tAmp;
    std::vector<int> tFrame;
    while (ixy < xylimit)
    {
        Point_2 p1 = blinkposn[ixy].posn;
        Point_2 p2 = blinkposn[ixy+1].posn;
        if (p2 != p1)
        {
            j++;
            blinkVertex v;
            v.setPosition(p1);
         v.setWeight(double(ndup+1));
         if(ndup > 0)
         {
            for(uint k = 0; k < tAmp.size(); k++)
            {
                v.addAmp(tAmp[k]);
                v.addFrame(tFrame[k]);
            }
            tAmp.clear();
            tFrame.clear();
         }
         v.addAmp(blinkposn[ixy].amp);
         v.addFrame(blinkposn[ixy].frame);
            bv.push_back(v);
            ndup = 0;
        }
        else
        {
           tAmp.push_back(blinkposn[ixy].amp);
           tFrame.push_back(blinkposn[ixy].frame);
            ndup++;
        }
        ixy++;
        if (ixy == xylimit)
        {
           blinkVertex v;
           if (ndup == 0)
           {
             j++;
             v.setPosition(p2);
           }
           else
            {
              v.setPosition(p1);
              for(uint k = 0; k < tAmp.size(); k++)
              {
                v.addAmp(tAmp[k]);
                v.addFrame(tFrame[k]);
              }
              tAmp.clear();
              tFrame.clear();
            }
            v.addAmp(blinkposn[ixy].amp);
            v.addFrame(blinkposn[ixy].frame);
            v.setWeight(double(ndup+1));
           bv.push_back(v);
        }
    }

    if(j !=  bv.size())
    {
       emit AlertMsg(QString("Mismatch : Member count = %1 of blinkvertex  does not match vector size = %2").arg(j).arg(bv.size()),'r');
       return false;
    }
    else
       emit AlertMsg(QString("%1 Points after multiple blinks assigned to points").arg(j),' ');

    decltype(bv.size()) ibv=0;

/*    std::stable_sort( bv.begin(), bv.end() );

    std::cerr << "Stable sort\n";
    for(uint k=0; k < 20; k++)
    {
        blinkVertex v =bv[k];
        std::cerr << "1st bv #" << k << " " << v.position(0) << " " << v.position(1) << "\n";
    }
*/
    CGAL::spatial_sort( bv.begin(), bv.end(), blinkVertex_sort_traits() );

    for(uint k=0; k < bv.size(); k++)
    {
        blinkVertex v =bv[k];
        v.setVertexNo(k);
        bv[k]=v;
    }

    j = 0;
    DT dt;
    Vertex_handle vh;     // vertex handle - points to each inserted point
    for (ibv=0; ibv!=bv.size(); ++ibv)
    {
        blinkVertex v = bv[ibv];
        Point r = Point(v.position(0), v.position(1));
        vh = dt.insert(r);
        j++;
        vh->info().setData(v);      // set vertex quantities
    }

    emit AlertMsg(QString("The triangulation has %1 vertices").arg(dt.number_of_vertices()),'b');

    double factor = 3;

    std::vector<distances>().swap(edgedistance);
    std::pair<int, double> edgeidsquare;
    std::vector<pair<int,double>> endpointdist;
    DT::Finite_vertices_iterator vIT;

    for (vIT = dt.finite_vertices_begin(); vIT != dt.finite_vertices_end(); ++vIT )
    {
        std::vector<Cell_handle> cells;
        DT::Face_circulator fc = dt.incident_faces( vIT );
        cells.push_back( fc++ );
        for (; fc!=cells[0]; ++fc)
            cells.push_back( fc );

        double vol = 0.;
        double dens = 0;
        bool infinite_volume = false;
        for ( std::vector< Cell_handle >::const_iterator itC = cells.begin(); itC!=cells.end(); ++itC )
        {
            if ( !dt.is_infinite(*itC) )
                  vol += dt.triangle( *itC  ).area();
            else
            {
                vIT->info().setDummyNeighbor();
                infinite_volume = true;
                break;
            }
        }

        double w = vIT->info().weight();
        uint pntno= vIT->info().VertexNo();
        // compute the density
        if ( !infinite_volume )
        {
            dens = w * factor /vol ;
            bv[pntno].setDensity(dens);
        }
        else
        {
            bv[pntno].setDensity( -7 );
        }

    }

    std::vector<distances> edgedistance_1stpass;
    for (vIT = dt.finite_vertices_begin(); vIT != dt.finite_vertices_end(); ++vIT )
    {
        uint pntno= vIT->info().VertexNo();
        double vxpos = vIT->point().x();
        double vypos = vIT->point().y();
        Point_2 st=Point_2(vxpos,vypos);

        endpointdist.clear();
        // Now look at edges coming from each vertex
        Vertex_handle ivs;
        DT::Vertex_circulator vc = dt.incident_vertices(vIT);
        DT::Vertex_circulator start = vc;

        do
        {
           ivs=vc;
           Point_2 edend=Point_2(ivs->point().x(),ivs->point().y());
           double sqdist = CGAL::squared_distance(st,edend);
           edgeidsquare = make_pair(ivs->info().VertexNo(),sqdist);
           endpointdist.push_back(edgeidsquare);
           vc++;
        } while(vc != start);

        edgedistance_1stpass.emplace_back(pntno,endpointdist);
    }
    GetPointsinClusters(30.0, BaseClusters, edgedistance_1stpass);  // identify single clusters use
    return true;
}

void STORMdensity:: GetPointsinClusters(double NLimit, std::vector<Cluster>& Clust, std::vector<distances>& edist)
{
   std::vector<distances>().swap(edgedistance);
   edgedistance = edist;
   double sqmax = NLimit*NLimit;
   std::vector<connections>().swap(edgedata);
   std::vector<int> validvertex;
   std::vector<int> endpoint;
   decltype(edgedistance.size()) ie=0;
   for(ie=0; ie != edgedistance.size(); ++ie)
   {
       endpoint.clear();
       std::vector<std::pair<int,double>> a=edgedistance[ie].edgelength;
       for(size_t j=0; j < a.size(); j++)
       {
           if(a[j].second <= sqmax)
              endpoint.push_back(a[j].first);
       }
        uint vITposn=edgedistance[ie].stpnt;
        bool bVisitpoint = true;
        if(!(endpoint.size() > 0 && bv[vITposn].density() > -7))
            bVisitpoint = false;
        edgedata.emplace_back(bVisitpoint,vITposn,endpoint);
        validvertex.push_back(vITposn);
    }
    CalculateClusters(validvertex, Clust );
}
void STORMdensity::CalculateClusterProperties_2ndpass()
{
    emit AlertMsg(" Beginning 2nd pass",'m');

    if(abs(NeighbourhoodLimit - 30.0) < DBL_EPSILON)
        secondPass_NL30();
    else
        secondPass();

    FindNearestCluster();

    if(bUnMatched)
      emit AlertMsg("Warning - unmatched clusters on 2nd pass - total cluster area will be inaccurate!",'r');

    double AreaFactor = 1.0;
    uint bdec = 0;
    QString AreaUnit = "nm<sup>2</sup>";
    if(TotalClusterArea > 300000)
    {
        AreaFactor = 1.0e6;
        AreaUnit = "&mu;m<sup>2</sup>";
        bdec = 3;
    }
    double ImageArea = Limits.width()*Limits.height();
    double Occupancy = 100.0*TotalClusterArea/ImageArea;
    uint BlinksIn = uint(dSTORMInput.size());
    double ClusterDensity = AreaFactor*Clusters.size()/ImageArea;
    double BlinkRatio = double(BlinksIn)/double(BlinkCoordinates);
    double BlinkDensity = AreaFactor*BlinksIn/ImageArea;
    TotalClusterArea /= AreaFactor;
    emit AlertMsg(QString("<br>Input blinks = %1; Blink coordinates = %2; <br> Ratio = %3; Blink Density = %4/%5").arg(BlinksIn).arg(BlinkCoordinates).arg(BlinkRatio,0,'f',2).arg(BlinkDensity,0,'f',2).arg(AreaUnit),'m');
    emit AlertMsg(QString("Image dimensions: X = %1 nm - Y = %2 nm").arg(Limits.width()).arg(Limits.height()),'m');
    emit AlertMsg(QString("Total Area in clusters = %1 %2 <br> Area/ImageArea = %3 %; Cluster Density = %4/%5").arg(TotalClusterArea,0,'f',bdec).arg(AreaUnit).arg(Occupancy,0,'f',2).arg(ClusterDensity,0,'f',2).arg(AreaUnit),'m');
}
void STORMdensity::secondPass()
{
    PassStage = Passtype::Second;
    DT dt2;
    Vertex_handle vh2;     // vertex handle - points to each inserted point

    std::vector<blinkVertex>().swap(bv);

    int vno = 0;
    for(uint i=0; i < BaseClusters.size(); i++)
    {
        if(BaseClusters[i].isValid())
        {
            std::vector<ClusterPoint> a = BaseClusters[i].getClustPoints();
            for(uint k=0; k < a.size(); k++)
            {
               blinkVertex v;
               v.setVertexNo(vno++);
               v.setWeight(a[k].NoBlinks);
               v.setDensity(a[k].density);
               v.setPosition(a[k].posn);
               v.setFrame(a[k].frame);
               v.setAmp(a[k].amp);
               bv.push_back(v);
            }
       }
    }

    decltype(bv.size()) ibv=0;

    for (ibv=0; ibv!=bv.size(); ++ibv)
    {
        blinkVertex v = bv[ibv];
        Point r = Point(v.position(0), v.position(1));
        vh2 = dt2.insert(r);
        vh2->info().setData(v);      // set vertex quantities
    }

    emit AlertMsg(QString("The 2nd pass triangulation has %1 vertices").arg(dt2.number_of_vertices()),'b');

    std::vector<int> validvertex;
    std::vector<distances>().swap(edgedistance_2ndpass);
    std::pair<int, double> edgeidsquare;
    std::vector<pair<int,double>> endpointdist;
    DT::Finite_vertices_iterator vIT2;

    for (vIT2 = dt2.finite_vertices_begin(); vIT2 != dt2.finite_vertices_end(); ++vIT2 )
    {
        double vxpos = vIT2->point().x();
        double vypos = vIT2->point().y();
        Point_2 st=Point_2(vxpos,vypos);
        int pntno= vIT2->info().VertexNo();
        validvertex.push_back(pntno);
        endpointdist.clear();
        // Now look at edges coming from each vertex
        Vertex_handle ivs;
        DT::Vertex_circulator vc = dt2.incident_vertices(vIT2);
        DT::Vertex_circulator start = vc;

        do
        {
           ivs=vc;
           Point_2 edend=Point_2(ivs->point().x(),ivs->point().y());
           double sqdist = CGAL::squared_distance(st,edend);
           edgeidsquare = make_pair(ivs->info().VertexNo(),sqdist);
           endpointdist.push_back(edgeidsquare);
           vc++;
        } while(vc != start);

        edgedistance_2ndpass.emplace_back(pntno,endpointdist);
   }
   GetPointsinClusters(NeighbourhoodLimit, Clusters, edgedistance_2ndpass);
   secondPass_Clusters();
}
void STORMdensity::secondPass_Clusters()
{
    CreateClusterMap(Clusters, CentroidMap, CentroidList);
    CalculateGroupAreas();
}
void STORMdensity::secondPass_NL30()
{
    TotalClusterArea = 0.;
    bUnMatched = false;
    uint j = 0;
    std::vector<Cluster>().swap(Clusters);
    for(uint i=0; i < BaseClusters.size(); i++)
    {
        Cluster a = BaseClusters[i];
        if(!a.isWrongSize() && a.isBoundaryValid())
        {
           uint bcno = a.getClustNo();
           double bcarea = a.getArea();
           a.addSubClusterInfo(bcno, bcarea);
           a.setClustNo(j++);
           a.setArea(bcarea);
           TotalClusterArea += bcarea;
           a.CalculateMinEllipse();
           a.CalculateLowestPoint();
           a.CalculateDensityMeanMedian();
           emit AlertMsg(QString("Cluster # %1 - Area = %2 nm<sup>2</sup> from basecluster %3 ").arg(j).arg(bcarea).arg(bcno),'b');
           Clusters.push_back(a);
        }
    }
    CreateClusterMap(Clusters, CentroidMap, CentroidList);
}
bool STORMdensity::CalculateGroupAreas()
{
    TotalClusterArea = 0.0;
    int maxx = int(Limits.right());
    std::list<uint> VBCNos; // short hand for ValidBaseClusterNos
    for(uint i=0; i < BaseClusters.size(); i++)
    {
        if(BaseClusters[i].isValid())
           VBCNos.push_back(BaseClusters[i].getClustNo());
    }
    uint nNoClusters = uint(Clusters.size());

// we want to sort the clusters by Convex Hull Area so that when can populate the small ones first
// this ensures that the hulls are are properly populated even if a small hull is conatined within a large one
    std::vector <std::pair<uint,double>> Hullsort;
    for(uint k=0; k < nNoClusters; k++)
      Hullsort.push_back(make_pair(Clusters[k].getClustNo(), Clusters[k].getConvexHullArea()));

    std::sort(Hullsort.begin(), Hullsort.end(), [](auto &left, auto &right) {return left.second < right.second;});


    //   emit AlertMsg(QString("Calculating Nearest Neighbour of Centroids...."),' ');
    //   emit AlertMsg(QString("Found %1 Centroids").arg(Centroid.size()),'g');

    // now examine the centroid set for nearest neighbours
    QString baselist;
    std::vector<std::pair<uint, double>> baseinfo;
    ITree2 itree_bc(BaseCentroidList.begin(), BaseCentroidList.end());

    emit AlertMsg("Cluster Areas computed from sum of subcluster areas ",'m');
    bUnMatched = false;
    for(uint l=0; l < nNoClusters; l++)
    {
        uint ClustNo = Hullsort[l].first;
        Cluster a = Clusters[ClustNo];
        Point_2 ClustCentre = a.getCentroid();
//        emit AlertMsg(QString("Cluster %1 : Centroid = (%2, %3); Area = %4 nm<sup>2</sup>").arg(a.getClustNo()).arg(ClustCentre.x()).arg(ClustCentre.y()).arg(a.getArea()),'m');
//        if(a.getArea() < 600)
//           continue;
        Polygon_2 cb = a.getConvexHull();
        APoint centroidquery = APoint(ClustCentre.x(), ClustCentre.y());
        double clusterarea = 0;
       // Find centroids within 250 nm of this one
        NN_incremental_search2 NN(itree_bc, centroidquery);
        NN_centroid_distance_iter cit(NN.end(), dist2_from_centre(), NN.begin()), end(NN.end(), dist2_from_centre());
        baselist.clear();
        baseinfo.clear();
        while(cit != end)
        {
            APoint NearbyCentroid = (*cit).first;
            uint cOffset = uint((NearbyCentroid.y()) * maxx + NearbyCentroid.x());
            int cno = BaseCentroidMap.value(cOffset, -1);
            if(cno == -1)
            {
                emit AlertMsg("Cannot find centroid in map",'r');
                return false;
            }
//            emit AlertMsg(QString("Nearby Cluster %1 : Centroid = (%2, %3)").arg(cno).arg(NearbyCentroid.x()).arg(NearbyCentroid.y()),'m');
            if(cb.has_on_bounded_side(NearbyCentroid))
            {
               auto iter = std::find(VBCNos.begin(), VBCNos.end(), cno);
               if(iter == VBCNos.end())
                 emit AlertMsg(QString("Cluster %1 - Base Cluster %2 already assigned").arg(ClustNo).arg(cno),'r');
               else
               {
                 VBCNos.erase(iter);
                 baselist.append(QString(" %1").arg(cno));
                 double bcarea = BaseClusters[uint(cno)].getArea();
                 baseinfo.push_back(std::make_pair(uint(cno),bcarea));
                 clusterarea += bcarea;
               }
            }
            cit++;
        }
        if(abs(clusterarea - 0.0) < DBL_EPSILON)
            emit AlertMsg(QString("Error - Cluster %1 - unable to match base clusters - area set to zero").arg(ClustNo),'r');
        else
        {
            Clusters[ClustNo].addSubClusterInfo(baseinfo);
            Clusters[ClustNo].setArea(clusterarea);
            TotalClusterArea += clusterarea;
            emit AlertMsg(QString("Cluster %1 has area = %2 nm<sup>2</sup> from baseclusters %3").arg(ClustNo).arg(clusterarea).arg(baselist),'b');
        }
    }
    uint vs = uint(VBCNos.size());
    if(vs > 0)
    {
       bUnMatched = true;
       emit AlertMsg("Area error - the following base clusters were not used",'r');
       for (uint i : VBCNos)
       {
         Cluster a = BaseClusters[i];
         emit AlertMsg(QString("Base Cluster %1 - area = %2 nm<sup>2</sup>").arg(i).arg(a.getArea()),'b');
       }
    }
    return true;
}
bool STORMdensity::writeClusterData()
{
    //Output data set
    QFileInfo fn(FileName);
    QString fnameb = fn.canonicalPath() + "/" + fn.completeBaseName();
    QString fnout = fnameb + "_NL" + QString().number(int(NeighbourhoodLimit)) + ".clustval";
    QFile rdata;
    rdata.setFileName(fnout);
    if(!rdata.open(QFile::WriteOnly))
    {
        emit AlertMsg(QString("Cannot open %1 for writing").arg(fnout),'r');
        return false;
    }
    emit AlertMsg("Writing cluster areas and nearest neighbour distances...",' ');

    size_t noC=Clusters.size();
    uint maxNoSub = 0;
    for(size_t i1=0; i1 < noC; i1++)
    {
        uint NumSubClust = Clusters[i1].getNoSubClusters();
        maxNoSub = (maxNoSub > NumSubClust) ? maxNoSub : NumSubClust;
    }

    QTextStream rout(&rdata);
    rout << "Neighbourhood Limit = " << NeighbourhoodLimit << " nm";
    rout << " - Minimum Blinks Per Cluster = " << MinimumBlinksPerCluster;
    rout << " - Minimum Tetramer Area = " << MinTetramerArea << " nm2";
    rout << " - Alpha = " << Alpha1stPass << "\n";
    rout << "No. of Clusters = " << noC << " - Max No of SubClusters = " << maxNoSub << "\n";

    for(size_t i1=0; i1 < noC; i1++)
    {
        Cluster c = Clusters[i1];
        uint NumSubClust = c.getNoSubClusters();
        std::vector<std::pair<uint, double>> sc = c.getBaseClusterInfo();
        std::tuple<double,double> el = c.getEllipseData();
        rout << c.getClustNo() << " " <<  c.getNoBlinks() << " " <<  c.getClustSize() << " " <<  c.getArea() << " "
             << c.getNND() << " " << c.getNearestArea() << " " << std::get<0>(el) << " " << std::get<1>(el) << " "
             << c.getConvexHullArea() << " " << NumSubClust;
        for (uint j = 0; j < NumSubClust; j++)
            rout << " " << sc[j].second;
        rout <<  "\n";
    }
    rdata.close();
    emit AlertMsg(QString("Wrote %1 to disk").arg(fnout),' ');

    return true;
}
void STORMdensity::writeROI(QString OutputFile, QRectF ROIXY)
{
    int ROIXStart = qRound(ROIXY.left()+XInputOffset);
    int ROIYStart = qRound(ROIXY.bottom()+YInputOffset);
    int ROIXEnd = qRound(ROIXY.right()+XInputOffset);
    int ROIYEnd = qRound(ROIXY.top()+YInputOffset);

    emit AlertMsg(QString("STORM:ROI (pixels) : ( %1, %2 ) ( %3 %4 )").arg(ROIXStart).arg(ROIYStart).arg(ROIXEnd).arg(ROIYEnd),'b');

    emit AlertMsg("Generating Blink output",' ');
    std::vector<InputBlinks> ROISTORM;
    for(uint j = 0; j < dSTORMInput.size(); j++)
    {
        InputBlinks b = dSTORMInput[j];
        if(b.x < ROIXStart || b.x > ROIXEnd || b.y < ROIYStart || b.y > ROIYEnd)
            continue;
        ROISTORM.push_back(b);
    }
    uint nROIBlinks = uint(ROISTORM.size());
    if(nROIBlinks == 0 )
    {
        emit AlertMsg(QString("Error - No blinks are within the ROI"),'r');
        return;
    }

    QFile qdata;
    qdata.setFileName(OutputFile);
    if(!qdata.open(QFile::WriteOnly))
    {
        emit AlertMsg(QString("Cannot open %1 for writing").arg(OutputFile),'r');
        std::vector<InputBlinks>().swap(ROISTORM);
        return;
    }
    emit AlertMsg(QString("Writing segment..."), ' ');

    QTextStream pout(&qdata);
    for(uint j = 0; j < nROIBlinks; j++)
    {
         InputBlinks a = ROISTORM[j];
         pout << a.x << "  " << a.y << "  " << a.z << "  " << a.dX << "  " << a.dY << "  " << a.dZ << " " << a.Amplitude << " " <<  a.pic_no << "\n";
    }
    qdata.close();
    emit AlertMsg(QString("Wrote %1 points - Segment saved as %2").arg(nROIBlinks).arg(OutputFile),'b');
    std::vector<InputBlinks>().swap(ROISTORM);

    QString infoFile = OutputFile.remove(".3dSTM");
    infoFile += "_info.txt";

    QFile idata;
    idata.setFileName(infoFile);
    if(!idata.open(QFile::WriteOnly))
    {
       emit AlertMsg(QString("Cannot open %1 for writing").arg(infoFile),'r');
       return;
    }
    emit AlertMsg(QString("Done!"),'b');

    QTextStream ofi(&idata);
    ofi << "Original file : " << qPrintable(FileName) << "\n";
    ofi << "Segmented in RyR_dSTORM2D\n";
    ofi << "ROI (nm) : (" << ROIXStart << ", " << ROIYStart << "), (" << ROIXEnd << ", "  << ROIYEnd << ")\n";
    idata.close();
}
void STORMdensity::CalculateClusters(std::vector<int>& validvertex, std::vector<Cluster>& Clust)
{
   int nv = int(validvertex.size());
   int cno = 0;
   int totalpnts = 0;
   std::vector<Cluster>().swap(Clust);

   QString ClusterType = "Base Cluster";
   if(PassStage == Passtype::Second)
       ClusterType = "Cluster";

   for(int k=0; k < nv; k++)
   {
     int i  = validvertex[k];
     if (!edgedata[i].visitflag)
       continue;

     std::vector<int> clpnts;
     clpnts=tracepoints(i);
     std::sort(clpnts.begin(),clpnts.end());
     uint npntsb=uint(clpnts.size());
     auto last=std::unique(clpnts.begin(),clpnts.end());
     clpnts.erase(last,clpnts.end());
     uint n=uint(clpnts.size());
     std::vector<ClusterPoint> cps;
     uint nNoBlinks = 0;
     for(uint m=0; m < n; m++)
     {
       uint pnt = clpnts[m];
       ClusterPoint cp;
       cp.posn=bv[pnt].position();
       cp.density = bv[pnt].density();
       cp.frame = bv[pnt].frameNos();
       cp.amp = bv[pnt].amplitude();
       cp.logdensity = log10(cp.density);
       cp.NoBlinks = bv[pnt].weight();
       nNoBlinks += uint(cp.NoBlinks);
       cps.push_back(cp);
     }
     CGAL::spatial_sort( cps.begin(), cps.end(), ClusterPoint_sort_traits() );
     Cluster a;
     a.setClustNo(cno);
     a.addClustPoints(cps,PassStage);
     a.setNoBlinks(nNoBlinks);
     Clust.push_back(a);
     totalpnts += n;
     if (npntsb != n)
        emit AlertMsg(QString("%1 %2 before dups: %3 - after: %4 points").arg(ClusterType).arg(cno).arg(npntsb).arg(n),'r');
     else
        emit AlertMsg(QString("%1 %2 contains %3 points").arg(ClusterType).arg(cno).arg(n),'b'); // - max. depth = " << tracecount << " repeat visits  = " << visitedcalls,' ');
     cno++;
   }
   BlinkCoordinates = totalpnts;
   emit AlertMsg(QString("Total = %1 points").arg(totalpnts),'b');
}
std::vector<int> STORMdensity::tracepoints(int i)
{
    // This is called recursively
    std::vector<int> B;
    if(!edgedata[i].visitflag)
        return B;

    edgedata[i].visitflag = false;

    B.push_back(i);
    uint ne=uint(edgedata[i].endpnts.size());
    if(ne == 0)
        return B;

    std::vector<int> C;
    for(uint j=0;j < ne; j++)
    {
       int l=edgedata[i].endpnts[j];
        C=tracepoints(l);
       if (C.size() > 0)
           B.insert(B.end(), C.begin(), C.end());
    }

    return B;
}
