/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, cluster.h, is a header file for cluster.cpp and is part of the RyR_STORM2D program.
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
#ifndef CLUSTER_H
#define CLUSTER_H

#include <algorithm>
#include <string>
#include <tuple>
#include <QVector2D>
#include <QString>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Timer.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_rational.h>

typedef CGAL::Cartesian<double> R;
typedef R::Point_2 Point_2;

typedef  CGAL::Exact_rational             NT;
typedef  CGAL::Cartesian<NT>              L;
typedef  CGAL::Point_2<L>                 EPoint;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Segment_2<K> Segment;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Point_2<K> APoint;

struct AreaInfo
{
   bool isSimple;
   bool isConvex;
   bool isClockwise;
   double Area;
   double Size;
};

struct ClusterPoint
{
    Point_2 posn;
    double density;
    double logdensity;
    int NoBlinks;
};
struct ClusterPoint_less_x
{
    bool operator()(ClusterPoint p, ClusterPoint q) const
    { return (p.posn.x() < q.posn.x()); }
};
struct ClusterPoint_less_y
{
    bool operator()(ClusterPoint p, ClusterPoint q) const
    { return (p.posn.y() < q.posn.y()); }
};
struct ClusterPoint_sort_traits
{
    typedef ClusterPoint Point_2;
    typedef ClusterPoint_less_x Less_x_2;
    typedef ClusterPoint_less_y Less_y_2;

    Less_x_2 less_x_2_object() const
    { return Less_x_2(); }
    Less_y_2 less_y_2_object() const
    { return Less_y_2(); }
    ClusterPoint_sort_traits() {}
};

struct Boundary
{
    Polygon_2 Bounds;
    int nNoPoints;
    AreaInfo BoundInfo;
};
struct PolygonData
{
    Polygon_2 P;
    int nNoPnts;
    double Area;
    bool isContained;
};
struct Objects
{
    std::vector<Segment> Perim;
    std::vector<Boundary> Shapes;
    uint nNoObjects;
    double totalArea;
    bool isValid;
};

enum Passtype{First, Second};

class Cluster
{
    uint clustno;
    bool NotSimple;
    bool TooSmall;
    bool TooFewBlinks;
    bool BadShape;

    std::vector<Point_2> MinEllipsePoints;
    std::vector<Point_2> SupportPts;
    uint NumEllipsePoints;
    double EllipseSemiMajor;
    double EllipseSemiMinor;
    double EllipseEccentricity;

    uint NumSubClusters;
    std::vector<std::pair<uint, double>> BaseClusterInfo;

    std::vector<ClusterPoint> CPoints;
    uint nPoints;

    uint NumHullPoints;
    std::vector<Point_2> ConvexHullPoints;
    Polygon_2 ConvexHull;
    AreaInfo ConvexHullAreaInfo;
    bool ConvexHullValid;

    QString Status;
    Point_2 Centroid;
    Objects Outer;
    Point_2 lowestPoint;
    uint nNoBlinks;
    double logDensityMedian;
    double logDensityMean;
    double clusterArea;
    double e_to_edist;
    double nearestClusterArea;
    double twopi;

    void ellipseparameter(double ec[]);
    void AnalyseMultiplePolygons(std::vector<Boundary>& Boundaries, double& totalarea);
    Objects CalculateAlphaShape(std::vector<APoint> EPoints, double alphaval);
    AreaInfo CalculatePolygonArea(Polygon_2 p);
public:
    Cluster();
    ~Cluster();
    bool isWrongSize() {return (TooSmall || TooFewBlinks || BadShape);}
    bool isValid() { return (Outer.isValid && !isWrongSize());}
    bool isTooSmall() {return TooSmall;}
    void setTooSmall(bool v) {TooSmall = v;}
    void setHasTooFewBlinks(bool v) {TooFewBlinks = v;}
    bool hasTooFewBlinks() {return TooFewBlinks;}
    void setHasBadShape(bool v) {BadShape = v;}
    bool hasBadShape() {return BadShape;}
    bool isBoundaryValid() {return Outer.isValid;}
    void setBoundaryValidity(bool v) {Outer.isValid =v;}

    void setClustNo(uint i){clustno = i;}
    uint getClustNo() {return clustno;}
    uint getClustSize() {return nPoints;}
    void setNoBlinks(uint nB){nNoBlinks = nB;}
    uint getNoBlinks() {return nNoBlinks;}
    std::vector<Boundary> getBoundary() {return Outer.Shapes;}

    QString getStatus() {return Status;}
    void setStatus(QString s) {Status =s;}

    void addClustPoints(std::vector<ClusterPoint> b, int PassStage = Passtype::First);
    std::vector<ClusterPoint> getClustPoints() {return CPoints;}

    void addSubClusterInfo(std::vector<std::pair<uint, double>> B) {BaseClusterInfo = B; NumSubClusters = uint(BaseClusterInfo.size());}
    void addSubClusterInfo(uint cno, double cArea) {BaseClusterInfo.push_back(std::make_pair(cno, cArea)); NumSubClusters = uint(BaseClusterInfo.size());}
    std::vector<std::pair<uint, double>> getBaseClusterInfo() {return BaseClusterInfo;}
    uint getNoSubClusters() {return uint(NumSubClusters);}

    void CalculateLowestPoint();
    Point_2 getLowestPoint() {return lowestPoint;}

    void CalculateDensityMeanMedian();
    std::pair<double,double> getDensityMeanMedian(){return (std::make_pair(logDensityMean, logDensityMedian));}

    void setEdgetoEdgeDist(double e){e_to_edist=e;}
    void setNearestClusterArea(double a ){nearestClusterArea = a;}
    double getNND() {return e_to_edist;}
    double getNearestArea() {return nearestClusterArea;}

    void CalculateConvexHull();
    std::vector<Point_2> getConvexHullPoints() {return ConvexHullPoints;}
    uint getNumHullPoints() {return NumHullPoints;}
    Polygon_2 getConvexHull() {return ConvexHull;}
    double getConvexHullArea() {return ConvexHullAreaInfo.Area;}
    bool isConvexHullValid() {return ConvexHullValid;}

    void CalculateMinEllipse();
    std::vector<Point_2> getEllipsePoints() {return MinEllipsePoints;}
    uint getNumEllipsePoints() {return NumEllipsePoints;}
    std::tuple<double, double> getEllipseData() {return std::make_tuple(EllipseSemiMajor, EllipseSemiMinor);}

    void CalculateCentroid();
    Point_2 getCentroid() {return Centroid;}

    double getArea() {return Outer.totalArea;}
    void setArea(double val) {Outer.totalArea = val;}
    Polygon_2 getOuterPolygon() {return Outer.Shapes[0].Bounds; }
    uint getNumBoundarySegmentPts() {return(uint(2*Outer.Perim.size()));}
    uint getNumOuterObjects() {return Outer.nNoObjects;}

    std::vector<QVector2D> getBoundarySegmentPts();

    bool CalculateBoundaries(double alphaval);

    static double MinTetramerArea;
    static double MinBlinksPerCluster;
//signals:
//    void AlertMsg(QString msg, char Color);
};

#endif // CLUSTER_H
