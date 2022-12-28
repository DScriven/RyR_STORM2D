/*********************************************************************************************
* RyR_STORM2D is a program designed to analyse blink data from the surface of ventricular
* myocytes.
*
* Copyright David Scriven, 2022.
*
* Moore Laboratory, Life Sciences Institute,  2350 Health Sciences
* Mall, University of British Columbia, Vancouver, Canada, V6T 1Z3
*
* This file, cluster.cpp, is part of the RyR_STORM2D program and characterises and stores information
* about the blink clusters
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

#include "stormdensity.h"
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/enum.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <algorithm>
#include <tuple>

typedef  CGAL::Min_ellipse_2_traits_2<L>  ETraits;
typedef  CGAL::Min_ellipse_2<ETraits>     Ellipse_2;

typedef K::FT FT;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vba;
typedef CGAL::Alpha_shape_face_base_2<K>  Fba;
typedef CGAL::Triangulation_data_structure_2<Vba,Fba> Tdsa;
typedef CGAL::Delaunay_triangulation_2<K,Tdsa> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;
typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
typedef Polygon_2::Vertex_iterator VertexIterator;

double Cluster::MinTetramerArea;
double Cluster::MinBlinksPerCluster;

bool compareArea(Boundary b1, Boundary b2)
{
    return (b1.BoundInfo.Area > b2.BoundInfo.Area);
}
Cluster::Cluster()
{
    TooSmall = false;
    TooFewBlinks = false;
    BadShape = false;
    NumHullPoints = 0;
    Outer.isValid = false;
    NotSimple = false;
    Status = "unknown";
    twopi=8.0*atan(1.0);
}
Cluster::~Cluster()
{
   std::vector<ClusterPoint>().swap(CPoints);
   std::vector<Segment>().swap(Outer.Perim);
   std::vector<Boundary>().swap(Outer.Shapes);
   std::vector<Point_2>().swap(MinEllipsePoints);
   std::vector<Point_2>().swap(ConvexHullPoints);
}
void Cluster::addClustPoints(std::vector<ClusterPoint> b, int PassStage)
{
    CPoints = b;
    nPoints=uint(CPoints.size());
    TooFewBlinks = false;
    if (nPoints < MinBlinksPerCluster)
    {
       Status = QString("has %1 blinks &lt; %2 required").arg(nPoints).arg(MinBlinksPerCluster);
       TooFewBlinks=true;
       NumHullPoints = 0;
       setArea(0.0);
       return;
    }
    if (PassStage == Passtype::Second)
    {
      CalculateLowestPoint();
      CalculateCentroid();
      CalculateConvexHull();
      CalculateDensityMeanMedian();
      CalculateMinEllipse();
    }
}
void Cluster::CalculateLowestPoint()
{
    Point_2 a = CPoints[0].posn;
    double xpos = a.x();
    double ypos = a.y();
    for(uint k=1; k < nPoints; k++)
    {
        Point_2 a1 = CPoints[k].posn;
        xpos += a1.x();
        ypos = (ypos < a1.y()) ? ypos : a1.y();
    }
    lowestPoint = Point_2(xpos/nPoints, ypos);
}
void Cluster::CalculateDensityMeanMedian()
{
    std::vector <double> ldens;
    double logdenssum=0;
    for(uint k=0; k < nPoints; k++)
    {
        double logdens = CPoints[k].logdensity;
        logdenssum += logdens;
        ldens.push_back(logdens);
    }
    logDensityMean = logdenssum/nPoints;
    std::sort (ldens.begin(), ldens.end());
    uint midpnt = nPoints/2;
    if (nPoints%2 == 0)
        logDensityMedian = (ldens[midpnt-1] + ldens[midpnt])/2;
    else
        logDensityMedian = ldens[midpnt];
}
void Cluster::CalculateConvexHull()
{
    ConvexHullPoints.clear();
    ConvexHull.clear();
    if(hasTooFewBlinks() || hasBadShape())
       return;
    std::vector <APoint> HPoints;
    std::vector <APoint> CHPoints;
    for(uint j = 0; j < nPoints; j++)
    {
        APoint plc = APoint(CPoints[j].posn.x(), CPoints[j].posn.y());
        HPoints.push_back(plc);
    }
    CGAL::convex_hull_2( HPoints.begin(), HPoints.end(), std::back_inserter(CHPoints));
// Convert from APoint to QVector2D for plotting
    for(uint i=0; i < CHPoints.size(); i++)
    {
       Point_2 ch = Point_2(CHPoints[i].x(),CHPoints[i].y());
       ConvexHull.push_back(CHPoints[i]);
       ConvexHullPoints.push_back(ch);
    }
    NumHullPoints = uint(ConvexHullPoints.size());
    ConvexHullAreaInfo=CalculatePolygonArea(ConvexHull);
    if(ConvexHullAreaInfo.isSimple && ConvexHullAreaInfo.isConvex)
       ConvexHullValid = true;
    else
       ConvexHullValid = false;
}
void Cluster::CalculateMinEllipse()
{
    TooFewBlinks = false;
    uint noPnts = uint(CPoints.size());
    std::vector<EPoint> EPoints;
    if (NumHullPoints == 0)
    {
        for(uint j = 0; j < noPnts; j++)
        {
            EPoint plc = EPoint(CPoints[j].posn.x(), CPoints[j].posn.y());
            EPoints.push_back(plc);
        }
    }
    else
    {
        for(uint j = 0; j < NumHullPoints; j++)
        {
            EPoint plc = EPoint(ConvexHullPoints[j].x(), ConvexHullPoints[j].y());
            EPoints.push_back(plc);
        }
    }
    Ellipse_2 MinEllipse(EPoints.begin(),EPoints.end(),true);
    if(MinEllipse.is_degenerate())
    {
       Status = QString("has a degenerate ellipse");
       BadShape = true;
       NumEllipsePoints = 0;
       setArea(0.0);
       return;
    }
    double ec[6];
    MinEllipse.ellipse().double_coefficients( ec[0], ec[2], ec[1], ec[3], ec[4], ec[5]);

    double minval = 1.e60;
    for(int j=0; j < 6; j++)
       minval = (abs(ec[j]) < minval) ? abs(ec[j]) : minval;

    for(int j=0; j < 6; j++)
       ec[j] /= minval;

    ellipseparameter(ec);
    double halftetramer = 12.25;
    if(EllipseSemiMinor < halftetramer)
    {
        Status = QString("has ellipse semi-minor axis = %1 &lt; %2 nm").arg(EllipseSemiMinor).arg(halftetramer);
        BadShape=true;
        MinEllipsePoints.clear();
        NumEllipsePoints = 0;
    }
    return;
}
void Cluster::ellipseparameter(double ec[])
{
    /* Conversion of Ellipse general equation to a trigonometric parameterization and calculation  of
     * major and minor axis from Wikipedia ellipse found at
     * https://www.wikipedia.org/wiki/Ellipse#General_ellipse */

    double phi = ec[1]*ec[1] - 4.0*ec[0]*ec[2];
    double det = 8.0*ec[0]*ec[2]*ec[5] + 2.0*ec[1]*ec[3]*ec[4] - 2.0*ec[1]*ec[1]*ec[5] - 2.0*ec[2]*ec[3]*ec[3]- 2.0*ec[0]*ec[4]*ec[4];
    double sigma = ec[0] + ec[2];
    double h = (2.0*ec[2]*ec[3] - ec[1]*ec[4])/phi;
    double k = (2.0*ec[0]*ec[4] - ec[1]*ec[3])/phi;
    double amc = ec[0] - ec[2];
    double rtamcec1 = sqrt(amc*amc + ec[1]*ec[1]);
    double theta;
    if(abs(ec[1] - 0.0) < DBL_EPSILON)
    {
       if(ec[0] < ec[2])
         theta = 0.0;
       else
         theta = twopi/4.0;
    }
    else
       theta=atan((-amc-rtamcec1)/ec[1]);

    double costheta = cos(theta);
    double sintheta = sin(theta);
    double alpha = -sqrt(-det*(sigma+rtamcec1))/phi;
    double beta = -sqrt(-det*(sigma-rtamcec1))/phi;
    double ratio;
    if (alpha > beta)
    {
       EllipseSemiMajor = alpha;
       EllipseSemiMinor = beta;
       ratio = beta*beta/(alpha*alpha);
    }
    else
    {
       EllipseSemiMajor = beta;
       EllipseSemiMinor = alpha;
       ratio = alpha*alpha/(beta*beta);
    }
    EllipseEccentricity = sqrt(1.0-(ratio));

    uint NumPoints = 120;
    MinEllipsePoints.clear();
    double delta = twopi/NumPoints;
    double ac = alpha*costheta;
    double bs = beta*sintheta;
    double as = alpha*sintheta;
    double bc = beta*costheta;
    for(uint i=0; i < NumPoints; i++)
    {
       double t = i*delta;
       double x = ac*cos(t) - bs*sin(t) + h;
       double y = as*cos(t) + bc*sin(t) + k;
       Point_2 ep = Point_2(x,y);
       MinEllipsePoints.push_back(ep);
    }
    NumEllipsePoints=uint(MinEllipsePoints.size());
}
bool Cluster::CalculateBoundaries(double alphaval)
{
    uint noPnts = uint(CPoints.size());
    if (hasTooFewBlinks() || hasBadShape())
        return false;

    std::vector <APoint> OPoints;
    for(uint j = 0; j < noPnts; j++)
    {
        APoint plc = APoint(CPoints[j].posn.x(), CPoints[j].posn.y());
        OPoints.push_back(plc);
    }
    Outer = CalculateAlphaShape(OPoints, alphaval);
    if(!Outer.isValid)
    {
        double xalpha = alphaval;
        std::vector<int> v = {2,3,4,5,6,7,8,9,10,20,50};
        for(auto &k : v)
        {
          xalpha = k*alphaval;
          Outer = CalculateAlphaShape(OPoints, xalpha);
          if(Outer.isValid)
             return true;
        }
        Outer.isValid = false;
        if(NotSimple)
           setStatus(QString("Outer Boundary - not simple - alpha = %1").arg(xalpha));
        else
           setStatus(QString("Outer Boundary - cannot draw alpha shape - alpha = %1").arg(xalpha));
        return false;
    }
    else
    {
        if(Outer.totalArea < MinTetramerArea )
        {
            Status=QString("Area = %1 nm<sup>2</sup> &lt; %2 nm<sup>2</sup>").arg(Outer.totalArea).arg(MinTetramerArea);
            TooSmall = true;
            NumHullPoints = 0;
            return false;
        }
    }

    return true;
}
Objects Cluster::CalculateAlphaShape(std::vector<APoint> EPoints, double alphaval)
{
    Objects alpha;
    alpha.isValid = false;
    alpha.Perim.clear();
    alpha.Shapes.clear();
    NotSimple = false;

    std::vector<Segment> segments;
    Alpha_shape_2 A(EPoints.begin(), EPoints.end(), FT(alphaval), Alpha_shape_2::GENERAL);
    for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
        it != A.alpha_shape_edges_end(); ++it){
        segments.push_back(A.segment(*it));
    }
    uint alphasegsize = uint(segments.size());

    if(alphasegsize == 0)
    {
        std::cerr << "Cluster " << getClustNo() << " Catastrophic error - alpha shape has no segments\n";
        return alpha;
    }

    for(uint kz=0; kz < alphasegsize; kz++)
        alpha.Perim.push_back(segments[kz]);

    std::vector<APoint> SegPnts;
    std::vector <bool> SegmentUsed;
    std::vector <APoint> SourceMap;
    for (uint k = 0; k < alphasegsize; k++)
    {
        SourceMap.push_back(segments[k].source());
        SegmentUsed.push_back(false);
    }

    SegmentUsed[0] = true;

    int NoObjects=0;
    int FoundVertex = -1;
    APoint startvertex = segments[0].source();
    APoint lastvertex = segments[0].target();
    SegPnts.push_back(startvertex);
    SegPnts.push_back(lastvertex);
    uint NoMatched = 1;
    uint VerticesPerObject = 2;
    bool bClosedAnObject = false;
    double totalarea = 0;
    std::vector<Boundary> Boundaries;
    while (NoMatched < alphasegsize)
    {
        std::vector<APoint>::iterator it2 = std::find(SourceMap.begin(), SourceMap.end(), lastvertex);
        if(it2 == SourceMap.end())
        {
            if(VerticesPerObject > 4)
            {
                SegPnts.push_back(startvertex); // reconnect to begin
                lastvertex = startvertex;
                bClosedAnObject = true;
            }
            else
            {
                SegPnts.clear();
                lastvertex = startvertex;
                bClosedAnObject = false;
            }
        }
        else
        {
            FoundVertex = std::distance(SourceMap.begin(), it2);
            if(segments[FoundVertex].source() != lastvertex)
            {
                Status = "error in segment match";
                return alpha;
            }
            //we're good - found the next segment
            SegmentUsed[FoundVertex] = true;
            NoMatched++;
            VerticesPerObject++;
            lastvertex = segments[FoundVertex].target();
            if(startvertex == lastvertex)
                bClosedAnObject = true;
        }

        if (startvertex == lastvertex)
        {
            if(bClosedAnObject) // we have a closed loop and thus an object
            {
                Boundary b;
                for(uint k = 0; k < SegPnts.size(); k++)
                    b.Bounds.push_back(SegPnts[k]);

                b.nNoPoints = uint(SegPnts.size());
                if(b.nNoPoints > 4)
                {
                  b.BoundInfo = CalculatePolygonArea(b.Bounds);
                  if(b.BoundInfo.isSimple)
                  {
                    Boundaries.push_back(b);
                    NoObjects++;
                  }
                  else
                    NotSimple=true;
                }
                SegPnts.clear();
            }
            uint nSegsLeft = alphasegsize - NoMatched;
            if (nSegsLeft > 3) // there's another object
            {
                std::vector<bool>::iterator it3 = std::find(SegmentUsed.begin(), SegmentUsed.end(),false);
                if(it3 == SegmentUsed.end())
                {
                    if (NoObjects > 0)
                    {
                        alpha.nNoObjects = NoObjects;
                        alpha.totalArea = totalarea;
                        alpha.isValid = true;
                        Status = "objects omitted - could not find unused segment in source map";
                    }
                    else
                    {
                        alpha.isValid = false;
                        Status = "no objects saved - could not find segment in source map";
                    }
                    return alpha;
                }
                uint NewVertexPosn = std::distance(SegmentUsed.begin(), it3);
                SegmentUsed[NewVertexPosn] = true;
                startvertex = segments[NewVertexPosn].source();
                lastvertex = segments[NewVertexPosn].target();
                SegPnts.push_back(startvertex);
                SegPnts.push_back(lastvertex);
                NoMatched++;
                VerticesPerObject = 2;
            }  // end - another object
        } // end - we have an object
        else
            SegPnts.push_back(lastvertex);
    } // end - while looking at all the segments

    uint segsleftsize=uint(SegPnts.size());

    if(Boundaries.size() == 0 || NoObjects == 0 || NotSimple)
    {
        if(NotSimple)
        {
            Status = QString(" Outer Boundary is not simple");
            alpha.nNoObjects = NoObjects;
        }
        else
        {
            Status = QString(" cannot calculate exterior boundary with %1 segments").arg(segsleftsize);
            alpha.nNoObjects = 0;
        }
        alpha.Shapes = Boundaries;
        alpha.totalArea = 0.0;
        alpha.isValid=false;
        return alpha;
    }
    else
    {
        if(NoObjects == 1 || segsleftsize == 0)
        {
            totalarea = Boundaries[0].BoundInfo.Area;
            std::cerr << "Cluster # " << getClustNo() << " -  External Boundary has " << Boundaries[0].nNoPoints << " points and area " << totalarea << " nm2\n";
        }
        else
            AnalyseMultiplePolygons(Boundaries, totalarea);
        alpha.isValid = true;
    }

    alpha.Shapes = Boundaries;
    alpha.nNoObjects = NoObjects;
    alpha.totalArea = totalarea;
    return alpha;
}
void Cluster::AnalyseMultiplePolygons(std::vector<Boundary>& Boundaries, double& totalarea)
{
    std::sort(Boundaries.begin(), Boundaries.end(), compareArea);
    int NoPolygons = int(Boundaries.size());
    std::vector<PolygonData> PolyData;

    Polygon_2 A = Boundaries[0].Bounds;

    CGAL::Orientation FirstOrientation = CGAL::CLOCKWISE;
    if(A.orientation() != FirstOrientation)
       FirstOrientation = CGAL::COUNTERCLOCKWISE;

    PolygonData PD;
    PD.P = A;
    PD.Area = Boundaries[0].BoundInfo.Area;
    PD.nNoPnts = Boundaries[0].nNoPoints;
    PD.isContained = false;
    PolyData.push_back(PD);
// make sure polygons are simple and are oriented the same as the largest
// so that we can determine whether they intersect

    for(int id=1; id < NoPolygons; id++)
    {
        Polygon_2 B = Boundaries[id].Bounds;
        if(!B.is_simple() || B.is_collinear_oriented())
           continue;
        if (B.orientation() != FirstOrientation)
        {
            B.reverse_orientation();
        }
        PolygonData PD;
        PD.P = B;
        PD.Area = Boundaries[id].BoundInfo.Area;
        PD.nNoPnts = Boundaries[id].nNoPoints;
        PD.isContained = false;
        PolyData.push_back(PD);
    }
    if(PolyData.size() < 2)
    {
        totalarea = PolyData[0].Area;
        return;
    }


    for(int id = int(PolyData.size()-1); id >= 0 ; id--)
    {
        Polygon_2 A = PolyData[id].P;
        if (PolyData[id].isContained)  // Assume polygon can't be contained in more than one Polygon
           continue;
        for(int ie = id-1; ie >=0 ; ie--)
        {
            Polygon_2 B = PolyData[ie].P;
            if (do_intersect(A,B))
            {
                PolygonData T = PolyData[id];
                T.Area = -T.Area;
                T.isContained = true;
                PolyData[id] = T;
                break;
            }
        }
    }

    uint bint = 0;
    uint bext = 0;
    for(uint id = 0; id < PolyData.size(); id++)
    {
        if(!PolyData[id].isContained)
            std::cerr << "Cluster # " << getClustNo() << " -  External Boundary " << (++bext) << " has " << PolyData[id].nNoPnts << " points and area " << PolyData[id].Area << " nm2\n";
        else
            std::cerr << "Cluster # " << getClustNo() << " -  Internal Boundary (Hole) " << (++bint) << " has " << PolyData[id].nNoPnts << " points and excludes area " << (-PolyData[id].Area) << " nm2\n";
        totalarea += PolyData[id].Area;
    }
}
AreaInfo Cluster::CalculatePolygonArea(Polygon_2 p)
{
    AreaInfo ai;
    ai.Area = 0.0;
    ai.isSimple = p.is_simple();
    if(!ai.isSimple)
       return ai;
    ai.isConvex = p.is_convex();
    ai.isClockwise = (p.orientation() == CGAL::CLOCKWISE);
    ai.Area = abs(p.area());

    return ai;
}
std::vector<QVector2D> Cluster::getBoundarySegmentPts()
{
    std::vector<QVector2D> a;
    if(!isBoundaryValid())
       return a;
    int segsize=int(Outer.Perim.size());
    for (int ia=0; ia < segsize; ia++)
    {
        Segment b = Outer.Perim[ia];
        a.push_back(QVector2D(b.source().x(),b.source().y()));
        a.push_back(QVector2D(b.target().x(),b.target().y()));
    }
    return a;
}
void Cluster::CalculateCentroid()
{
    double xclust = 0;
    double yclust = 0;
    int nBlinks = 0;
    for(uint i=0; i < nPoints; i++)
    {
        Point_2 p = CPoints[i].posn;
        int n=CPoints[i].NoBlinks;
        nBlinks += n;
        xclust += n*p.x();
        yclust += n*p.y();
    }
    int xp = qRound(xclust/nBlinks);
    int yp = qRound(yclust/nBlinks);
    Centroid = Point_2(xp, yp);
}
